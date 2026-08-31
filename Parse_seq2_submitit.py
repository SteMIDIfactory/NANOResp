import sys, os, time, glob, socket
import json
import subprocess
import submitit
from datetime import datetime
from Bio import SeqIO
from collections import defaultdict
import argparse

## HELPER FUNCTIONS
def get_startingpoint(input_folder):
    """
    Scans the input folder to find the earliest FASTQ file timestamp.
    Uses 'ontime' as the primary method, falling back to OS modification time.
    """
    inputs = glob.glob(os.path.join(input_folder, "**", "*.fastq*"), recursive=True)
    earliest = 0
    for i in inputs:
        try:
            out = subprocess.check_output(["ontime", "-s", i], text=True).strip()
            start = out.split("\n")[0].replace("Earliest:", "").split(".")[0]
            start = datetime.strptime(start, " %Y-%m-%dT%H:%M:%S")
            if earliest == 0 or start < earliest:
                earliest = start
        except Exception:
            continue
    if earliest == 0:
        if inputs:
            return min(os.path.getmtime(f) for f in inputs)
        else:
            raise RuntimeError("Unable to determine starting point: no FASTQ files found.")
    return earliest.timestamp() if isinstance(earliest, datetime) else earliest

def extract_data(input_folder, tt, all_start, tmp_dir):
    """
    Used EXCLUSIVELY in ONTIME mode (static files).
    Extracts reads generated up to the specified timepoint (tt) using 'ontime'.
    """
    inputs = glob.glob("%s/*.fastq" % (input_folder))
    start_dt = datetime.fromtimestamp(all_start) if isinstance(all_start, (int, float)) else all_start
    start_str = start_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    for g in inputs:
        out_name = os.path.basename(g)
        cmd = f"ontime -o {tmp_dir}/{out_name} -f {start_str} -t {tt}m {g}"
        os.system(cmd)

def concatenate_live_barcodes(input_dir, output_dir):
    """
    Used EXCLUSIVELY in LIVE mode:
    Scans subdirectories (1 subdir = 1 sample/barcode) and concatenates
    all fastq/fastq.gz files accumulated up to that moment into a single FASTQ.
    """
    os.makedirs(output_dir, exist_ok=True)
    subdirs = [d for d in glob.glob(os.path.join(input_dir, "*")) if os.path.isdir(d)]

    concat_files = []
    for s_dir in subdirs:
        sample_name = os.path.basename(s_dir)
        fastqs = glob.glob(os.path.join(s_dir, "*.fastq*"))
        if not fastqs:
            continue

        out_fq = os.path.join(output_dir, f"{sample_name}.fastq")
        cmd = f"zcat -f {' '.join(fastqs)} > {out_fq}"
        subprocess.run(cmd, shell=True)
        concat_files.append(out_fq)

    return concat_files

def get_best_mlst_scheme(global_kma_res):
    """Parses KMA results to identify the most likely MLST scheme based on alignment scores."""
    if not os.path.exists(global_kma_res):
        return None

    scheme_scores = {}
    with open(global_kma_res, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                template_name = parts[0]
                if "___" in template_name:
                    scheme = template_name.split("___")[0]
                    try:
                        coverage = float(parts[5])
                        identity = float(parts[4])
                        if coverage >= 50.0 and identity >= 70.0:
                            scheme_scores[scheme] = scheme_scores.get(scheme, 0.0) + float(parts[3])
                    except ValueError:
                        continue

    return max(scheme_scores, key=scheme_scores.get) if scheme_scores else None

def run_gene_consensus_mlst(fastq_file, best_scheme, threads, local_scratch, shared_dir):
    """Maps reads to the chosen MLST scheme, builds a consensus, and blasts to identify alleles."""
    schemes_base = os.path.join(shared_dir, "KMA_DB/MLST_system/schemes")
    scheme_dir = os.path.join(schemes_base, best_scheme)
    scheme_fasta = os.path.join(scheme_dir, "all_alleles.fasta")

    if not os.path.exists(scheme_fasta):
        return "NA"

    paf_file = os.path.join(local_scratch, f"mlst_{best_scheme}.paf")
    os.system(f"minimap2 -x map-ont -t {threads} {scheme_fasta} {fastq_file} > {paf_file} 2>/dev/null")

    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return "NA"

    gene_to_reads = defaultdict(set)
    gene_allele_counts = defaultdict(lambda: defaultdict(int))

    with open(paf_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                read_id = parts[0]
                allele_id = parts[5]
                gene_name = allele_id.rsplit("_", 1)[0] if "_" in allele_id else allele_id

                gene_to_reads[gene_name].add(read_id)
                gene_allele_counts[gene_name][allele_id] += 1

    if not gene_to_reads:
        return "NA"

    fastq_dict = SeqIO.to_dict(SeqIO.parse(fastq_file, "fastq"))
    ref_dict = SeqIO.to_dict(SeqIO.parse(scheme_fasta, "fasta"))

    blast_db = os.path.join(local_scratch, f"blast_db_{best_scheme}")
    os.system(f"makeblastdb -in {scheme_fasta} -dbtype nucl -out {blast_db} 2>/dev/null")

    final_alleles = []

    for gene_name, reads_set in sorted(gene_to_reads.items()):
        if len(reads_set) < 3:
            continue

        gene_scratch = os.path.join(local_scratch, f"gene_{gene_name}")
        os.makedirs(gene_scratch, exist_ok=True)

        gene_fastq = os.path.join(gene_scratch, "reads.fastq")
        with open(gene_fastq, "w") as out_fq:
            for r_id in reads_set:
                if r_id in fastq_dict:
                    SeqIO.write(fastq_dict[r_id], out_fq, "fastq")

        best_guide_allele = max(gene_allele_counts[gene_name], key=gene_allele_counts[gene_name].get)
        guide_fasta = os.path.join(gene_scratch, "guide.fasta")
        if best_guide_allele in ref_dict:
            with open(guide_fasta, "w") as out_g:
                SeqIO.write(ref_dict[best_guide_allele], out_g, "fasta")
        else:
            continue

        bam_file = os.path.join(gene_scratch, "aligned.bam")
        consensus_fasta = os.path.join(gene_scratch, "consensus.fasta")
        cmd_consensus = (
            f"minimap2 -ax map-ont -t {threads} {guide_fasta} {gene_fastq} 2>/dev/null | "
            f"samtools view -bS - | samtools sort -o {bam_file} - && "
            f"samtools index {bam_file} 2>/dev/null && "
            f"samtools consensus -f fasta {bam_file} -o {consensus_fasta} 2>/dev/null"
        )
        os.system(cmd_consensus)

        if os.path.exists(consensus_fasta) and os.path.getsize(consensus_fasta) > 0:
            blast_out = os.path.join(gene_scratch, "blast.tsv")
            os.system(
                f"blastn -query {consensus_fasta} -db {blast_db} "
                f"-outfmt '6 sseqid pident length slen' -perc_identity 80 -out {blast_out} 2>/dev/null"
            )

            best_match = "NA"
            if os.path.exists(blast_out):
                with open(blast_out, "r") as bf:
                    for b_line in bf:
                        b_parts = b_line.strip().split("\t")
                        if len(b_parts) >= 4:
                            hit_allele = b_parts[0]
                            pident = float(b_parts[1])
                            length = int(b_parts[2])
                            slen = int(b_parts[3])

                            hit_gene = hit_allele.rsplit("_", 1)[0] if "_" in hit_allele else hit_allele
                            if hit_gene == gene_name:
                                if pident == 100.0 and length >= int(slen * 0.99):
                                    best_match = hit_allele
                                    break
                                elif best_match == "NA":
                                    best_match = f"{hit_allele}({pident:.1f}%)"

            if best_match != "NA":
                final_alleles.append(best_match)

    return ";".join(sorted(final_alleles)) if final_alleles else "NA"

def run_two_step_kma_mlst(fastq_file, threads, local_scratch, shared_dir):
    """Executes the two-step MLST profiling using KMA."""
    global_db = os.path.join(shared_dir, "KMA_DB/MLST_system/global/pubmlst_global")
    prefix_global = os.path.join(local_scratch, "kma_mlst_global")
    os.system(f"kma -i {fastq_file} -t_db {global_db} -o {prefix_global} -t {threads} -1t1 2>{local_scratch}/kma_mlst_global.log")

    best_scheme = get_best_mlst_scheme(f"{prefix_global}.res")
    if not best_scheme:
        return {"detected_scheme": "NA", "alleles": "NA"}

    alleles = run_gene_consensus_mlst(fastq_file, best_scheme, threads, local_scratch, shared_dir)
    return {"detected_scheme": best_scheme, "alleles": alleles}

## WORKER MACRO-FUNCTION (Submitit Job)
def process_single_sample(sample_path, threads, tt, manager_node, assemblies_dir, shared_dir, env_config):
    """
    Main analysis pipeline for a single sample. 
    Runs on worker nodes via Slurm.
    """
    worker_node = socket.gethostname()
    sample_name = os.path.basename(sample_path)

    local_scratch = os.path.join("/scratch/sgaiarsa", f"worker_T{tt}_{sample_name}")
    os.makedirs(local_scratch, exist_ok=True)
    local_fastq = os.path.join(local_scratch, sample_name)

    if worker_node != manager_node:
        os.system(f"scp -q {manager_node}:{sample_path} {local_fastq}")
    else:
        os.system(f"cp {sample_path} {local_fastq}")

    res = {}

    try:
        if not os.path.exists(local_fastq) or os.path.getsize(local_fastq) == 0:
            return {sample_name: {"num_reads": 0, "status": "Empty or Missing FASTQ"}}

        line_count = os.popen(f"wc -l {local_fastq}").read().strip().split()
        if not line_count or int(line_count[0]) == 0:
            return {sample_name: {"num_reads": 0, "status": "Empty FASTQ"}}

        checknum = int(line_count[0]) / 4

        # 1. Porechop
        porechop_bin = os.path.join(shared_dir, "Porechop/porechop-runner.py")
        cut_file = f"{local_fastq}_cut"
        os.system(f"python {porechop_bin} -i {local_fastq} -o {cut_file} --threads {threads}")

        if not os.path.exists(cut_file) or os.path.getsize(cut_file) == 0:
            return {sample_name: {"num_reads": 0, "status": "Porechop failed or empty output"}}

        dict_seqs = {}
        cut2_file = f"{local_fastq}_cut2"
        with open(cut2_file, "w") as outf:
            for seq in SeqIO.parse(cut_file, "fastq"):
                dict_seqs[seq.id] = dict_seqs.get(seq.id, 0) + 1
                seq.id = f"{seq.id}_{dict_seqs[seq.id]}"
                SeqIO.write(seq, outf, "fastq")

        os.rename(cut2_file, cut_file)
        res["num_reads"] = int(os.popen(f"wc -l {cut_file}").read().strip().split()[0]) / 4

        # 2. Host Filtering (minimap2)
        hg38_ref = os.path.join(shared_dir, "Hg38.fna")
        maplog = f"{local_fastq}.maplog"
        map_out = []

        for _ in range(10):
            map_out = os.popen(f"minimap2 -t {threads} -ax map-ont {hg38_ref} -K 200M {cut_file} 2>{maplog} | grep -v '^@'").read().strip("\n").split("\n")
            numap = os.popen(f"grep '] mapped ' {maplog}").read().strip("\n").split("\n")
            if numap and numap[0]:
                try:
                    mapped = int(numap[0].split("]")[1].strip().split()[1])
                    if mapped >= float(checknum) * 0.9 and mapped <= float(checknum) * 1.1:
                        break
                except Exception:
                    continue

        chopselected = f"{local_fastq}_chopselected.fastq"
        included, bases = [], 0
        with open(chopselected, "w") as oF:
            for m in map_out:
                if m and m[0] != "@":
                    m_split = m.split("\t")
                    if len(m_split) >= 11 and m_split[2] == "*" and m_split[0] not in included:
                        oF.write(f"@{m_split[0]}\n{m_split[9]}\n+\n{m_split[10]}\n")
                        included.append(m_split[0])
                        bases += len(m_split[9])

        res["nohuman_bases"] = bases
        res["num_nohuman_reads"] = len(included)

        if len(included) == 0:
            return {sample_name: res}

        # 3. Kraken2
        kraken_db = os.path.join(shared_dir, "krakenDB/standard16")
        kraken_out = f"{local_fastq}_kraken"
        kraken_report = f"{local_fastq}.report"

        for _ in range(10):
            os.system(f"kraken2 --db {kraken_db} --threads {threads} --output {kraken_out} --memory-mapping --use-names --report {kraken_report} --use-mpa-style {chopselected}")
            if os.path.exists(kraken_report) and os.stat(kraken_report).st_size != 0:
                break

        taxlist, taxhash = {}, {}
        if os.path.exists(kraken_report):
            with open(kraken_report, "r") as f:
                for s in f:
                    parts = s.split("\t")
                    if len(parts) > 1 and parts[0].split("|")[-1][0] == "s":
                        taxlist[parts[0].split("|")[-1].split("_")[-1]] = int(parts[1])

        if os.path.exists(kraken_out):
            tax = os.popen(f"cut -f3 {kraken_out} | sort | uniq").read().strip("\n").split("\n")
            for t in tax:
                t_split = t.split(" (taxid ")
                if len(t_split) > 1 and t_split[0] in taxlist:
                    taxhash[t_split[0]] = t_split[1].split(")")[0]

        res["taxa"] = taxlist
        if not taxlist:
            return {sample_name: res}

        main_tax = max(taxlist, key=taxlist.get)
        res["main"] = main_tax

        kraken_tools = os.path.join(shared_dir, "KrakenTools/extract_kraken_reads.py")
        mlst_taxa, assembly_taxa = [], []

        for ttax, count in taxlist.items():
            if count >= 500:
                mlst_taxa.append(ttax)
                if count >= 1000:
                    assembly_taxa.append(ttax)
                reads_out = f"{local_fastq}_{ttax.replace(' ', '_')}_reads.fastq"
                os.system(f"python {kraken_tools} -k {kraken_out} -s {chopselected} -o {reads_out} -t {taxhash[ttax]} -r {kraken_report} --fastq-output --include-children")

        res["mlst"] = mlst_taxa if mlst_taxa else "NA"
        res["assembly"] = assembly_taxa if assembly_taxa else "NA"

        if not assembly_taxa and main_tax in taxhash:
            main_reads_out = f"{local_fastq}_{main_tax.replace(' ', '_')}_reads.fastq"
            os.system(f"python {kraken_tools} -k {kraken_out} -s {chopselected} -o {main_reads_out} -t {taxhash[main_tax]} -r {kraken_report} --fastq-output --include-children")

        # 4. Abricate AMR (CARD) on reads
        chopselected_fasta = f"{chopselected}.fasta"
        os.system(f"seqtk seq -a {chopselected} > {chopselected_fasta}")

        genes_overall = os.popen(f"conda run -n {env_config['abricate']} abricate --db card --threads {threads} {chopselected_fasta} | cut -f6,10,11 | tail -n +2").read().strip("\n").split("\n")
        genepresence_overall = []
        if genes_overall == ['']:
            res["AMR_overall"] = "NA"
        else:
            for g in genes_overall:
                g_parts = g.strip().split()
                if len(g_parts) >= 3 and float(g_parts[1]) >= 80 and float(g_parts[2]) >= 95:
                    genepresence_overall.append(g_parts[0])
            res["AMR_overall"] = ";".join(list(set(genepresence_overall))) if genepresence_overall else "NA"

        prevalent_fastq = f"{local_fastq}_{main_tax.replace(' ', '_')}_reads.fastq"
        if os.path.exists(prevalent_fastq):
            prevalent_fasta = f"{prevalent_fastq}.fasta"
            os.system(f"seqtk seq -a {prevalent_fastq} > {prevalent_fasta}")

            genes_main = os.popen(f"conda run -n {env_config['abricate']} abricate --db card --threads {threads} {prevalent_fasta} | cut -f6,10,11 | tail -n +2").read().strip("\n").split("\n")
            genepresence_main = []
            if genes_main == ['']:
                res["AMR_main"] = "NA"
            else:
                for g in genes_main:
                    g_parts = g.strip().split()
                    if len(g_parts) >= 3 and float(g_parts[1]) >= 80 and float(g_parts[2]) >= 95:
                        genepresence_main.append(g_parts[0])
            res["AMR_main"] = ";".join(list(set(genepresence_main))) if genepresence_main else "NA"
        else:
            res["AMR_main"] = "NA"

        # 5. Two-step KMA MLST directly on reads
        res["MLST_kma_reads"] = run_two_step_kma_mlst(chopselected, threads, local_scratch, shared_dir)

        # 6. Flye + BUSCO
        res["assembly_stats"] = {}

        for sp in assembly_taxa:
            sp_clean = sp.replace(" ", "_")
            sp_res = {}
            fastq_file = f"{local_fastq}_{sp_clean}_reads.fastq"
            flye_out = f"{local_fastq}_out_flye_{sp_clean}"

            try:
                os.system(f"conda run -n {env_config['flye']} flye --nano-raw {fastq_file} -o {flye_out} -t {threads}")
                asm_fasta = os.path.join(flye_out, "assembly.fasta")

                if os.path.exists(asm_fasta):
                    stats = os.popen(f'assembly-stats -s {asm_fasta} | cut -f2,3 | grep -e "total_length" -e "n50" -e "number" | cut -f2 | head -n3').read().strip().split("\n")
                    sp_res["genome-size"] = int(stats[0]) if len(stats) > 0 else "NA"
                    sp_res["num_contigs"] = int(stats[1]) if len(stats) > 1 else "NA"
                    sp_res["N50"] = int(stats[2]) if len(stats) > 2 else "NA"

                    sp_res["mlst"] = os.popen(f"conda run -n {env_config['mlst']} mlst {asm_fasta} | cut -f2,3").read().strip().replace("\t", " ")

                    genes = os.popen(f"conda run -n {env_config['abricate']} abricate --db card --threads {threads} {asm_fasta} | cut -f6,10,11 | tail -n +2").read().strip().split("\n")
                    g_pres = [g.split()[0] for g in genes if len(g.split()) >= 3 and float(g.split()[1]) >= 80 and float(g.split()[2]) >= 95]
                    sp_res["AMR_assembly"] = ";".join(set(g_pres)) if g_pres else "NA"

                    busco_name = f"busco_{sp_clean}"
                    busco_out_dir = os.path.join(local_scratch, busco_name)
                    busco_downloads = os.path.join(shared_dir, "busco_downloads/lineages/bacteria_odb12")

                    cmd_busco = (
                        f"cd {local_scratch} && "
                        f"busco -i {asm_fasta} -o {busco_name} -m geno "
                        f"--offline -l {busco_downloads} "
                        f"-c {threads} -f > busco_{sp_clean}.log 2>&1"
                    )
                    os.system(cmd_busco)

                    b_sums = glob.glob(f"{busco_out_dir}/short_summary.*.txt")

                    sp_res["BUSCO"] = "NA"
                    if b_sums:
                        with open(b_sums[0], "r") as bf:
                            for line in bf:
                                if "C:" in line and "[S:" in line:
                                    sp_res["BUSCO"] = line.strip().split("\t")[0]
                                    break

                    remote_assembly_path = os.path.join(assemblies_dir, f"T{tt}_{sample_name}_{sp_clean}.fasta")
                    if worker_node != manager_node:
                        os.system(f"scp -q {asm_fasta} {manager_node}:{remote_assembly_path}")
                    else:
                        os.system(f"cp {asm_fasta} {remote_assembly_path}")

            except Exception:
                sp_res.update({"genome-size": "NA", "num_contigs": "NA", "N50": "NA", "AMR_assembly": "NA", "mlst": "NA", "BUSCO": "NA"})

            res["assembly_stats"][sp] = sp_res

    finally:
        os.system(f"rm -rf {local_scratch}")

    return {sample_name: res}

## MAIN MANAGER
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="NANOResp Manager Pipeline")
    parser.add_argument("input_folder", help="Input folder containing FASTQ files")
    parser.add_argument("--mode", choices=["ontime", "live"], default="live")
    parser.add_argument("--timepoints", default="10,60,1440")
    parser.add_argument("--max-cpus", type=int, default=100)
    parser.add_argument("--max-mem", type=int, default=600)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--log-dir", required=True)
    parser.add_argument("--shared-dir", required=True)
    
    # Conda environment names parameters
    parser.add_argument("--conda-profile", default="/home/sgaiarsa/miniconda3/etc/profile.d/conda.sh")
    parser.add_argument("--env-main", default="ONT")
    parser.add_argument("--env-abricate", default="abricate")
    parser.add_argument("--env-flye", default="flye")
    parser.add_argument("--env-mlst", default="MLST")
    
    args = parser.parse_args()

    env_config = {
        'abricate': args.env_abricate,
        'flye': args.env_flye,
        'mlst': args.env_mlst
    }

    manager_node = socket.gethostname()
    base_dir = os.path.abspath(os.getcwd())

    run_num = "".join(filter(str.isdigit, os.path.basename(args.input_folder)))
    if not run_num: run_num = "X"
    json_out = f"results_R{run_num}.json"

    assemblies_dir = os.path.join(args.out_dir, f"assemblies_R{run_num}")
    os.makedirs(assemblies_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)

    MIN_CPUS_PER_WORKER = 12
    MIN_MEM_GB_PER_WORKER = 32

    max_workers_by_cpu = args.max_cpus // MIN_CPUS_PER_WORKER
    max_workers_by_mem = args.max_mem // MIN_MEM_GB_PER_WORKER
    max_parallel_workers = max(1, min(max_workers_by_cpu, max_workers_by_mem))

    threads_per_worker = MIN_CPUS_PER_WORKER
    mem_per_worker = f"{MIN_MEM_GB_PER_WORKER}GB"

    print(f"[MANAGER] Main Node: {manager_node}")
    print(f"[MANAGER] Allocation: {threads_per_worker} CPU, {mem_per_worker} RAM per Worker (Max Parallel Jobs: {max_parallel_workers}).")

    executor = submitit.AutoExecutor(folder=args.log_dir)
    executor.update_parameters(
        slurm_job_name="NanoWorker",
        slurm_mem=mem_per_worker,
        slurm_cpus_per_task=threads_per_worker,
        slurm_time="INFINITE",
        slurm_exclude="node1",
        slurm_setup=[
            f"source {args.conda_profile}",
            f"conda activate {args.env_main}"
        ]
    )

    timepoints = sorted([int(x.strip()) for x in args.timepoints.split(",")])
    results = {str(tp): {} for tp in timepoints}

    # =========================================================================
    # ONTIME MODE (Static files, 1 folder = N sample files)
    # =========================================================================
    if args.mode == "ontime":
        T0 = get_startingpoint(args.input_folder)
        all_tasks = []

        for tp in timepoints:
            tt_tmp_dir = os.path.join(base_dir, f"tmp_T{tp}")
            os.makedirs(tt_tmp_dir, exist_ok=True)

            extract_data(args.input_folder, tp, T0, tt_tmp_dir)

            samples = glob.glob(os.path.join(tt_tmp_dir, "*.fastq"))
            for sample_path in samples:
                all_tasks.append((sample_path, threads_per_worker, str(tp), manager_node, assemblies_dir, args.shared_dir, env_config))

        if all_tasks:
            executor.update_parameters(slurm_array_parallelism=max_parallel_workers)
            samples_l, threads_l, tt_l, manager_nodes_l, asm_dir_l, shared_dir_l, env_config_l = zip(*all_tasks)
            jobs = executor.map_array(process_single_sample, samples_l, threads_l, tt_l, manager_nodes_l, asm_dir_l, shared_dir_l, env_config_l)

            for job, tp_str in zip(jobs, tt_l):
                try:
                    res = job.result()
                    if isinstance(res, dict): results[tp_str].update(res)
                except Exception as e:
                    print(f"Job Error: {e}")

        for tp in timepoints:
            os.system(f"rm -rf {os.path.join(base_dir, f'tmp_T{tp}')}")

        with open(json_out, "w") as js:
            json.dump(results, js)

    # =========================================================================
    # LIVE MODE (1 folder = N subfolders with zcat + throttling)
    # =========================================================================
    elif args.mode == "live":
        tp_status = {tp: False for tp in timepoints}
        active_arrays = {} # format: array_id -> {'jobs': list_of_jobs, 'tp': tp}

        T0 = get_startingpoint(args.input_folder)
        start_time_str = datetime.fromtimestamp(T0).strftime("%Y-%m-%d %H:%M:%S")
        print(f"[MANAGER] T0 (Sequencing Start): {start_time_str}")

        try:
            while True:
                # 1. Fetch results from completed jobs for real-time JSON updates
                for arr_id, arr_data in list(active_arrays.items()):
                    tp = arr_data['tp']
                    still_active_jobs = []

                    for job, sample_path in arr_data['jobs']:
                        if job.done():
                            try:
                                res = job.result()
                                if isinstance(res, dict): results[str(tp)].update(res)
                            except Exception as e:
                                print(f"[MANAGER] Job Error for {sample_path}: {e}")
                        else:
                            still_active_jobs.append((job, sample_path))

                    if not still_active_jobs:
                        del active_arrays[arr_id]
                    else:
                        active_arrays[arr_id]['jobs'] = still_active_jobs

                with open(json_out, "w") as js:
                    json.dump(results, js)

                # 2. Exit if all timepoints are processed and completed
                if all(tp_status.values()) and not active_arrays:
                    print("[MANAGER] All timepoints processed and jobs concluded.")
                    break

                # 3. Trigger Timepoint based on elapsed minutes from T0
                elapsed_mins = (time.time() - T0) / 60
                for tp in timepoints:
                    if elapsed_mins >= tp and not tp_status[tp]:
                        print(f"[MANAGER] Reached Timepoint {tp} min. Concatenating sample subfolders...")
                        tp_status[tp] = True
                        out_tmp = os.path.join(base_dir, f"tmp_LIVE_T{tp}")

                        concat_files = concatenate_live_barcodes(args.input_folder, out_tmp)

                        if concat_files:
                            samples_l = concat_files
                            threads_l = [threads_per_worker] * len(concat_files)
                            tt_l = [str(tp)] * len(concat_files)
                            manager_nodes_l = [manager_node] * len(concat_files)
                            asm_dir_l = [assemblies_dir] * len(concat_files)
                            shared_dir_l = [args.shared_dir] * len(concat_files)
                            env_config_l = [env_config] * len(concat_files)

                            executor.update_parameters(slurm_array_parallelism=1)
                            jobs = executor.map_array(process_single_sample, samples_l, threads_l, tt_l, manager_nodes_l, asm_dir_l, shared_dir_l, env_config_l)

                            try:
                                array_id = jobs[0].job_id.split('_')[0]
                                active_arrays[array_id] = {'jobs': list(zip(jobs, samples_l)), 'tp': tp}
                                print(f"[MANAGER] Submitted Slurm Array ID: {array_id} for Timepoint T{tp}.")
                            except Exception as e:
                                print(f"[MANAGER] Error extracting array ID: {e}")

                # 4. Dynamic Throttling via squeue and scontrol
                if active_arrays:
                    total_running = 0
                    pending_arrays = []

                    for arr_id in active_arrays.keys():
                        res = subprocess.run(["squeue", "-j", str(arr_id), "-h", "-o", "%T"], capture_output=True, text=True)
                        states = res.stdout.strip().split('\n')

                        running = states.count("RUNNING")
                        pending = states.count("PENDING")
                        total_running += running

                        if pending > 0:
                            pending_arrays.append(arr_id)

                    free_slots = max(0, max_parallel_workers - total_running)

                    for arr_id in pending_arrays:
                        if free_slots > 0:
                            subprocess.run(["scontrol", "update", f"JobID={arr_id}", f"ArrayTaskThrottle={free_slots}"])
                            free_slots = 0
                        else:
                            subprocess.run(["scontrol", "update", f"JobID={arr_id}", f"ArrayTaskThrottle=0"])

                time.sleep(30)

        except KeyboardInterrupt:
            print("\n[MANAGER] LIVE Watchdog interrupted by user.")

    print(f"[MANAGER] Pipeline completed for Run {run_num}.")
