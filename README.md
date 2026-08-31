# NANOResp
A pipeline for real time analysis of bacterial Nanopore metagenomics of clinical samples. The pipeline was designed for the analysis of respiratory samples and optimized for HPC clusters (slurm manager).

Scripts and configuration files used for manuscript XXXX

## Instructions:

1. Install che conda environments using the txt configuration files
```
conda create --name ONT --file env_ONT_md5.txt
conda create --name flye --file env_flye_md5.txt
conda create --name MLST --file env_MLST_md5.txt
conda create --name abricate --file env_abricate_md5.txt

```
2. Download the human genome reference:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna Hg38.fna
```
3. Download and set the Kraken DB:
```
mkdir -p krakenDB/standard16/
cd krakenDB/standard16/
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20260626.tar.gz
tar -zxf k2_standard_16_GB_20260626.tar.gz
rm k2_standard_16_GB_20260626.tar.gz
```
4. Download additional files and folders:
```
conda activate ONT
gdown https://drive.google.com/file/d/10avumi4817w3iahEadj6AfkM_Z_jWtGt/view?usp=drive_link
tar -zxf Additional_folders.tar.gz
```
5. Set all variables in the manager.sh script header. Most importantly, you can decide to run the script during live sequencing (basecalling not included) or on after sequencing and basecalling using the ontime tool. Please note: the variables are currently set to work on our cluster, please take some time to adapt it to yours.
6. Usage: `sbatch manager.sh Run1 [Run2 Run3 Run4]`
