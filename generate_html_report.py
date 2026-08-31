import json
import sys
import os

def build_html(json_file, output_html):
    """Parses JSON results and generates an interactive, localized English HTML report."""
    with open(json_file, 'r') as f:
        data = json.load(f)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NANOResp - Results Viewer</title>
    <!-- Chart.js CDN for donut charts -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        :root {{
            --bg: #f8fafc;
            --card-bg: #ffffff;
            --text: #0f172a;
            --border: #e2e8f0;
            --primary: #2563eb;
            --primary-light: #eff6ff;
            --accent: #059669;
            --warning: #d97706;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg);
            color: var(--text);
            margin: 0;
            padding: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        header {{
            background: var(--card-bg);
            padding: 20px 24px;
            border-radius: 12px;
            border: 1px solid var(--border);
            margin-bottom: 24px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        }}
        h1 {{ margin: 0 0 16px 0; font-size: 1.5rem; color: #1e293b; }}
        .controls {{
            display: flex;
            gap: 16px;
            flex-wrap: wrap;
        }}
        .control-group {{
            display: flex;
            flex-direction: column;
            gap: 6px;
        }}
        label {{ font-size: 0.85rem; font-weight: 600; color: #64748b; }}
        select {{
            padding: 8px 12px;
            border-radius: 6px;
            border: 1px solid var(--border);
            background: #fff;
            font-size: 0.95rem;
            min-width: 200px;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .card {{
            background: var(--card-bg);
            padding: 16px 20px;
            border-radius: 10px;
            border: 1px solid var(--border);
            box-shadow: 0 1px 2px rgba(0,0,0,0.03);
        }}
        .card-title {{ font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.05em; color: #64748b; font-weight: 700; }}
        .card-value {{ font-size: 1.5rem; font-weight: 700; color: #0f172a; margin-top: 4px; }}
        .section-title {{
            font-size: 1.1rem;
            font-weight: 700;
            margin: 28px 0 14px 0;
            color: #334155;
            border-bottom: 2px solid var(--border);
            padding-bottom: 6px;
        }}
        .taxa-container {{
            display: grid;
            grid-template-columns: 320px 1fr;
            gap: 20px;
            background: var(--card-bg);
            padding: 20px;
            border-radius: 10px;
            border: 1px solid var(--border);
            margin-bottom: 24px;
            align-items: center;
        }}
        .chart-box {{
            position: relative;
            width: 100%;
            max-width: 300px;
            margin: 0 auto;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: var(--card-bg);
            border-radius: 8px;
            overflow: hidden;
            border: 1px solid var(--border);
            margin-bottom: 20px;
        }}
        th, td {{
            padding: 10px 14px;
            text-align: left;
            border-bottom: 1px solid var(--border);
            font-size: 0.88rem;
        }}
        th {{ background: #f1f5f9; font-weight: 600; color: #475569; }}
        tr:last-child td {{ border-bottom: none; }}
        .badge {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 4px;
            font-size: 0.78rem;
            font-weight: 600;
            background: #e2e8f0;
            color: #334155;
            word-break: break-all;
            margin: 1px;
        }}
        .badge-green {{ background: #d1fae5; color: #065f46; }}
        .badge-blue {{ background: #dbeafe; color: #1e40af; }}
        .badge-orange {{ background: #fef3c7; color: #92400e; }}
        .badge-purple {{ background: #f3e8ff; color: #6b21a8; }}
        .progress-bar {{
            height: 6px;
            background: #e2e8f0;
            border-radius: 4px;
            overflow: hidden;
            margin-top: 4px;
        }}
        .progress-fill {{ height: 100%; background: var(--primary); }}
        @media (max-width: 768px) {{
            .taxa-container {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>NANOResp Analysis Report</h1>
            <div class="controls">
                <div class="control-group">
                    <label for="timepointSelect">Timepoint (Minutes)</label>
                    <select id="timepointSelect" onchange="onTimepointChange()"></select>
                </div>
                <div class="control-group">
                    <label for="sampleSelect">Sample</label>
                    <select id="sampleSelect" onchange="renderDetails()"></select>
                </div>
            </div>
        </header>

        <div id="content">
            <!-- Metric Cards -->
            <div class="grid">
                <div class="card">
                    <div class="card-title">Total Reads</div>
                    <div class="card-value" id="valTotalReads">-</div>
                </div>
                <div class="card">
                    <div class="card-title">Filtered Human Reads</div>
                    <div class="card-value" id="valHumanPct">-</div>
                </div>
                <div class="card">
                    <div class="card-title">Non-Human Reads</div>
                    <div class="card-value" id="valNonHumanReads">-</div>
                </div>
                <div class="card">
                    <div class="card-title">Prevalent Species</div>
                    <div class="card-value" id="valMainTaxon" style="font-size: 1.1rem;">-</div>
                </div>
            </div>

            <!-- Host Filtering Details -->
            <div class="section-title">1. Host Filtering (Human vs Bacterial)</div>
            <table>
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                    </tr>
                </thead>
                <tbody id="hostTable"></tbody>
            </table>

            <!-- Kraken Taxonomy -->
            <div class="section-title">2. Taxonomic Profile (Kraken2)</div>
            <div class="taxa-container">
                <div class="chart-box">
                    <canvas id="taxaChart"></canvas>
                </div>
                <div>
                    <table>
                        <thead>
                            <tr>
                                <th>Species (Top Hits)</th>
                                <th>Reads</th>
                                <th>% Non-Human</th>
                            </tr>
                        </thead>
                        <tbody id="taxaTable"></tbody>
                    </table>
                </div>
            </div>

            <!-- MLST Comparison -->
            <div class="section-title">3. MLST Typing (Direct on Reads)</div>
            <table>
                <thead>
                    <tr>
                        <th>Method</th>
                        <th>Detected Scheme</th>
                        <th>Identified Alleles</th>
                    </tr>
                </thead>
                <tbody id="mlstTable"></tbody>
            </table>

            <!-- AMR -->
            <div class="section-title">4. Antimicrobial Resistance from Reads (CARD)</div>
            <table>
                <thead>
                    <tr>
                        <th>Target Alignment</th>
                        <th>Detected Resistance Genes</th>
                    </tr>
                </thead>
                <tbody id="amrTable"></tbody>
            </table>

            <!-- Assembly & Quality -->
            <div class="section-title">5. Assembly Features & Quality (Flye + BUSCO + MLST + AMR)</div>
            <table>
                <thead>
                    <tr>
                        <th>Assembled Species</th>
                        <th>Size (bp)</th>
                        <th>Contigs</th>
                        <th>N50 (bp)</th>
                        <th>BUSCO Score</th>
                        <th>Assembly MLST</th>
                        <th>AMR Genes (ABRicate/CARD)</th>
                    </tr>
                </thead>
                <tbody id="assemblyTable"></tbody>
            </table>
        </div>
    </div>

    <script>
        const rawData = {json.dumps(data)};
        let chartInstance = null;

        function init() {{
            const tpSelect = document.getElementById('timepointSelect');
            const timepoints = Object.keys(rawData);

            tpSelect.innerHTML = '';
            timepoints.forEach(tp => {{
                const opt = document.createElement('option');
                opt.value = tp;
                opt.textContent = `Timepoint: T${{tp}}m`;
                tpSelect.appendChild(opt);
            }});

            if (timepoints.length > 0) {{
                onTimepointChange();
            }}
        }}

        function onTimepointChange() {{
            const tp = document.getElementById('timepointSelect').value;
            const sampleSelect = document.getElementById('sampleSelect');
            const samples = Object.keys(rawData[tp] || {{}});

            sampleSelect.innerHTML = '';
            samples.forEach(s => {{
                const opt = document.createElement('option');
                opt.value = s;
                opt.textContent = s;
                sampleSelect.appendChild(opt);
            }});

            if (samples.length > 0) {{
                renderDetails();
            }}
        }}

        function renderDetails() {{
            const tp = document.getElementById('timepointSelect').value;
            const sample = document.getElementById('sampleSelect').value;
            const d = rawData[tp]?.[sample];

            if (!d) return;

            // 1. Metric Cards
            const totalReads = d.num_reads || 0;
            const nonHumanReads = d.num_nohuman_reads || 0;
            const humanReads = Math.max(0, totalReads - nonHumanReads);
            const humanPct = totalReads > 0 ? ((humanReads / totalReads) * 100).toFixed(1) : 0;

            document.getElementById('valTotalReads').textContent = totalReads.toLocaleString();
            document.getElementById('valHumanPct').textContent = `${{humanPct}}%`;
            document.getElementById('valNonHumanReads').textContent = nonHumanReads.toLocaleString();
            document.getElementById('valMainTaxon').textContent = d.main || 'N/A';

            // 2. Host Table
            document.getElementById('hostTable').innerHTML = `
                <tr><td>Raw Total Reads</td><td><b>${{totalReads.toLocaleString()}}</b></td></tr>
                <tr><td>Filtered Human Reads</td><td>${{humanReads.toLocaleString()}} (${{humanPct}}%)</td></tr>
                <tr><td>Retained Bacterial/Non-Human Reads</td><td>${{nonHumanReads.toLocaleString()}} (${{(100 - humanPct).toFixed(1)}}%)</td></tr>
                <tr><td>Total Non-Human Bases</td><td>${{(d.nohuman_bases || 0).toLocaleString()}} bp</td></tr>
            `;

            // 3. Taxa Table & Donut Chart
            const taxa = d.taxa || {{}};
            const sortedTaxa = Object.entries(taxa).sort((a,b) => b[1] - a[1]);

            const topN = sortedTaxa.slice(0, 7);
            const others = sortedTaxa.slice(7);
            const othersCount = others.reduce((acc, curr) => acc + curr[1], 0);

            let chartLabels = topN.map(x => x[0]);
            let chartData = topN.map(x => x[1]);

            if (othersCount > 0) {{
                chartLabels.push('Other minor species');
                chartData.push(othersCount);
            }}

            renderTaxaChart(chartLabels, chartData);

            let taxaRows = '';
            topN.forEach(([sp, count]) => {{
                const pct = nonHumanReads > 0 ? ((count / nonHumanReads) * 100).toFixed(1) : 0;
                taxaRows += `
                    <tr>
                        <td><b>${{sp}}</b></td>
                        <td>${{count.toLocaleString()}}</td>
                        <td>
                            ${{pct}}%
                            <div class="progress-bar"><div class="progress-fill" style="width: ${{pct}}%"></div></div>
                        </td>
                    </tr>
                `;
            }});
            if (othersCount > 0) {{
                const othersPct = nonHumanReads > 0 ? ((othersCount / nonHumanReads) * 100).toFixed(1) : 0;
                taxaRows += `
                    <tr>
                        <td><i>Other ${{others.length}} minor species...</i></td>
                        <td>${{othersCount.toLocaleString()}}</td>
                        <td>${{othersPct}}%</td>
                    </tr>
                `;
            }}
            document.getElementById('taxaTable').innerHTML = taxaRows || '<tr><td colspan="3">No taxa identified</td></tr>';

            // 4. MLST Table (Reads only)
            const mlstReads = d.MLST_kma_reads || {{}};
            const scheme = mlstReads.detected_scheme || 'NA';
            const alleles = mlstReads.alleles || 'NA';

            document.getElementById('mlstTable').innerHTML = `
                <tr>
                    <td><span class="badge badge-blue">KMA 2-step (Reads)</span></td>
                    <td><b>${{scheme}}</b></td>
                    <td><span class="badge">${{alleles.replace(/;/g, ', ')}}</span></td>
                </tr>
            `;

            // 5. AMR Table (Reads only)
            document.getElementById('amrTable').innerHTML = `
                <tr>
                    <td>Reads Overall (All reads)</td>
                    <td>${{formatBadges(d.AMR_overall)}}</td>
                </tr>
                <tr>
                    <td>Prevalent Reads (${{d.main || 'N/A'}})</td>
                    <td>${{formatBadges(d.AMR_main)}}</td>
                </tr>
            `;

            // 6. Unified Assembly Section
            const asmStats = d.assembly_stats || {{}};
            let asmRows = '';
            Object.entries(asmStats).forEach(([sp, stats]) => {{
                if (!stats || Object.keys(stats).length === 0) {{
                    asmRows += `
                        <tr>
                            <td><b>${{sp}}</b></td>
                            <td colspan="6"><i>Assembly failed or not executed</i></td>
                        </tr>
                    `;
                    return;
                }}
                asmRows += `
                    <tr>
                        <td><b>${{sp}}</b></td>
                        <td>${{stats['genome-size'] && stats['genome-size'] !== 'NA' ? stats['genome-size'].toLocaleString() : 'NA'}}</td>
                        <td>${{stats.num_contigs || 'NA'}}</td>
                        <td>${{stats.N50 && stats.N50 !== 'NA' ? stats.N50.toLocaleString() : 'NA'}}</td>
                        <td><span class="badge badge-orange">${{stats.BUSCO || 'NA'}}</span></td>
                        <td><span class="badge badge-blue">${{stats.mlst || 'NA'}}</span></td>
                        <td>${{formatBadges(stats.AMR_assembly)}}</td>
                    </tr>
                `;
            }});
            document.getElementById('assemblyTable').innerHTML = asmRows || '<tr><td colspan="7">No assembly available</td></tr>';
        }}

        function renderTaxaChart(labels, dataCounts) {{
            const ctx = document.getElementById('taxaChart').getContext('2d');
            if (chartInstance) {{
                chartInstance.destroy();
            }}
            chartInstance = new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: labels,
                    datasets: [{{
                        data: dataCounts,
                        backgroundColor: [
                            '#2563eb', '#059669', '#d97706', '#7c3aed',
                            '#dc2626', '#0891b2', '#4b5563', '#cbd5e1'
                        ]
                    }}]
                }},
                options: {{
                    responsive: true,
                    plugins: {{
                        legend: {{
                            display: false
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const val = context.raw;
                                    const pct = ((val / total) * 100).toFixed(1);
                                    return ` ${{context.label}}: ${{val.toLocaleString()}} (${{pct}}%)`;
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }}

        function formatBadges(str) {{
            if (!str || str === 'NA') return '<span class="badge">None</span>';
            return str.split(';').map(g => `<span class="badge badge-green">${{g}}</span>`).join(' ');
        }}

        window.onload = init;
    </script>
</body>
</html>
"""
    with open(output_html, 'w') as out:
        out.write(html_content)
    print(f"[REPORT] HTML Report successfully generated: {output_html}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_html_report.py results_R1.json [output.html]")
        sys.exit(1)

    json_in = sys.argv[1]
    html_out = sys.argv[2] if len(sys.argv) > 2 else json_in.replace('.json', '.html')
    build_html(json_in, html_out)
