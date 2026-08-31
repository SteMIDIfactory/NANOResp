#!/bin/bash
#SBATCH --job-name=ONT_Manager
#SBATCH --output=ONT_Manager.out
#SBATCH --error=ONT_Manager.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --exclude=node1

# =========================================================================
# GLOBAL VARIABLES CONFIGURATION
# =========================================================================
# Slurm Orchestration Resources (Worker)
MAX_CPUS=240
MAX_MEM_GB=1200

# Pipeline Parameters
MODE="ontime"             # Options: "ontime" or "live"
TIMEPOINTS="10,20,30,60,1440" # Used in BOTH modes

# SYSTEM PATHS & CONDA ENVIRONMENTS (Modify these for deployment)
CONDA_PROFILE_PATH="/home/sgaiarsa/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_MAIN="ONT"

BASE_SCRIPT_DIR="/home/sgaiarsa/NANOResp"
INPUT_BASE_DIR="/home/sgaiarsa/NANOResp_inputs"
RESULTS_DIR="/home/sgaiarsa/NANORespRESULTS"
LOGS_DIR="/home/sgaiarsa/submitit_logs"
SHARED_DIR="/home/sgaiarsa/NANOResp"
# =========================================================================

if [ "$#" -eq 0 ]; then
    echo "[ERROR] Please specify at least one Run folder as an argument."
    echo "Usage: sbatch manager.sh Run1 [Run2 Run3 ...]"
    exit 1
fi

source "$CONDA_PROFILE_PATH"
conda activate "$CONDA_ENV_MAIN"

mkdir -p "$RESULTS_DIR"
mkdir -p "$LOGS_DIR"

for RUN in "$@"; do
    echo "=================================================="
    echo "[MANAGER] Starting processing for: $RUN (Mode: $MODE)"
    echo "=================================================="

    RUN_NUM=$(echo "$RUN" | sed 's/[^0-9]//g')
    SOURCE_DIR="$INPUT_BASE_DIR/$RUN"

    if [ ! -d "$SOURCE_DIR" ]; then
        echo "[ERROR] Source folder $SOURCE_DIR does not exist. Skipping."
        continue
    fi

    RUN_SCRATCH="$SCRATCH/NANOResp_$RUN"
    rm -rf "$RUN_SCRATCH"
    mkdir -p "$RUN_SCRATCH/$RUN"

    # Transfer raw data to Manager's scratch
    cp -r "$SOURCE_DIR"/* "$RUN_SCRATCH/$RUN/" 2>/dev/null || true
    cd "$RUN_SCRATCH"

    # Execute Python Pipeline with dynamic parameters
    python "$BASE_SCRIPT_DIR/Parse_seq2_submitit.py" "$RUN_SCRATCH/$RUN" \
        --mode "$MODE" \
        --timepoints "$TIMEPOINTS" \
        --max-cpus "$MAX_CPUS" \
        --max-mem "$MAX_MEM_GB" \
        --out-dir "$RESULTS_DIR" \
        --log-dir "$LOGS_DIR" \
        --shared-dir "$SHARED_DIR" \
        --conda-profile "$CONDA_PROFILE_PATH" \
        --env-main "$CONDA_ENV_MAIN"

    # Generate HTML Report from JSON
    JSON_FILE="$RUN_SCRATCH/results_R${RUN_NUM}.json"
    HTML_FILE="$RUN_SCRATCH/report_R${RUN_NUM}.html"

    if [ -f "$JSON_FILE" ]; then
        python "$BASE_SCRIPT_DIR/generate_html_report.py" "$JSON_FILE" "$HTML_FILE"
        cp "$JSON_FILE" "$RESULTS_DIR/"
        cp "$HTML_FILE" "$RESULTS_DIR/"
        echo "[MANAGER] Results and HTML Report saved to $RESULTS_DIR"
    fi

    # Cleanup scratch
    rm -rf "$RUN_SCRATCH"
    echo "[MANAGER] Completed processing for $RUN"
    echo ""
done

echo "[MANAGER] Processing completed for all inputs!"
