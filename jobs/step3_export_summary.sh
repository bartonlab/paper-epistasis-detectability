#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --job-name="export-summary"
##SBATCH -p <partition>        # set your cluster partition
#SBATCH --array=1-20

# Step 3: Aggregate per-(sigma,rec) accuracy CSVs into one summary CSV per (N,mu).
# Reads accuracy output from step 2a and FDR output from step 2b.
# Writes accuracy_perRep.csv and physics-metrics.csv to OUT_DIR/N-{N}_mu-{mu}/.
#
# Array mapping (one job per (N, mu) combination, 5 N × 4 mu = 20 jobs):
#   ID  1- 4: N=100,   mu=1e-4/5e-4/1e-3/5e-3
#   ID  5- 8: N=500,   mu=1e-4/5e-4/1e-3/5e-3
#   ID  9-12: N=1000,  mu=1e-4/5e-4/1e-3/5e-3
#   ID 13-16: N=5000,  mu=1e-4/5e-4/1e-3/5e-3
#   ID 17-20: N=10000, mu=1e-4/5e-4/1e-3/5e-3

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
JULIA=/path/to/julia/bin/julia
SRC_DIR=/path/to/paper-epistasis-detectability/src
ACC_DIR=/path/to/accuracy_output   # input:  accuracy CSVs from step 2a
INF_DIR=/path/to/inference_output  # input:  inference output from step 1 (for physics-metrics)
OUT_DIR=/path/to/output            # output: aggregated summary CSVs
# ─────────────────────────────────────────────────────────────────────────────

date
hostname

i_N=$(( ($SLURM_ARRAY_TASK_ID - 1) / 4 + 1 ))
i_mu=$(( ($SLURM_ARRAY_TASK_ID - 1) % 4 + 1 ))

echo "array_id=${SLURM_ARRAY_TASK_ID}  i_N=${i_N}  i_mu=${i_mu}"

# Patch hardcoded paths in a temporary copy of the script
TMPSCRIPT=$(mktemp /tmp/export_summary_XXXX.jl)
sed -e "s|__ACC_DIR__|${ACC_DIR}|g" \
    -e "s|__INF_DIR__|${INF_DIR}|g" \
    -e "s|__OUT_DIR__|${OUT_DIR}|g" \
    $SRC_DIR/export_summary.jl > $TMPSCRIPT

$JULIA $TMPSCRIPT $i_N $i_mu

rm -f $TMPSCRIPT
date
exit 0
