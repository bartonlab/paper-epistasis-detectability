#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name="accuracy"
##SBATCH -p <partition>        # set your cluster partition
#SBATCH --array=1-5440%500

# Step 2a: Compute AUROC and related accuracy metrics for positive epistasis.
# Reads WF ground-truth files and inference output from step 1.
# Writes per-(sigma,rec) accuracy summary CSVs used by step 3.
#
# Array dimensions: 5 N × 4 mu × 16 sigma × 17 rec = 5440 jobs
# Columns of job_params.txt: job_id  i_N  i_mu  i_sigma  i_rec

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
JULIA=/path/to/julia/bin/julia
SRC_DIR=/path/to/paper-epistasis-detectability/src
WF_DIR=/path/to/wf_output          # input:  WF simulation trajectories (ground truth)
INF_DIR=/path/to/inference_output  # input:  inference results from step 1
ACC_DIR=/path/to/accuracy_output   # output: per-param accuracy CSVs
# ─────────────────────────────────────────────────────────────────────────────

date
hostname

params=$(dirname "$0")/job_params.txt
i_N=$(  awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $2}' $params)
i_mu=$( awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $3}' $params)
i_sig=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $4}' $params)
i_rec=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $5}' $params)

echo "array_id=${SLURM_ARRAY_TASK_ID}  i_N=${i_N}  i_mu=${i_mu}  i_sig=${i_sig}  i_rec=${i_rec}"

$JULIA $SRC_DIR/get_accuracy.jl \
    $i_N $i_mu $i_sig $i_rec $WF_DIR $INF_DIR $ACC_DIR

date
exit 0
