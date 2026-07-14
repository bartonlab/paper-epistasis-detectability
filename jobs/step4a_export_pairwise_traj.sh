#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name="pairwise-traj"
##SBATCH -p <partition>        # set your cluster partition

# Step 4a: Export per-replicate pairwise co-frequency trajectories x2[i,j].
#
# Fixed case: N=1000, mu=0.001. Sigma and rec are set below.
# Output: OUT_DIR/N-1000_mu-0.00100/pairwise_freq_traj_v4_sigma-{SIGMA}_rec-{REC}.csv
#
# Columns: i, j, id, time, x2, GT
#   diagonal (i==j): x2 = allele frequency x1[i],       GT = S_GT[i]
#   off-diagonal(i<j): x2 = double-mutant co-frequency,  GT = J_GT[i,j]
#
# Run BEFORE step4b. Must be run after step1 (inference output required).
#
# To run interactively (no SBATCH):
#   bash step4a_export_pairwise_traj.sh
#
# To submit to SLURM:
#   sbatch step4a_export_pairwise_traj.sh

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
JULIA=/path/to/julia/bin/julia
SRC_DIR=/path/to/paper-epistasis-detectability/src
INF_DIR=/path/to/inference_output   # input: inference output from step 1
WF_DIR=/path/to/wf_output           # input: WF simulation output (ground truth)
OUT_DIR=/path/to/output             # output: directory for the trajectory CSV

SIGMA="0.10000"   # epistasis strength (5 decimal places)
REC="0.00000"     # recombination rate (5 decimal places)
# ─────────────────────────────────────────────────────────────────────────────

# Patch hardcoded paths in a temporary copy of the script
TMPSCRIPT=$(mktemp /tmp/export_pairwise_traj_XXXX.jl)
sed -e "s|__INF_DIR__|${INF_DIR}|g" \
    -e "s|__WF_DIR__|${WF_DIR}|g" \
    -e "s|__OUT_DIR__|${OUT_DIR}|g" \
    $SRC_DIR/export_pairwise_traj.jl > $TMPSCRIPT

date
echo "sigma=${SIGMA}  rec=${REC}"
echo "INF_DIR=${INF_DIR}"
echo "WF_DIR=${WF_DIR}"
echo "OUT_DIR=${OUT_DIR}"

$JULIA $TMPSCRIPT $SIGMA $REC

rm -f $TMPSCRIPT
date
exit 0
