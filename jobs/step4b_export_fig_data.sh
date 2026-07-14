#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --job-name="fig-data"
##SBATCH -p <partition>        # set your cluster partition

# Step 4b: Export figure data CSVs for the 4-panel pairwise frequency figure.
#
# Fixed case: N=1000, mu=0.001. Sigma and rec are set below.
# Reads:
#   DATA_DIR/accuracy_perRep.csv                              (from step 3)
#   DATA_DIR/pairwise_freq_traj_v4_sigma-{SIGMA}_rec-{REC}.csv  (from step 4a)
# Writes:
#   DATA_DIR/fig_auroc_sigma-{SIGMA}_rec-{REC}.csv
#   DATA_DIR/fig_pos_pair_traj_sigma-{SIGMA}_rec-{REC}.csv
#   DATA_DIR/fig_x2_group_summary_sigma-{SIGMA}_rec-{REC}.csv
#
# Run AFTER step3 and step4a.
#
# To run interactively (no SBATCH):
#   bash step4b_export_fig_data.sh
#
# To submit to SLURM:
#   sbatch step4b_export_fig_data.sh

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
JULIA=/path/to/julia/bin/julia
SRC_DIR=/path/to/paper-epistasis-detectability/src
DATA_DIR=/path/to/output/N-1000_mu-0.00100  # directory containing accuracy_perRep.csv
                                              # and pairwise_freq_traj_v4_*.csv

SIGMA="0.10000"   # epistasis strength (5 decimal places)
REC="0.00000"     # recombination rate (5 decimal places)
# ─────────────────────────────────────────────────────────────────────────────

# Patch hardcoded path in a temporary copy of the script
TMPSCRIPT=$(mktemp /tmp/export_fig_data_XXXX.jl)
sed -e "s|__DATA_DIR__|${DATA_DIR}|g" \
    $SRC_DIR/export_fig_data.jl > $TMPSCRIPT

date
echo "sigma=${SIGMA}  rec=${REC}"
echo "DATA_DIR=${DATA_DIR}"

$JULIA $TMPSCRIPT $SIGMA $REC

rm -f $TMPSCRIPT
date
exit 0
