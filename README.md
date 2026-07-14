
# Overview

This repository contains scripts and notebooks for reproducing the results accompanying the manuscript

### Title of paper
Kai Shimagaki<sup>1, 2, 3</sup> and John P. Barton<sup>1,2,#</sup>

<sup>1</sup> Department of Physics and Astronomy, University of California, Riverside  
<sup>2</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine  
<sup>3</sup> Michigan Center for Applied and Interdisciplinary Mathematics, College of Literature, Science, and the Arts, University of Michigan, Michigan, USA.  
<sup>#</sup> correspondence to [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)



# Directory structure

```
src/    Julia scripts implementing the analysis pipeline
jobs/   SLURM batch scripts and parameter table for running src/ scripts on a cluster
note/   Jupyter notebooks for figure generation
data/   Output CSVs produced by the pipeline (populated after running the scripts)
```

# Pipeline

The analysis runs in four sequential steps. Each step is executed by submitting the
corresponding script in `jobs/` via `sbatch`, or running it directly with `bash` for
single-job steps. Before submitting, set the path variables in the `USER CONFIGURATION`
block at the top of each job script.

### Step 1 — Inference (`jobs/step1_run_inference.sh`)

Runs epistasis inference on Wright–Fisher simulation trajectories.  
Julia script: `src/run_inference.jl`  
Input: WF simulation output (allele-frequency trajectories, ground-truth coupling files)  
Output: per-replicate inferred epistasis and selection matrices for all methods (LD, KNS, KNSGC, UFE, MPL, SL), fitness-wave files, and diversity time courses

```bash
sbatch jobs/step1_run_inference.sh   # SLURM array: 5440 jobs (5 N × 4 μ × 16 σ × 17 rec)
```

### Step 2a — Accuracy / AUROC (`jobs/step2a_get_accuracy.sh`)

Computes AUROC and related accuracy metrics for positive epistasis detection.  
Julia script: `src/get_accuracy.jl`  
Input: WF ground-truth files, inference output from step 1  
Output: per-(σ, rec) accuracy summary CSVs

```bash
sbatch jobs/step2a_get_accuracy.sh   # SLURM array: 5440 jobs
```

### Step 2b — Statistical power / TPR (`jobs/step2b_get_fdr_auprc.sh`)

Computes TPR at 5% FPR (statistical power) for epistasis detection.  
Julia script: `src/get_fdr_auprc.jl`  
Input: WF ground-truth files, inference output from step 1  
Output: per-(σ, rec) FDR/AUPRC summary CSVs

```bash
sbatch jobs/step2b_get_fdr_auprc.sh  # SLURM array: 5440 jobs
```

### Step 3 — Aggregation (`jobs/step3_export_summary.sh`)

Aggregates per-(σ, rec) CSVs into one summary file per (N, μ).  
Julia script: `src/export_summary.jl`  
Input: accuracy output from step 2a, inference output from step 1  
Output (written to `data/N-{N}_mu-{mu}/`):

- `accuracy_perRep.csv` — AUROC and accuracy metrics for all methods, times, and parameter combinations
- `physics-metrics.csv` — physical observables (fitness wave, LD, diversity)

```bash
sbatch jobs/step3_export_summary.sh  # SLURM array: 20 jobs (5 N × 4 μ)
```

### Step 4a — Pairwise frequency trajectories (`jobs/step4a_export_pairwise_traj.sh`)

Exports per-replicate pairwise co-frequency trajectories x₂[i,j] for a single parameter
combination (default: N=1000, μ=10⁻³, σ=0.1, r=0).  
Julia script: `src/export_pairwise_traj.jl`  
Input: inference output from step 1, WF ground-truth files  
Output (written to `data/N-1000_mu-0.00100/`):

- `pairwise_freq_traj_v4_sigma-{σ}_rec-{r}.csv` — i, j, replicate id, time, x₂, ground-truth coupling

```bash
bash jobs/step4a_export_pairwise_traj.sh   # single job; can also use sbatch
```

### Step 4b — Figure data (`jobs/step4b_export_fig_data.sh`)

Computes AUROC time series, positive-pair trajectories, and z-scores from the pairwise
trajectory file produced in step 4a.  
Julia script: `src/export_fig_data.jl`  
Input: `accuracy_perRep.csv` from step 3, pairwise trajectory CSV from step 4a  
Output (written to `data/N-1000_mu-0.00100/`):

- `fig_auroc_sigma-{σ}_rec-{r}.csv` — AUROC time series (mean ± error)
- `fig_pos_pair_traj_sigma-{σ}_rec-{r}.csv` — per-pair x₂ trajectories and z-scores
- `fig_x2_group_summary_sigma-{σ}_rec-{r}.csv` — group mean ± std for positive and non-positive pairs

```bash
bash jobs/step4b_export_fig_data.sh        # single job; can also use sbatch
```

# Notebooks

`note/figure.ipynb` — generates manuscript figures from the CSV files in `data/`

# Software dependencies

- [Julia](https://julialang.org/) 1.8 or later  
  Packages: `CSV`, `DataFrames`, `DelimitedFiles`, `Statistics`, `LinearAlgebra`, `CodecZlib`, `Printf`, `StatsBase`
- [Jupyter](https://jupyter.org/) with a Julia kernel for running the notebooks

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
