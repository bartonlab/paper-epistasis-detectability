using CSV, DataFrames, Printf, Statistics

# Exports figure data for pairwise co-frequency analysis.
# Case: N=1000, mu=0.001, sigma=0.1, r=0
#
# Output files in data/out_manuscript_expanded/N-1000_mu-0.00100/:
#   fig_auroc_sigma-{σ}_rec-{r}.csv          — time, Mean, Error_L, Error_H
#   fig_pos_pair_traj_sigma-{σ}_rec-{r}.csv  — i, j, time, x2_mean, z_score
#   fig_x2_group_summary_sigma-{σ}_rec-{r}.csv — group, time, x2_mean, x2_std
#
# Usage: julia export_fig_data_pairwise_Jun2026.jl [sigma_str rec_str]

σ_str   = length(ARGS) >= 1 ? ARGS[1] : "0.10000"
rec_str = length(ARGS) >= 2 ? ARGS[2] : "0.00000"

const data_dir = "__DATA_DIR__"

function main()
    # ── 1. AUROC: filter from accuracy_perRep.csv ─────────────────────────────
    acc = CSV.read(joinpath(data_dir, "accuracy_perRep.csv"), DataFrame)
    auroc = filter(
        r -> r.method == "LD_Temp" && r.metric == "auroc" &&
             isapprox(r.sigma, parse(Float64, σ_str), atol=1e-6) &&
             isapprox(r.rec,   parse(Float64, rec_str), atol=1e-6),
        acc)
    sort!(auroc, :time)
    auroc_out = select(auroc, :time, :Mean, :Error_L, :Error_H)

    # ── 2. Trajectory file ────────────────────────────────────────────────────
    traj = CSV.read(
        joinpath(data_dir, "pairwise_freq_traj_v4_sigma-$(σ_str)_rec-$(rec_str).csv"),
        DataFrame)

    off = filter(r -> r.i != r.j, traj)

    # Positive pairs: rep-averaged x2
    pos_raw = filter(r -> r.GT > 0.5, off)
    pos_avg = combine(groupby(pos_raw, [:i, :j, :time]),
        :x2 => mean => :x2_mean)
    sort!(pos_avg, [:i, :j, :time])
    pos_pairs = unique(pos_avg[:, [:i, :j]])
    @printf("Positive pairs: %d\n", nrow(pos_pairs))

    # Non-positive pairs: mean and std (over all pairs × replicates)
    nonpos_raw = filter(r -> r.GT <= 0.5, off)
    nonpos_summary = combine(groupby(nonpos_raw, :time),
        :x2 => mean => :x2_mean,
        :x2 => std  => :x2_std)
    sort!(nonpos_summary, :time)

    # ── 3. Z-score per positive pair ─────────────────────────────────────────
    nonpos_ref = rename(copy(nonpos_summary),
        :x2_mean => :mu_nonpos, :x2_std => :sig_nonpos)
    pos_z = innerjoin(pos_avg, nonpos_ref, on=:time)
    pos_z.z_score = (pos_z.x2_mean .- pos_z.mu_nonpos) ./ pos_z.sig_nonpos
    sort!(pos_z, [:i, :j, :time])

    pair_traj_out = select(pos_z, :i, :j, :time, :x2_mean, :z_score)

    # ── 4. Group summary: positive and non-positive in same format ────────────
    # Positive group: mean/std over all (pairs × reps) at each time
    pos_group = combine(groupby(pos_raw, :time),
        :x2 => mean => :x2_mean,
        :x2 => std  => :x2_std)
    sort!(pos_group, :time)
    pos_group.group .= "positive"

    nonpos_group = select(nonpos_summary, :time, :x2_mean, :x2_std)
    nonpos_group.group .= "non-positive"

    group_summary_out = vcat(
        select(pos_group,    :group, :time, :x2_mean, :x2_std),
        select(nonpos_group, :group, :time, :x2_mean, :x2_std))
    sort!(group_summary_out, [:group, :time])

    # ── Write ─────────────────────────────────────────────────────────────────
    f_auroc   = joinpath(data_dir, "fig_auroc_sigma-$(σ_str)_rec-$(rec_str).csv")
    f_pair    = joinpath(data_dir, "fig_pos_pair_traj_sigma-$(σ_str)_rec-$(rec_str).csv")
    f_summary = joinpath(data_dir, "fig_x2_group_summary_sigma-$(σ_str)_rec-$(rec_str).csv")

    CSV.write(f_auroc,   auroc_out)
    CSV.write(f_pair,    pair_traj_out)
    CSV.write(f_summary, group_summary_out)

    @printf("Written:\n")
    @printf("  %d rows  →  %s\n", nrow(auroc_out),       basename(f_auroc))
    @printf("  %d rows  →  %s\n", nrow(pair_traj_out),   basename(f_pair))
    @printf("  %d rows  →  %s\n", nrow(group_summary_out), basename(f_summary))
end

main()
