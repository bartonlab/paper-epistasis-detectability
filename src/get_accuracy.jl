using Random, Statistics, LinearAlgebra, Printf, CSV, DataFrames, DelimitedFiles, StatsBase, CodecZlib

include("epistasis_detectability.jl")

# ── I/O helpers: auto-detect plain .txt or gzip .txt.gz ──────────────────────
function find_output_file(base_no_ext)
    isfile(base_no_ext * ".txt.gz") && return base_no_ext * ".txt.gz"
    isfile(base_no_ext * ".txt")    && return base_no_ext * ".txt"
    return nothing
end

function readdlm_file(path)
    if endswith(path, ".gz")
        open(GzipDecompressorStream, path) do io
            readdlm(io)
        end
    else
        readdlm(path)
    end
end

Random.seed!(1234)

# ── Arguments ─────────────────────────────────────────────────────────────────
# Usage: julia get_accuracy_expanded_grid_Jun2026_v5.jl \
#            i_N i_mu i_σ i_rec wf_base_dir inf_base_dir acc_base_dir
#
# Changes from v4:
#   1. α_locus threshold removed for selection: idx_s = collect(1:L)
#      (all loci included in Selection_MPL/Selection_SL; epistasis thresholds unchanged)
#   2. Output dir: accuracy_output_expanded_v5_Jun2026
i_N          = parse(Int, ARGS[1])
i_mu         = parse(Int, ARGS[2])
i_σ          = parse(Int, ARGS[3])
i_rec        = parse(Int, ARGS[4])
wf_base_dir  = ARGS[5]
inf_base_dir = ARGS[6]
acc_base_dir = ARGS[7]

# ── Parameter grids ───────────────────────────────────────────────────────────
list_N   = [100, 500, 1000, 5000, 10000]
list_mu  = [1e-4, 5e-4, 1e-3, 5e-3]
list_σ   = collect(range(0.001, 0.1, length=16))
list_rec = vcat([0.0, 0.0001], collect(range(0.0, 0.01, length=16))[2:end])

N   = list_N[i_N]
mu  = list_mu[i_mu]
σ   = list_σ[i_σ]
rec = list_rec[i_rec]
L   = 50
rho = 0.1

σ_str   = @sprintf("%.5f", σ)
rec_str = @sprintf("%.5f", rec)

@printf("N=%d  mu=%.5f  σ=%.5f  rec=%.5f\n", N, mu, σ, rec)

# ── Paths ─────────────────────────────────────────────────────────────────────
N_mu_dir  = @sprintf("N-%d_mu-%.5f", N, mu)
param_dir = @sprintf("epistasis-strength-%s_recombination_strengh-%s_density-%.2f",
                     σ_str, rec_str, rho)

coupling_dir = joinpath(wf_base_dir,  "coupling_localfield")
inf_dir      = joinpath(inf_base_dir, N_mu_dir, param_dir)
out_dir      = joinpath(acc_base_dir, N_mu_dir)
mkpath(out_dir)
out_csv      = joinpath(out_dir, "accuracy_summary_" * param_dir * ".csv")

# ── Ground truth ──────────────────────────────────────────────────────────────
J_GT = readdlm(joinpath(coupling_dir,
    @sprintf("coupling_epistasis-strength-%s_recombination_strengh-%s.txt", σ_str, rec_str)))
S_GT = vec(readdlm(joinpath(coupling_dir,
    @sprintf("local-field_epistasis-strength-%s_recombination_strengh-%s.txt", σ_str, rec_str))))

# ── Constants ─────────────────────────────────────────────────────────────────
const id_ensemble_max        = 100  # attempt to load all available replicates up to 100
const time_estimate_epistasis = [30, 60, 90, 150, 210, 300, 450, 600, 810, 990]
const α_pair  = 0.0001   # pairwise: x_i * x_j > 0.0001 (reduced threshold)

# Upper-triangle mask (j > i), used to select unique pairs
const upper_tri_mask = (repeat(collect(1:L), 1, L) .- repeat(collect(1:L)', L, 1)) .< 0

# ── Metric helpers ────────────────────────────────────────────────────────────
function compute_auroc(scores::Vector, bool_GT::AbstractVector{Bool})
    any(bool_GT) && any(.!bool_GT) || return NaN
    idx_sort = sortperm(scores, rev=true)
    TPv, FPv = get_TP_FP(bool_GT, idx_sort)
    return get_ROC_AUC(TPv, FPv)
end

function compute_metrics(scores::Vector, GT::Vector, bool_GT::AbstractVector{Bool}, σ_scale::Float64)
    auroc = compute_auroc(scores, bool_GT)
    nm    = norm(scores .- σ_scale .* GT)
    cr    = my_cor(scores, σ_scale .* GT)
    return (auroc, nm, cr)
end

# ── Output DataFrame ──────────────────────────────────────────────────────────
results = DataFrame(
    Category = String[], Metric = String[], Time = Int[],
    Mean = Float64[], Error_L = Float64[], Error_H = Float64[])

# ── Main loop ─────────────────────────────────────────────────────────────────
function run_accuracy()
    for t_temp in time_estimate_epistasis

        # ── Load all replicates for this time point ───────────────────────────
        loaded = Dict{Int, NamedTuple}()

        for id_ens in 1:id_ensemble_max
            fkey = @sprintf("_sigma-%s_rec-%s_id-%d_time-%d_density-%.2f",
                            σ_str, rec_str, id_ens, t_temp, rho)

            f_det_a = find_output_file(joinpath(inf_dir, "detectability", "detectability-accum" * fkey))
            f_det_t = find_output_file(joinpath(inf_dir, "detectability", "detectability-temp"  * fkey))
            f_LD_a  = find_output_file(joinpath(inf_dir, "LD",    "LD11-accum"    * fkey))
            f_LD_t  = find_output_file(joinpath(inf_dir, "LD",    "LD11-temp"     * fkey))
            f_KNS_a = find_output_file(joinpath(inf_dir, "KNS",   "KNS-accum"     * fkey))
            f_GC_a  = find_output_file(joinpath(inf_dir, "KNSGC", "KNSGC-accum"   * fkey))
            f_UFE_a = find_output_file(joinpath(inf_dir, "UFE",   "UFE-accum"     * fkey))
            f_UFE_t = find_output_file(joinpath(inf_dir, "UFE",   "UFE-temp"      * fkey))
            f_MPL   = find_output_file(joinpath(inf_dir, "MPL",   "epis_MPL"      * fkey))
            f_SL    = find_output_file(joinpath(inf_dir, "SL",    "epis_SL"       * fkey))
            f_MPLs  = find_output_file(joinpath(inf_dir, "MPL",   "additive_MPL"  * fkey))
            f_SLs   = find_output_file(joinpath(inf_dir, "SL",    "additive_SL"   * fkey))

            flist = [f_det_a, f_det_t, f_LD_a, f_LD_t, f_KNS_a, f_GC_a,
                     f_UFE_a, f_UFE_t, f_MPL, f_SL, f_MPLs, f_SLs]
            if any(isnothing, flist)
                @warn @sprintf("Files missing for id-%d time-%d", id_ens, t_temp)
                continue
            end

            det_a = readdlm_file(f_det_a)
            det_t = readdlm_file(f_det_t)

            idx_a = findall((det_a .> α_pair) .* upper_tri_mask)
            idx_t = findall((det_t .> α_pair) .* upper_tri_mask)

            idx_s = collect(1:L)   # all loci — no frequency threshold for selection

            loaded[id_ens] = (
                LD_a  = vec(readdlm_file(f_LD_a)[idx_a]),
                LD_t  = vec(readdlm_file(f_LD_t)[idx_t]),
                KNS_a = vec(readdlm_file(f_KNS_a)[idx_a]),
                GC_a  = vec(readdlm_file(f_GC_a)[idx_a]),
                UFE_a = vec(readdlm_file(f_UFE_a)[idx_a]),
                UFE_t = vec(readdlm_file(f_UFE_t)[idx_t]),
                MPL   = vec(readdlm_file(f_MPL)[idx_a]),
                SL    = vec(readdlm_file(f_SL)[idx_a]),
                S_MPL = vec(readdlm_file(f_MPLs))[idx_s],
                S_SL  = vec(readdlm_file(f_SLs))[idx_s],
                GT_J_a = vec(J_GT[idx_a]),
                GT_J_t = vec(J_GT[idx_t]),
                GT_S   = S_GT[idx_s],
            )
        end

        avail = sort(collect(keys(loaded)))
        length(avail) < 2 && continue

        # ── True LOO: held-out replicate cycles over available replicates ─────
        auroc_LD_a  = Float64[]; auroc_LD_t  = Float64[]
        auroc_KNS_a = Float64[]
        auroc_GC_a  = Float64[]
        auroc_UFE_a = Float64[]; auroc_UFE_t = Float64[]
        auroc_MPL   = Float64[]; auroc_SL    = Float64[]
        auroc_S_MPL = Float64[]; auroc_S_SL  = Float64[]

        norm_LD_a  = Float64[]; norm_LD_t  = Float64[]
        norm_KNS_a = Float64[]
        norm_GC_a  = Float64[]
        norm_UFE_a = Float64[]; norm_UFE_t = Float64[]
        norm_MPL   = Float64[]; norm_SL    = Float64[]
        norm_S_MPL = Float64[]; norm_S_SL  = Float64[]

        cor_LD_a  = Float64[]; cor_LD_t  = Float64[]
        cor_KNS_a = Float64[]
        cor_GC_a  = Float64[]
        cor_UFE_a = Float64[]; cor_UFE_t = Float64[]
        cor_MPL   = Float64[]; cor_SL    = Float64[]
        cor_S_MPL = Float64[]; cor_S_SL  = Float64[]

        for held_out in avail
            incl = [id for id in avail if id != held_out]
            isempty(incl) && continue

            p(f) = vcat([getproperty(loaded[id], f) for id in incl]...)

            GT_J_a = p(:GT_J_a);  GT_J_t = p(:GT_J_t);  GT_S = p(:GT_S)
            bool_J_a = GT_J_a .> 1e-10
            bool_J_t = GT_J_t .> 1e-10
            bool_S   = GT_S   .> 1e-10

            for (al, nl, cl, scores, GT, bgt) in [
                (auroc_LD_a,  norm_LD_a,  cor_LD_a,  p(:LD_a),  GT_J_a, bool_J_a),
                (auroc_LD_t,  norm_LD_t,  cor_LD_t,  p(:LD_t),  GT_J_t, bool_J_t),
                (auroc_KNS_a, norm_KNS_a, cor_KNS_a, p(:KNS_a), GT_J_a, bool_J_a),
                (auroc_GC_a,  norm_GC_a,  cor_GC_a,  p(:GC_a),  GT_J_a, bool_J_a),
                (auroc_UFE_a, norm_UFE_a, cor_UFE_a, p(:UFE_a), GT_J_a, bool_J_a),
                (auroc_UFE_t, norm_UFE_t, cor_UFE_t, p(:UFE_t), GT_J_t, bool_J_t),
                (auroc_MPL,   norm_MPL,   cor_MPL,   p(:MPL),   GT_J_a, bool_J_a),
                (auroc_SL,    norm_SL,    cor_SL,    p(:SL),    GT_J_a, bool_J_a),
            ]
                a, n, c = compute_metrics(scores, GT, bgt, σ)
                push!(al, a); push!(nl, n); push!(cl, c)
            end

            a, n, c = compute_metrics(p(:S_MPL), GT_S, bool_S, σ)
            push!(auroc_S_MPL, a); push!(norm_S_MPL, n); push!(cor_S_MPL, c)
            a, n, c = compute_metrics(p(:S_SL),  GT_S, bool_S, σ)
            push!(auroc_S_SL,  a); push!(norm_S_SL,  n); push!(cor_S_SL,  c)
        end

        categories = [
            ("LD_Accum",     "auroc",       auroc_LD_a),
            ("KNS_Accum",    "auroc",       auroc_KNS_a),
            ("KNSGC_Accum",  "auroc",       auroc_GC_a),
            ("UFE_Accum",    "auroc",       auroc_UFE_a),
            ("LD_Temp",      "auroc",       auroc_LD_t),
            ("UFE_Temp",     "auroc",       auroc_UFE_t),
            ("MPL",          "auroc",       auroc_MPL),
            ("SL",           "auroc",       auroc_SL),
            ("Selection_MPL","auroc",       auroc_S_MPL),
            ("Selection_SL", "auroc",       auroc_S_SL),
            ("LD_Accum",     "norm",        norm_LD_a),
            ("KNS_Accum",    "norm",        norm_KNS_a),
            ("KNSGC_Accum",  "norm",        norm_GC_a),
            ("UFE_Accum",    "norm",        norm_UFE_a),
            ("LD_Temp",      "norm",        norm_LD_t),
            ("UFE_Temp",     "norm",        norm_UFE_t),
            ("MPL",          "norm",        norm_MPL),
            ("SL",           "norm",        norm_SL),
            ("Selection_MPL","norm",        norm_S_MPL),
            ("Selection_SL", "norm",        norm_S_SL),
            ("LD_Accum",     "correlation", cor_LD_a),
            ("KNS_Accum",    "correlation", cor_KNS_a),
            ("KNSGC_Accum",  "correlation", cor_GC_a),
            ("UFE_Accum",    "correlation", cor_UFE_a),
            ("LD_Temp",      "correlation", cor_LD_t),
            ("UFE_Temp",     "correlation", cor_UFE_t),
            ("MPL",          "correlation", cor_MPL),
            ("SL",           "correlation", cor_SL),
            ("Selection_MPL","correlation", cor_S_MPL),
            ("Selection_SL", "correlation", cor_S_SL),
        ]

        for (cat, metric, data) in categories
            isempty(data) && continue
            valid = filter(!isnan, data)
            isempty(valid) && continue
            mean_v, err_L, err_H = mean_errorbars(data)
            push!(results, (cat, metric, t_temp, mean_v, err_L, err_H))
        end

        # ── Per-replicate analysis ─────────────────────────────────────────────
        pr_LD_a  = Float64[]; pr_LD_t  = Float64[]
        pr_KNS_a = Float64[]
        pr_GC_a  = Float64[]
        pr_UFE_a = Float64[]; pr_UFE_t = Float64[]
        pr_MPL   = Float64[]; pr_SL    = Float64[]
        pr_S_MPL = Float64[]; pr_S_SL  = Float64[]

        pr_nm_LD_a  = Float64[]; pr_nm_LD_t  = Float64[]
        pr_nm_KNS_a = Float64[]
        pr_nm_GC_a  = Float64[]
        pr_nm_UFE_a = Float64[]; pr_nm_UFE_t = Float64[]
        pr_nm_MPL   = Float64[]; pr_nm_SL    = Float64[]
        pr_nm_S_MPL = Float64[]; pr_nm_S_SL  = Float64[]

        pr_cr_LD_a  = Float64[]; pr_cr_LD_t  = Float64[]
        pr_cr_KNS_a = Float64[]
        pr_cr_GC_a  = Float64[]
        pr_cr_UFE_a = Float64[]; pr_cr_UFE_t = Float64[]
        pr_cr_MPL   = Float64[]; pr_cr_SL    = Float64[]
        pr_cr_S_MPL = Float64[]; pr_cr_S_SL  = Float64[]

        for id in avail
            d = loaded[id]
            bJa = d.GT_J_a .> 1e-10
            bJt = d.GT_J_t .> 1e-10
            bS  = d.GT_S   .> 1e-10

            for (al, nl, cl, scores, GT, bgt) in [
                (pr_LD_a,  pr_nm_LD_a,  pr_cr_LD_a,  d.LD_a,  d.GT_J_a, bJa),
                (pr_LD_t,  pr_nm_LD_t,  pr_cr_LD_t,  d.LD_t,  d.GT_J_t, bJt),
                (pr_KNS_a, pr_nm_KNS_a, pr_cr_KNS_a, d.KNS_a, d.GT_J_a, bJa),
                (pr_GC_a,  pr_nm_GC_a,  pr_cr_GC_a,  d.GC_a,  d.GT_J_a, bJa),
                (pr_UFE_a, pr_nm_UFE_a, pr_cr_UFE_a, d.UFE_a, d.GT_J_a, bJa),
                (pr_UFE_t, pr_nm_UFE_t, pr_cr_UFE_t, d.UFE_t, d.GT_J_t, bJt),
                (pr_MPL,   pr_nm_MPL,   pr_cr_MPL,   d.MPL,   d.GT_J_a, bJa),
                (pr_SL,    pr_nm_SL,    pr_cr_SL,    d.SL,    d.GT_J_a, bJa),
            ]
                a, n, c = compute_metrics(scores, GT, bgt, σ)
                push!(al, a); push!(nl, n); push!(cl, c)
            end

            a, n, c = compute_metrics(d.S_MPL, d.GT_S, bS, σ)
            push!(pr_S_MPL, a); push!(pr_nm_S_MPL, n); push!(pr_cr_S_MPL, c)
            a, n, c = compute_metrics(d.S_SL,  d.GT_S, bS, σ)
            push!(pr_S_SL,  a); push!(pr_nm_S_SL,  n); push!(pr_cr_S_SL,  c)
        end

        pr_categories = [
            ("LD_Accum_perRep",     "auroc",       pr_LD_a),
            ("KNS_Accum_perRep",    "auroc",       pr_KNS_a),
            ("KNSGC_Accum_perRep",  "auroc",       pr_GC_a),
            ("UFE_Accum_perRep",    "auroc",       pr_UFE_a),
            ("LD_Temp_perRep",      "auroc",       pr_LD_t),
            ("UFE_Temp_perRep",     "auroc",       pr_UFE_t),
            ("MPL_perRep",          "auroc",       pr_MPL),
            ("SL_perRep",           "auroc",       pr_SL),
            ("Selection_MPL_perRep","auroc",       pr_S_MPL),
            ("Selection_SL_perRep", "auroc",       pr_S_SL),
            ("LD_Accum_perRep",     "norm",        pr_nm_LD_a),
            ("KNS_Accum_perRep",    "norm",        pr_nm_KNS_a),
            ("KNSGC_Accum_perRep",  "norm",        pr_nm_GC_a),
            ("UFE_Accum_perRep",    "norm",        pr_nm_UFE_a),
            ("LD_Temp_perRep",      "norm",        pr_nm_LD_t),
            ("UFE_Temp_perRep",     "norm",        pr_nm_UFE_t),
            ("MPL_perRep",          "norm",        pr_nm_MPL),
            ("SL_perRep",           "norm",        pr_nm_SL),
            ("Selection_MPL_perRep","norm",        pr_nm_S_MPL),
            ("Selection_SL_perRep", "norm",        pr_nm_S_SL),
            ("LD_Accum_perRep",     "correlation", pr_cr_LD_a),
            ("KNS_Accum_perRep",    "correlation", pr_cr_KNS_a),
            ("KNSGC_Accum_perRep",  "correlation", pr_cr_GC_a),
            ("UFE_Accum_perRep",    "correlation", pr_cr_UFE_a),
            ("LD_Temp_perRep",      "correlation", pr_cr_LD_t),
            ("UFE_Temp_perRep",     "correlation", pr_cr_UFE_t),
            ("MPL_perRep",          "correlation", pr_cr_MPL),
            ("SL_perRep",           "correlation", pr_cr_SL),
            ("Selection_MPL_perRep","correlation", pr_cr_S_MPL),
            ("Selection_SL_perRep", "correlation", pr_cr_S_SL),
        ]

        for (cat, metric, data) in pr_categories
            isempty(data) && continue
            valid = filter(!isnan, data)
            isempty(valid) && continue
            mean_v, err_L, err_H = mean_errorbars(data)
            push!(results, (cat, metric, t_temp, mean_v, err_L, err_H))
        end

    end # time loop

    CSV.write(out_csv, results)
    @printf("Written: %s\n", out_csv)
end

run_accuracy()
