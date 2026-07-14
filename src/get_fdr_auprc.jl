using Random, Statistics, LinearAlgebra, Printf, CSV, DataFrames, DelimitedFiles, CodecZlib

include("epistasis_detectability.jl")

# ── I/O helpers ───────────────────────────────────────────────────────────────
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
# Usage: julia get_fdr_auprc_expanded_grid_Jun2026_v3.jl \
#            i_N i_mu i_σ i_rec wf_base_dir inf_base_dir fdr_base_dir
#
# Positive epistasis + selection, perRep only.
# Metric: tpr_at5pct_fpr (TPR at 5% FPR = statistical power).
# Selection uses all L loci (no detectability threshold).
# list_rec: 17-value base grid.
# Output: fdr_auprc_output_expanded_v3_Jun2026/
i_N          = parse(Int, ARGS[1])
i_mu         = parse(Int, ARGS[2])
i_σ          = parse(Int, ARGS[3])
i_rec        = parse(Int, ARGS[4])
wf_base_dir  = ARGS[5]
inf_base_dir = ARGS[6]
fdr_base_dir = ARGS[7]

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
out_dir      = joinpath(fdr_base_dir, N_mu_dir)
mkpath(out_dir)
out_csv      = joinpath(out_dir, "fdr_auprc_summary_" * param_dir * ".csv")

# ── Ground truth ──────────────────────────────────────────────────────────────
J_GT = readdlm(joinpath(coupling_dir,
    @sprintf("coupling_epistasis-strength-%s_recombination_strengh-%s.txt", σ_str, rec_str)))
S_GT = vec(readdlm(joinpath(coupling_dir,
    @sprintf("local-field_epistasis-strength-%s_recombination_strengh-%s.txt", σ_str, rec_str))))

# ── Constants ─────────────────────────────────────────────────────────────────
const id_ensemble_max         = 100
const time_estimate_epistasis = [30, 60, 90, 150, 210, 300, 450, 600, 810, 990]
const α_pair = 0.0001

const upper_tri_mask = (repeat(collect(1:L), 1, L) .- repeat(collect(1:L)', L, 1)) .< 0

# ── Metric helper ─────────────────────────────────────────────────────────────
function compute_tpr_at_fpr(scores::Vector, bool_GT::AbstractVector{Bool}, target_fpr::Float64=0.05)
    any(bool_GT) && any(.!bool_GT) || return NaN
    n_pos = sum(bool_GT)
    n_neg = sum(.!bool_GT)
    tp = 0; fp = 0
    for i in sortperm(scores, rev=true)
        bool_GT[i] ? (tp += 1) : (fp += 1)
        fp / n_neg >= target_fpr && return tp / n_pos
    end
    return tp / n_pos
end

# ── Method tables ─────────────────────────────────────────────────────────────
const METHODS_EPIS_A = [
    ("LD_Accum",    :LD_a),
    ("KNS_Accum",   :KNS_a),
    ("KNSGC_Accum", :GC_a),
    ("UFE_Accum",   :UFE_a),
    ("MPL",         :MPL),
    ("SL",          :SL),
]
const METHODS_EPIS_T = [
    ("LD_Temp",  :LD_t),
    ("UFE_Temp", :UFE_t),
]
const METHODS_SEL = [
    ("Selection_MPL", :S_MPL),
    ("Selection_SL",  :S_SL),
]

const METRIC_NAMES = ["tpr_at5pct_fpr"]

# ── Output DataFrame ──────────────────────────────────────────────────────────
results = DataFrame(
    Category = String[], Metric = String[], Time = Int[],
    Mean = Float64[], Error_L = Float64[], Error_H = Float64[])

# ── Main loop ─────────────────────────────────────────────────────────────────
function run_fdr_auprc()
    for t_temp in time_estimate_epistasis

        loaded = Dict{Int, NamedTuple}()

        for id_ens in 1:id_ensemble_max
            fkey = @sprintf("_sigma-%s_rec-%s_id-%d_time-%d_density-%.2f",
                            σ_str, rec_str, id_ens, t_temp, rho)

            f_det_a = find_output_file(joinpath(inf_dir, "detectability", "detectability-accum" * fkey))
            f_det_t = find_output_file(joinpath(inf_dir, "detectability", "detectability-temp"  * fkey))
            f_LD_a  = find_output_file(joinpath(inf_dir, "LD",    "LD11-accum"   * fkey))
            f_LD_t  = find_output_file(joinpath(inf_dir, "LD",    "LD11-temp"    * fkey))
            f_KNS_a = find_output_file(joinpath(inf_dir, "KNS",   "KNS-accum"    * fkey))
            f_GC_a  = find_output_file(joinpath(inf_dir, "KNSGC", "KNSGC-accum"  * fkey))
            f_UFE_a = find_output_file(joinpath(inf_dir, "UFE",   "UFE-accum"    * fkey))
            f_UFE_t = find_output_file(joinpath(inf_dir, "UFE",   "UFE-temp"     * fkey))
            f_MPL   = find_output_file(joinpath(inf_dir, "MPL",   "epis_MPL"     * fkey))
            f_SL    = find_output_file(joinpath(inf_dir, "SL",    "epis_SL"      * fkey))
            f_MPLs  = find_output_file(joinpath(inf_dir, "MPL",   "additive_MPL" * fkey))
            f_SLs   = find_output_file(joinpath(inf_dir, "SL",    "additive_SL"  * fkey))

            flist = [f_det_a, f_det_t, f_LD_a, f_LD_t, f_KNS_a, f_GC_a,
                     f_UFE_a, f_UFE_t, f_MPL, f_SL, f_MPLs, f_SLs]
            any(isnothing, flist) && continue

            det_a = readdlm_file(f_det_a)
            det_t = readdlm_file(f_det_t)
            idx_a = findall((det_a .> α_pair) .* upper_tri_mask)
            idx_t = findall((det_t .> α_pair) .* upper_tri_mask)
            idx_s = collect(1:L)

            loaded[id_ens] = (
                LD_a   = vec(readdlm_file(f_LD_a)[idx_a]),
                LD_t   = vec(readdlm_file(f_LD_t)[idx_t]),
                KNS_a  = vec(readdlm_file(f_KNS_a)[idx_a]),
                GC_a   = vec(readdlm_file(f_GC_a)[idx_a]),
                UFE_a  = vec(readdlm_file(f_UFE_a)[idx_a]),
                UFE_t  = vec(readdlm_file(f_UFE_t)[idx_t]),
                MPL    = vec(readdlm_file(f_MPL)[idx_a]),
                SL     = vec(readdlm_file(f_SL)[idx_a]),
                S_MPL  = vec(readdlm_file(f_MPLs))[idx_s],
                S_SL   = vec(readdlm_file(f_SLs))[idx_s],
                GT_J_a = vec(J_GT[idx_a]),
                GT_J_t = vec(J_GT[idx_t]),
                GT_S   = S_GT[idx_s],
            )
        end

        avail = sort(collect(keys(loaded)))
        isempty(avail) && continue

        # ── Per-replicate metrics ──────────────────────────────────────────────
        pr_data = Dict{Tuple{String,String}, Vector{Float64}}()
        for (cat, _) in vcat(METHODS_EPIS_A, METHODS_EPIS_T, METHODS_SEL)
            for m in METRIC_NAMES
                pr_data[(cat, m)] = Float64[]
            end
        end

        for id in avail
            d   = loaded[id]
            bJa = d.GT_J_a .> 1e-10
            bJt = d.GT_J_t .> 1e-10
            bS  = d.GT_S   .> 1e-10

            for (methods, bgt_sym) in [
                    (METHODS_EPIS_A, :bJa),
                    (METHODS_EPIS_T, :bJt),
                    (METHODS_SEL,    :bS)]
                bgt = bgt_sym === :bJa ? bJa : bgt_sym === :bJt ? bJt : bS
                for (cat, sym) in methods
                    scores = getproperty(d, sym)
                    push!(pr_data[(cat, "tpr_at5pct_fpr")], compute_tpr_at_fpr(scores, bgt))
                end
            end
        end

        for (cat, _) in vcat(METHODS_EPIS_A, METHODS_EPIS_T, METHODS_SEL)
            for m in METRIC_NAMES
                data  = pr_data[(cat, m)]
                isempty(data) && continue
                valid = filter(!isnan, data)
                isempty(valid) && continue
                mean_v, err_L, err_H = mean_errorbars(data)
                push!(results, (cat, m, t_temp, mean_v, err_L, err_H))
            end
        end

    end  # time loop

    CSV.write(out_csv, results)
    @printf("Written: %s\n", out_csv)
end

run_fdr_auprc()
