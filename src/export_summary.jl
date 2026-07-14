using CSV, DataFrames, DelimitedFiles, Printf, Statistics

# ── Parameter grids ───────────────────────────────────────────────────────────
const list_N   = [100, 500, 1000, 5000, 10000]
const list_mu  = [1e-4, 5e-4, 1e-3, 5e-3]
const list_σ   = collect(range(0.001, 0.1, length=16))
const list_rec = vcat([0.0, 0.0001], collect(range(0.0, 0.01, length=16))[2:end], [0.02, 0.05])

const n_σ_g   = length(list_σ)
const n_rec_g = length(list_rec)
const n_reps  = 30
const rho     = 0.1

const acc_base = "__ACC_DIR__"
const inf_base = "__INF_DIR__"
const out_base = "__OUT_DIR__"

# ── 3×3 box filter over (σ, rec) grid ────────────────────────────────────────
function smooth2d(mat::Matrix{Float64}, w::Int=1)
    nr, nc = size(mat)
    out = fill(NaN, nr, nc)
    for i in 1:nr, j in 1:nc
        vals = Float64[]
        for di in -w:w, dj in -w:w
            ni, nj = i+di, j+dj
            if 1 <= ni <= nr && 1 <= nj <= nc && !isnan(mat[ni,nj])
                push!(vals, mat[ni,nj])
            end
        end
        isempty(vals) || (out[i,j] = mean(vals))
    end
    out
end

# ── Dict accumulator helper ───────────────────────────────────────────────────
function push_d!(d::Dict{Int, Vector{Float64}}, t::Int, v::Float64)
    haskey(d, t) || (d[t] = Float64[])
    push!(d[t], v)
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 1 — Accuracy (perRep + 2D-smoothed at all time points)
# Output: out_base/N-{N}_mu-{mu}/accuracy_perRep_v4.csv
# Reads from: accuracy_output_expanded_v5_Jun2026 (α_locus=0 for selection)
# Columns: N, mu, sigma, rec, method, metric, time, Mean, Error_L, Error_H
#
# Raw rows: all _perRep categories, all 10 time points.
# Smoothed rows: 3×3 box filter over (σ,rec) at every time point;
#   method name suffixed with "_smoothed".
# ══════════════════════════════════════════════════════════════════════════════

function export_accuracy_Nmu(N::Int, mu::Float64)
    N_mu_dir = @sprintf("N-%d_mu-%.5f", N, mu)
    out_dir  = joinpath(out_base, N_mu_dir)
    mkpath(out_dir)

    rows = DataFrame(
        N=Int[], mu=Float64[], sigma=Float64[], rec=Float64[],
        method=String[], metric=String[], time=Int[],
        Mean=Float64[], Error_L=Float64[], Error_H=Float64[])

    n_missing = 0
    for σ in list_σ, rec in list_rec
        σ_str   = @sprintf("%.5f", σ)
        rec_str = @sprintf("%.5f", rec)
        param_dir = @sprintf(
            "epistasis-strength-%s_recombination_strengh-%s_density-%.2f",
            σ_str, rec_str, rho)
        fpath = joinpath(acc_base, N_mu_dir, "accuracy_summary_" * param_dir * ".csv")
        if !isfile(fpath)
            n_missing += 1
            continue
        end
        df     = CSV.read(fpath, DataFrame)
        df_rep = df[endswith.(string.(df.Category), "_perRep"), :]
        for row in eachrow(df_rep)
            method = replace(string(row.Category), "_perRep" => "")
            push!(rows, (N, mu, σ, rec, method, string(row.Metric),
                         row.Time, row.Mean, row.Error_L, row.Error_H))
        end
    end

    # ── 2D smoothing at every time point ─────────────────────────────────────
    smooth_rows = DataFrame(
        N=Int[], mu=Float64[], sigma=Float64[], rec=Float64[],
        method=String[], metric=String[], time=Int[],
        Mean=Float64[], Error_L=Float64[], Error_H=Float64[])

    for method in unique(rows.method)
        sub_m = rows[rows.method .== method, :]
        for metric in unique(sub_m.metric)
            sub_mm = sub_m[sub_m.metric .== metric, :]
            for t_val in sort(unique(sub_mm.time))
                sub = sub_mm[sub_mm.time .== t_val, :]
                isempty(sub) && continue

                mat_M = fill(NaN, n_σ_g, n_rec_g)
                mat_L = fill(NaN, n_σ_g, n_rec_g)
                mat_H = fill(NaN, n_σ_g, n_rec_g)
                for row in eachrow(sub)
                    i_σ = findfirst(s -> abs(s - row.sigma) < 1e-9, list_σ)
                    i_r = findfirst(r -> abs(r - row.rec)   < 1e-9, list_rec)
                    (isnothing(i_σ) || isnothing(i_r)) && continue
                    mat_M[i_σ, i_r] = row.Mean
                    mat_L[i_σ, i_r] = row.Error_L
                    mat_H[i_σ, i_r] = row.Error_H
                end

                sm_M = smooth2d(mat_M)
                sm_L = smooth2d(mat_L)
                sm_H = smooth2d(mat_H)

                for i_σ in 1:n_σ_g, i_r in 1:n_rec_g
                    isnan(sm_M[i_σ, i_r]) && continue
                    push!(smooth_rows, (N, mu, list_σ[i_σ], list_rec[i_r],
                                        method * "_smoothed", metric, t_val,
                                        sm_M[i_σ, i_r], sm_L[i_σ, i_r], sm_H[i_σ, i_r]))
                end
            end
        end
    end

    all_rows = vcat(rows, smooth_rows)
    CSV.write(joinpath(out_dir, "accuracy_perRep.csv"), all_rows)
    @printf("  N=%5d mu=%.1e: %d raw + %d smoothed rows  (%d files missing)\n",
            N, mu, nrow(rows), nrow(smooth_rows), n_missing)
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 2 — Physics metrics
# Output: out_base/N-{N}_mu-{mu}/physics-metrics.csv
# (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

function export_physics_Nmu(N::Int, mu::Float64)
    N_mu_dir = @sprintf("N-%d_mu-%.5f", N, mu)
    out_dir  = joinpath(out_base, N_mu_dir)
    mkpath(out_dir)

    rows = DataFrame(
        N=Int[], mu=Float64[], sigma=Float64[], rec=Float64[],
        time=Int[], metric=String[], mean=Float64[], std=Float64[])

    n_combos  = 0
    n_missing = 0

    for σ in list_σ, rec in list_rec
        σ_str   = @sprintf("%.5f", σ)
        rec_str = @sprintf("%.5f", rec)
        param_dir = @sprintf(
            "epistasis-strength-%s_recombination_strengh-%s_density-%.2f",
            σ_str, rec_str, rho)
        inf_dir = joinpath(inf_base, N_mu_dir, param_dir)
        if !isdir(inf_dir)
            n_missing += 1
            continue
        end

        het_d    = Dict{Int, Vector{Float64}}()
        hetmax_d = Dict{Int, Vector{Float64}}()
        npoly_d  = Dict{Int, Vector{Float64}}()
        nhigh_d  = Dict{Int, Vector{Float64}}()
        nlow_d   = Dict{Int, Vector{Float64}}()
        fmean_d  = Dict{Int, Vector{Float64}}()
        fvar_d   = Dict{Int, Vector{Float64}}()
        pr_d     = Dict{Int, Vector{Float64}}()
        qlen_d   = Dict{Int, Vector{Float64}}()
        qleD_d   = Dict{Int, Vector{Float64}}()

        for id in 1:n_reps
            key = @sprintf("_sigma-%s_rec-%s_id-%d", σ_str, rec_str, id)

            f = joinpath(inf_dir, "diversity_time_course", "diversity" * key * ".txt")
            if isfile(f)
                d = readdlm(f)
                for ri in 1:size(d, 1)
                    t = round(Int, d[ri, 1])
                    t % 30 == 0 || continue
                    push_d!(het_d,    t, Float64(d[ri, 2]))
                    push_d!(hetmax_d, t, Float64(d[ri, 7]))
                    push_d!(npoly_d,  t, Float64(d[ri, 4]))
                    push_d!(nhigh_d,  t, Float64(d[ri, 5]))
                    push_d!(nlow_d,   t, Float64(d[ri, 6]))
                end
            end

            f = joinpath(inf_dir, "fitness_wave", "Fitness_evo" * key * ".txt")
            if isfile(f)
                d = readdlm(f)
                for ri in 1:size(d, 1)
                    t = round(Int, d[ri, 1])
                    t % 30 == 0 || continue
                    push_d!(fmean_d, t, Float64(d[ri, 2]))
                    push_d!(fvar_d,  t, Float64(d[ri, 6]))
                end
            end

            f = joinpath(inf_dir, "participation_rate", "participation_rate" * key * ".txt")
            if isfile(f)
                d = readdlm(f)
                for ri in 1:size(d, 1)
                    t = round(Int, d[ri, 1])
                    t % 30 == 0 || continue
                    push_d!(pr_d, t, Float64(d[ri, 2]))
                end
            end

            f = joinpath(inf_dir, "QLE_order_parameter", "QLE_order_parameter" * key * ".txt")
            if isfile(f)
                d = readdlm(f)
                for ri in 1:size(d, 1)
                    t = round(Int, d[ri, 1])
                    t % 30 == 0 || continue
                    push_d!(qlen_d, t, Float64(d[ri, 2]))
                    push_d!(qleD_d, t, Float64(d[ri, 3]))
                end
            end
        end

        function push_metric_d!(dct, mname)
            for t in sort(collect(keys(dct)))
                vs = dct[t]
                isempty(vs) && continue
                push!(rows, (N, mu, σ, rec, t, mname,
                             mean(vs), length(vs) > 1 ? std(vs) : 0.0))
            end
        end
        push_metric_d!(het_d,    "heterozygosity")
        push_metric_d!(hetmax_d, "heterozygosity_max")
        push_metric_d!(npoly_d,  "n_poly")
        push_metric_d!(nhigh_d,  "n_high")
        push_metric_d!(nlow_d,   "n_low")
        push_metric_d!(fmean_d,  "fitness_mean")
        push_metric_d!(fvar_d,   "fitness_variance")
        push_metric_d!(pr_d,     "participation_rate")
        push_metric_d!(qlen_d,   "QLE_numer")
        push_metric_d!(qleD_d,   "QLE_D")

        f_qle2 = joinpath(inf_dir, "QLE_ensemble_v2",
                          @sprintf("QLE_ensemble_v2_sigma-%s_rec-%s.txt", σ_str, rec_str))
        if isfile(f_qle2)
            d = readdlm(f_qle2)
            for ri in 1:size(d, 1)
                t_qle  = round(Int, d[ri, 1])
                t_qle % 30 == 0 || continue
                q2_val = Float64(d[ri, 5])
                q1q2_v = q2_val > 0.0 ? Float64(d[ri, 4]) / q2_val : NaN
                push!(rows, (N, mu, σ, rec, t_qle, "QLE_ensemble_Q",          Float64(d[ri, 2]), NaN))
                push!(rows, (N, mu, σ, rec, t_qle, "QLE_ensemble_D",          Float64(d[ri, 3]), NaN))
                push!(rows, (N, mu, σ, rec, t_qle, "QLE_ensemble_q1",         Float64(d[ri, 4]), NaN))
                push!(rows, (N, mu, σ, rec, t_qle, "QLE_ensemble_q2",         q2_val,            NaN))
                push!(rows, (N, mu, σ, rec, t_qle, "QLE_ensemble_q3",         Float64(d[ri, 6]), NaN))
                push!(rows, (N, mu, σ, rec, t_qle, "QLE_ensemble_q1_over_q2", q1q2_v,            NaN))
            end
        else
            for mname in ["QLE_ensemble_Q", "QLE_ensemble_D", "QLE_ensemble_q1",
                          "QLE_ensemble_q2", "QLE_ensemble_q3", "QLE_ensemble_q1_over_q2"]
                push!(rows, (N, mu, σ, rec, 990, mname, NaN, NaN))
            end
        end

        n_combos += 1
    end

    smooth_rows = DataFrame(
        N=Int[], mu=Float64[], sigma=Float64[], rec=Float64[],
        time=Int[], metric=String[], mean=Float64[], std=Float64[])

    for mname in unique(rows.metric)
        sub_m = rows[rows.metric .== mname, :]
        for t_val in sort(unique(sub_m.time))
            sub = sub_m[sub_m.time .== t_val, :]
            isempty(sub) && continue

            mat_mean = fill(NaN, n_σ_g, n_rec_g)
            mat_std  = fill(NaN, n_σ_g, n_rec_g)
            for row in eachrow(sub)
                i_σ = findfirst(s -> abs(s - row.sigma) < 1e-9, list_σ)
                i_r = findfirst(r -> abs(r - row.rec)   < 1e-9, list_rec)
                (isnothing(i_σ) || isnothing(i_r)) && continue
                mat_mean[i_σ, i_r] = row.mean
                mat_std[i_σ, i_r]  = row.std
            end

            sm_mean = smooth2d(mat_mean)
            sm_std  = smooth2d(mat_std)

            for i_σ in 1:n_σ_g, i_r in 1:n_rec_g
                isnan(sm_mean[i_σ, i_r]) && continue
                push!(smooth_rows, (N, mu, list_σ[i_σ], list_rec[i_r],
                                    t_val, mname * "_smoothed",
                                    sm_mean[i_σ, i_r], sm_std[i_σ, i_r]))
            end
        end
    end

    all_rows = vcat(rows, smooth_rows)
    CSV.write(joinpath(out_dir, "physics-metrics.csv"), all_rows)
    @printf("  N=%5d mu=%.1e: %d raw + %d smoothed rows, %d combos  (%d missing)\n",
            N, mu, nrow(rows), nrow(smooth_rows), n_combos, n_missing)
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
mkpath(out_base)
@printf("Output directory: %s\n", out_base)
@printf("Grid: %d N × %d μ × %d σ × %d rec\n\n",
        length(list_N), length(list_mu), length(list_σ), length(list_rec))

run_N  = length(ARGS) >= 1 ? [list_N[parse(Int, ARGS[1])]]  : list_N
run_mu = length(ARGS) >= 2 ? [list_mu[parse(Int, ARGS[2])]] : list_mu

@printf("=== Part 1: Accuracy (perRep + smoothed, α_locus=0) ===\n")
for N in run_N, mu in run_mu
    export_accuracy_Nmu(N, mu)
end

@printf("\n=== Part 2: Physics metrics (all t%%30==0) ===\n")
for N in run_N, mu in run_mu
    export_physics_Nmu(N, mu)
    GC.gc()
end

@printf("\nDone. Files written to %s\n", out_base)
