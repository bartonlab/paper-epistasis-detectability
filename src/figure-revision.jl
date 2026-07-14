# ── Basic helpers ─────────────────────────────────────────────────────────────
scaled_fontsize(fig_width; scale_factor=20/600) = fig_width * scale_factor

FS_L = 18

# ── Data loading ──────────────────────────────────────────────────────────────
# load_data(N, mu) returns (df_acc, df_phys). Requires data_base to be defined.
function load_data(N::Int, mu::Float64)
    dir = joinpath(data_base, @sprintf("N-%d_mu-%.5f", N, mu))
    df_acc  = CSV.read(joinpath(dir, "accuracy_perRep.csv"),  DataFrame)
    df_phys = CSV.read(joinpath(dir, "physics-metrics.csv"), DataFrame)
    @printf("Loaded N=%d μ=%.1e: %d accuracy rows, %d physics rows\n",
            N, mu, nrow(df_acc), nrow(df_phys))
    @printf("  accuracy methods : %s\n", join(sort(unique(df_acc.method)),  ", "))
    @printf("  accuracy metrics : %s\n", join(sort(unique(df_acc.metric)),  ", "))
    @printf("  accuracy times   : %s\n", join(string.(sort(unique(df_acc.time))), ", "))
    @printf("  physics metrics  : %s\n", join(sort(unique(df_phys.metric)), ", "))
    return df_acc, df_phys
end

# ── get_acc_matrix ────────────────────────────────────────────────────────────
# Returns (mat_mean, mat_errL, mat_errH), each (n_σ × n_rec).
# Requires list_σ, list_rec, n_σ, n_rec to be defined globally.
#
# Usage:
#   mat, _, _ = get_acc_matrix(df_acc, "KNS_Accum_smoothed", "auroc", 990)

function get_acc_matrix(df::DataFrame, method::String, metric::String, time::Int)
    sub = df[(df.method .== method) .&
             (df.metric .== metric) .&
             (df.time   .== time),  :]
    if isempty(sub)
        @warn "No rows: method=$method  metric=$metric  time=$time"
        return fill(NaN, n_σ, n_rec), fill(NaN, n_σ, n_rec), fill(NaN, n_σ, n_rec)
    end
    mat_M = fill(NaN, n_σ, n_rec)
    mat_L = fill(NaN, n_σ, n_rec)
    mat_H = fill(NaN, n_σ, n_rec)
    map_sigma = Dict(zip(list_σ, collect(1:n_σ)))
    map_rec   = Dict(zip(list_rec, collect(1:n_rec)))
    for row in eachrow(sub)
        i_s = map_sigma[row.sigma]
        i_r = map_rec[row.rec]
        mat_M[i_s, i_r] = row.Mean
        mat_L[i_s, i_r] = row.Error_L
        mat_H[i_s, i_r] = row.Error_H
    end
    return mat_M, mat_L, mat_H
end

# ── get_acc_timeseries ────────────────────────────────────────────────────────
# Returns (time_vec, mean_vec, errL_vec, errH_vec) for a fixed
# (sigma, rec, method, metric).  sigma and rec are matched to nearest values.
#
# Usage (by value):
#   get_acc_timeseries(df_acc, 0.047, 0.0, "KNS_Accum_smoothed", "auroc")
# Usage (by index):
#   get_acc_timeseries_idx(df_acc, 8, 1, "KNS_Accum_smoothed", "auroc")

function get_acc_timeseries(df::DataFrame, sigma::Float64, rec::Float64,
                            method::String, metric::String)
    s_v   = list_σ[argmin(abs.(list_σ   .- sigma))]
    rec_v = list_rec[argmin(abs.(list_rec .- rec))]
    sub = df[(abs.(df.sigma .- s_v)   .< 1e-9) .&
             (abs.(df.rec   .- rec_v) .< 1e-9) .&
             (df.method .== method)             .&
             (df.metric .== metric),  :]
    if isempty(sub)
        @warn "No rows: sigma=$s_v  rec=$rec_v  method=$method  metric=$metric"
        return Int[], Float64[], Float64[], Float64[]
    end
    sub = sort(sub, :time)
    return sub.time, sub.Mean, sub.Error_L, sub.Error_H
end

function get_acc_timeseries_idx(df::DataFrame, i_sig::Int, i_rec_idx::Int,
                                method::String, metric::String)
    get_acc_timeseries(df, list_σ[i_sig], list_rec[i_rec_idx], method, metric)
end

# ── get_phys_matrix ───────────────────────────────────────────────────────────
# Returns (mat_mean, mat_std), each (n_σ × n_rec), for a given physics metric.
#
# Usage:
#   mat, mat_std = get_phys_matrix(df_phys, "fitness_variance", 990)

function get_phys_matrix(df::DataFrame, metric::String, time::Int)
    sub = df[(df.metric .== metric) .& (df.time .== time), :]
    if isempty(sub)
        @warn "No rows: metric=$metric"
        return fill(NaN, n_σ, n_rec), fill(NaN, n_σ, n_rec)
    end
    mat_M = fill(NaN, n_σ, n_rec)
    mat_S = fill(NaN, n_σ, n_rec)
    map_sigma = Dict(zip(list_σ, collect(1:n_σ)))
    map_rec   = Dict(zip(list_rec, collect(1:n_rec)))
    for row in eachrow(sub)
        i_s = map_sigma[row.sigma]
        i_r = map_rec[row.rec]
        mat_M[i_s, i_r] = row.mean
        mat_S[i_s, i_r] = row.std
    end
    return mat_M, mat_S
end

# ── get_phys_timeseries ───────────────────────────────────────────────────────
# Returns (time_vec, mean_vec, std_vec) for a fixed (sigma, rec, metric).
# Requires list_σ, list_rec to be defined globally.
#
# Usage:  get_phys_timeseries_idx(df_phys, 8, 1, "fitness_mean")

function get_phys_timeseries(df::DataFrame, sigma::Float64, rec::Float64,
                              metric::String)
    s_v   = list_σ[argmin(abs.(list_σ   .- sigma))]
    rec_v = list_rec[argmin(abs.(list_rec .- rec))]
    sub = df[(abs.(df.sigma .- s_v)   .< 1e-9) .&
             (abs.(df.rec   .- rec_v) .< 1e-9) .&
             (df.metric .== metric),  :]
    if isempty(sub)
        @warn "No rows: sigma=$s_v  rec=$rec_v  metric=$metric"
        return Int[], Float64[], Float64[]
    end
    sub = sort(sub, :time)
    return sub.time, sub.mean, sub.std
end

function get_phys_timeseries_idx(df::DataFrame, i_sig::Int, i_rec::Int,
                                  metric::String)
    get_phys_timeseries(df, list_σ[i_sig], list_rec[i_rec], metric)
end

# ── get_auroc_pt ──────────────────────────────────────────────────────────────
# Extract AUROC (mean, errL, errH) at a single (sigma_idx, rec_idx, method, t)
# from a pre-loaded accuracy DataFrame.
function get_auroc_pt(df::DataFrame, sigma_idx::Int, rec_idx::Int,
                      method::String, t::Int=990)
    σ_v = list_σ[sigma_idx]
    r_v = list_rec[rec_idx]
    sub = df[(df.method .== method) .& (df.metric .== "auroc") .&
             (df.time   .== t) .&
             (abs.(df.sigma .- σ_v) .< 1e-9) .&
             (abs.(df.rec   .- r_v) .< 1e-9), :]
    isempty(sub) && return NaN, NaN, NaN
    return sub.Mean[1], sub.Error_L[1], sub.Error_H[1]
end
