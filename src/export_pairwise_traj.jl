using DelimitedFiles, CSV, DataFrames, Printf, CodecZlib

# Exports per-replicate pairwise co-frequency trajectories x2[i,j].
# Covers diagonal (i==j) and upper triangle (i<j).
# No averaging — one row per (i, j, replicate id, time).
#
# Output columns: i, j, id, time, x2, GT
#   diagonal (i==j): x2 = allele frequency x1[i],        GT = S_GT[i]
#   off-diagonal (i<j): x2 = double-mutant co-frequency,  GT = J_GT[i,j]
#
# Usage: julia export_pairwise_freq_traj_v4_Jun2026.jl [sigma_str rec_str]
# Default: sigma=0.10000, rec=0.00000

const N      = 1000
const mu     = 0.001
const L      = 50
const rho    = 0.1
const n_reps = 30
const time_checkpoints = [30, 60, 90, 150, 210, 300, 450, 600, 810, 990]

const inf_base = "__INF_DIR__"
const wf_base  = "__WF_DIR__"
const out_base = "__OUT_DIR__"

function read_x2(path)
    open(GzipDecompressorStream, path) do io
        readdlm(io)
    end
end

function build_df(per_rep, J_GT, S_GT)
    n_rows = (L * (L + 1) ÷ 2) * n_reps * length(time_checkpoints)
    col_i    = Vector{Int}(undef, n_rows)
    col_j    = Vector{Int}(undef, n_rows)
    col_id   = Vector{Int}(undef, n_rows)
    col_time = Vector{Int}(undef, n_rows)
    col_x2   = Vector{Float64}(undef, n_rows)
    col_gt   = Vector{Float64}(undef, n_rows)

    row = 0
    for t in time_checkpoints
        reps = per_rep[t]
        for id in sort(collect(keys(reps)))
            mat = reps[id]
            for i in 1:L
                # diagonal: x2[i,i] = x1[i],  GT = S_GT[i]
                row += 1
                col_i[row]    = i
                col_j[row]    = i
                col_id[row]   = id
                col_time[row] = t
                col_x2[row]   = mat[i, i]
                col_gt[row]   = S_GT[i]

                # upper triangle: x2[i,j],  GT = J_GT[i,j]
                for j in (i+1):L
                    row += 1
                    col_i[row]    = i
                    col_j[row]    = j
                    col_id[row]   = id
                    col_time[row] = t
                    col_x2[row]   = mat[i, j]
                    col_gt[row]   = J_GT[i, j]
                end
            end
        end
    end

    DataFrame(
        i    = col_i[1:row],
        j    = col_j[1:row],
        id   = col_id[1:row],
        time = col_time[1:row],
        x2   = col_x2[1:row],
        GT   = col_gt[1:row])
end

function main()
    σ_str   = length(ARGS) >= 1 ? ARGS[1] : "0.10000"
    rec_str = length(ARGS) >= 2 ? ARGS[2] : "0.00000"

    N_mu_dir  = @sprintf("N-%d_mu-%.5f", N, mu)
    param_dir = @sprintf("epistasis-strength-%s_recombination_strengh-%s_density-%.2f",
                         σ_str, rec_str, rho)
    inf_dir      = joinpath(inf_base, N_mu_dir, param_dir)
    coupling_dir = joinpath(wf_base, "coupling_localfield")
    out_dir      = joinpath(out_base, N_mu_dir)
    mkpath(out_dir)

    @printf("N=%d  mu=%.4f  σ=%s  rec=%s\n", N, mu, σ_str, rec_str)

    # Ground truth
    J_GT = readdlm(joinpath(coupling_dir,
        @sprintf("coupling_epistasis-strength-%s_recombination_strengh-%s.txt", σ_str, rec_str)))
    S_GT = vec(readdlm(joinpath(coupling_dir,
        @sprintf("local-field_epistasis-strength-%s_recombination_strengh-%s.txt", σ_str, rec_str))))

    # Load x2-temp matrices
    per_rep = Dict{Int, Dict{Int, Matrix{Float64}}}()
    for t in time_checkpoints
        per_rep[t] = Dict{Int, Matrix{Float64}}()
        n_loaded = 0
        for id in 1:n_reps
            fpath = joinpath(inf_dir, "x2",
                @sprintf("x2-temp_sigma-%s_rec-%s_id-%d_time-%d_density-%.2f.txt.gz",
                         σ_str, rec_str, id, t, rho))
            if isfile(fpath)
                per_rep[t][id] = read_x2(fpath)
                n_loaded += 1
            else
                @warn @sprintf("Missing: id=%d t=%d", id, t)
            end
        end
        @printf("  t=%4d: loaded %d/%d replicates\n", t, n_loaded, n_reps)
    end

    df = build_df(per_rep, J_GT, S_GT)

    out_csv = joinpath(out_dir,
        @sprintf("pairwise_freq_traj_v4_sigma-%s_rec-%s.csv", σ_str, rec_str))
    CSV.write(out_csv, df)
    @printf("\nWritten %d rows  →  %s\n", nrow(df), out_csv)
    @printf("  diagonal (i==j): %d rows per (rep, time)  x2=x1[i]    GT=S_GT[i]\n", L)
    @printf("  off-diag (i< j): %d rows per (rep, time)  x2=x2[i,j]  GT=J_GT[i,j]\n", L*(L-1)÷2)
end

main()
