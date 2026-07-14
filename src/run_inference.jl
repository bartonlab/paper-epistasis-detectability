using Statistics
using DelimitedFiles
using Distributions
using Random
using LinearAlgebra
using Printf

# ── Helper functions inlined from systematic_performance_evaluation_Feb21_2026.jl
# ── and function_for_wave.jl ──────────────────────────────────────────────────

stepfunc(x) = x > 0 ? 1 : 0
G(i, j, L)  = stepfunc(j - i) * Int((i - 1) * (L - i / 2.0) + j - i)

function get_sample_at_t(data, t_get)
    n_t = []
    sample_t = []
    for n in 1:size(data, 1)
        if Int(data[n, 1]) == t_get
            if length(sample_t) > 0
                sample_t = hcat(sample_t, data[n, 4:end])
                n_t = vcat(n_t, data[n, 2])
            end
            if length(sample_t) == 0
                sample_t = copy(data[n, 4:end])
                n_t = data[n, 2]
            end
        end
        if data[n, 1] > t_get
            break
        end
    end
    return (n_t, sample_t)
end

function flat_J_to_mat_J(s_epis, L)
    J = zeros(L, L)
    n = 1
    for i in 1:L
        for j in (i+1):L
            J[i, j] = s_epis[L + n]
            J[j, i] = s_epis[L + n]
            n += 1
        end
    end
    return J
end

function get_x_1w2(seq_in, L, LLhalf)
    x1 = seq_in[1:L]
    x2_vec = zeros(LLhalf)
    for i in 1:L
        for j in (i+1):L
            ξ = G(i, j, L)
            x2_vec[ξ] = seq_in[i] * seq_in[j]
        end
    end
    return [x1; x2_vec]
end

function get_d_mu(x_1w2, L, LLhalf)
    d_mu_temp = zeros(L + LLhalf)
    for i in 1:L
        d_mu_temp[i] = 1 - 2 * x_1w2[i]
        for j in (i+1):L
            ξ = G(i, j, L) + L
            d_mu_temp[ξ] = x_1w2[i] + x_1w2[j] - 4 * x_1w2[ξ]
        end
    end
    return d_mu_temp
end

function get_d_rec(x_1w2, L, LLhalf)
    d_rec_temp = zeros(L + LLhalf)
    scale_L = 1.0 / L
    for i in 1:L
        for j in (i+1):L
            ξ = G(i, j, L) + L
            d_rec_temp[ξ] = ((j - i) * scale_L) * (x_1w2[ξ] - x_1w2[i] * x_1w2[j])
        end
    end
    return d_rec_temp
end

function get_x_1w2_set_and_d_mu_random_compression(sample_t, n_t, L, rank_x)
    Neff = Int(sum(n_t))
    LLhalf = Int(L * (L - 1) / 2)
    n_species = size(sample_t, 2)
    x_1w2 = zeros(L + LLhalf)
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        x_1w2_m = get_x_1w2(sample_t[1:L, m], L, LLhalf)
        x_1w2 += n_t_scale * x_1w2_m
    end
    d_mu = get_d_mu(x_1w2, L, LLhalf)
    d_rec = get_d_rec(x_1w2, L, LLhalf)
    x_1w2_set = zeros(n_species, rank_x)
    for m in 1:n_species
        x_1w2_m = sqrt(n_t[m] / Neff) * (get_x_1w2(sample_t[1:L, m], L, LLhalf) - x_1w2)
        x_1w2_set[m, :] = copy(x_1w2_m)
    end
    return (x_1w2_set, x_1w2, d_mu, d_rec)
end

function get_averaged_Z(L, x1, x2, threshold)
    Z_value = zeros(L, L)
    Z_av_value = 0
    n_count = 0
    for i in 1:L
        for j in (i+1):L
            x11 = x2[i, j]
            x01 = x1[j] - x2[i, j]
            x10 = x1[i] - x2[i, j]
            x00 = 1 + x2[i, j] - x1[i] - x1[j]
            if x00 > threshold && x01 > threshold && x10 > threshold && x11 > threshold
                Z_temp = (x11 / x00) / (x01 * x10)
                Z_value[i, j] = Z_temp; Z_value[j, i] = Z_temp
                Z_av_value += Z_temp
                n_count += 1
            end
        end
    end
    if n_count > 0; Z_av_value /= n_count; end
    return (Z_value, Z_av_value)
end

function get_averaged_LD(L, x1, x2, J_GT, threshold)
    LD_av_value = 0
    r2_tot, r2_neu, r2_pos, r2_neg = 0.0, 0.0, 0.0, 0.0
    n_count = 0
    n_neu, n_pos, n_neg = 0, 0, 0
    C = x2 - x1 * x1'
    LD_value = zeros(L, L)
    for i in 1:L
        for j in (i+1):L
            if (threshold < x1[i] < 1 - threshold) && (threshold < x1[j] < 1 - threshold)
                LD_temp = C[i, j] / sqrt(x1[i] * (1 - x1[i]) * x1[j] * (1 - x1[j]))
                LD_value[i, j] = LD_temp; LD_value[j, i] = LD_temp
                LD_av_value += LD_temp
                n_count += 1
                r2_tot += LD_temp^2
                if J_GT[i, j] > 0;  r2_pos += LD_temp^2; n_pos += 1; end
                if J_GT[i, j] == 0; r2_neu += LD_temp^2; n_neu += 1; end
                if J_GT[i, j] < 0;  r2_neg += LD_temp^2; n_neg += 1; end
            end
        end
    end
    if n_count > 0; LD_av_value /= n_count; r2_tot /= n_count; end
    if n_pos > 0; r2_pos /= n_pos; end
    if n_neu > 0; r2_neu /= n_neu; end
    if n_neg > 0; r2_neg /= n_neg; end
    return (LD_value, LD_av_value, r2_tot, r2_neu, r2_pos, r2_neg)
end

function get_UFE_2nd(L, x1, x2, threshold)
    UFE = zeros(L, L)
    psd_count = 1e-4
    for i in 1:L
        for j in (i+1):L
            x11 = x2[i, j] + psd_count
            x01 = x1[j] - x2[i, j] + psd_count
            x10 = x1[i] - x2[i, j] + psd_count
            x00 = 1 + x2[i, j] - x1[i] - x1[j] + psd_count
            if x00 > threshold && x01 > threshold && x10 > threshold && x11 > threshold
                UFE[i, j] = 1 - log(x11 / x00) / log((x01 * x10) / (x00 * x00))
                UFE[j, i] = UFE[i, j]
            end
        end
    end
    return UFE
end

function get_expected_changes_in_moment(LLhalf, L, S_GT, J_GT, x1, x2, x_1w2_set, mu, d_mu, rec, d_rec)
    S_J_true_flat = zeros(L + LLhalf)
    for i in 1:L
        S_J_true_flat[i] = S_GT[i]
        for j in (i+1):L
            ξ = G(i, j, L) + L
            S_J_true_flat[ξ] = J_GT[i, j]
        end
    end
    Cs_projected_temp = zeros(LLhalf + L)
    ΞS_projected = (x_1w2_set * S_J_true_flat)
    for m in 1:size(x_1w2_set[:, 1], 1)
        Cs_projected_temp += x_1w2_set[m, :] * ΞS_projected[m]
    end
    expected_changes_in_1st_moment = Cs_projected_temp[1:L] - mu * d_mu[1:L] - rec * d_rec[1:L]
    expected_changes_in_2nd_moment = Cs_projected_temp[(L+1):end] - mu * d_mu[(L+1):end] - rec * d_rec[(L+1):end]
    denominator_dot_x_ij_temp = 0; numerator_dot_x_ij_temp = 0
    denominator_dot_x_ij_scaled_temp = 0; numerator_dot_x_ij_scaled_temp = 0
    for i in 1:L
        for j in (i+1):L
            ξ = G(i, j, L)
            denominator_dot_x_ij_temp += (expected_changes_in_1st_moment[i] * expected_changes_in_1st_moment[j])^2
            numerator_dot_x_ij_temp   += (expected_changes_in_2nd_moment[ξ] - (expected_changes_in_1st_moment[i] * x1[j] + expected_changes_in_1st_moment[j] * x1[i]))^2
            denominator_dot_x_ij_scaled_temp += x2[i, j] * (expected_changes_in_1st_moment[i] * expected_changes_in_1st_moment[j])^2
            numerator_dot_x_ij_scaled_temp   += x2[i, j] * (expected_changes_in_2nd_moment[ξ] - (expected_changes_in_1st_moment[i] * x1[j] + expected_changes_in_1st_moment[j] * x1[i]))^2
        end
    end
    return (sqrt(denominator_dot_x_ij_temp), sqrt(numerator_dot_x_ij_temp),
            sqrt(denominator_dot_x_ij_scaled_temp), sqrt(numerator_dot_x_ij_scaled_temp))
end

function get_fitness_at_t(x1_temp, x2_temp, S_GT, J_GT)
    fitness = 0
    for i in 1:L
        fitness += S_GT[i] * x1_temp[i]
        for j in (i+1):L
            fitness += J_GT[i, j] * x2_temp[i, j]
        end
    end
    return fitness + 1
end

function get_fitness_at_t_additive(x1_temp, S_GT)
    fitness = 0
    for i in 1:L
        fitness += S_GT[i] * x1_temp[i]
    end
    return fitness + 1
end

function get_x1_x2_wave(L, data, time_unique)
    n_time_step = length(time_unique)
    x1_t = zeros(L, n_time_step)
    x2_t = zeros(L, L, n_time_step)
    for i_t in 1:n_time_step
        t = time_unique[i_t]
        (n_t, sample_t) = get_sample_at_t(data, t)
        n_scale = sum(n_t)
        for k in 1:length(n_t)
            x1_t[:, i_t]    += (n_t[k] / n_scale) * sample_t[:, k]
            x2_t[:, :, i_t] += (n_t[k] / n_scale) * (sample_t[:, k] * sample_t[:, k]')
        end
    end
    return (x1_t, x2_t)
end

function get_fitness_evolution(L, data, time_unique, S_GT, J_GT)
    n_time_step = length(time_unique)
    Δt = time_unique[2]
    var_fitness_set          = zeros(n_time_step)
    var_fitness_set_additive = zeros(n_time_step)
    change_in_fitness_set    = zeros(n_time_step)
    Fav_1st_2nd_set          = zeros(2, n_time_step)
    Fav_1st_2nd_set_additive = zeros(2, n_time_step)
    f_old = 1
    for i_t in 1:n_time_step
        t = time_unique[i_t]
        (n_t, sample_t) = get_sample_at_t(data, t)
        n_scale = sum(n_t)
        f_temp_2nd, f_temp_1st = 0, 0
        f_temp_2nd_additive, f_temp_1st_additive = 0, 0
        for k in 1:length(n_t)
            f_temp = get_fitness_at_t(sample_t[:, k], sample_t[:, k] * sample_t[:, k]', S_GT, J_GT)
            f_temp_2nd += (n_t[k] / n_scale) * f_temp * f_temp
            f_temp_1st += (n_t[k] / n_scale) * f_temp
            f_temp_additive = get_fitness_at_t_additive(sample_t[:, k], S_GT)
            f_temp_2nd_additive += (n_t[k] / n_scale) * f_temp_additive * f_temp_additive
            f_temp_1st_additive += (n_t[k] / n_scale) * f_temp_additive
        end
        Fav_1st_2nd_set[1, i_t] = f_temp_1st
        Fav_1st_2nd_set[2, i_t] = f_temp_2nd
        Fav_1st_2nd_set_additive[1, i_t] = f_temp_1st_additive
        Fav_1st_2nd_set_additive[2, i_t] = f_temp_2nd_additive
        var_fitness_set[i_t]          = f_temp_2nd - f_temp_1st * f_temp_1st
        var_fitness_set_additive[i_t] = f_temp_2nd_additive - f_temp_1st_additive * f_temp_1st_additive
        change_in_fitness_set[i_t]    = (f_temp_1st - f_old) / Δt
        f_old = f_temp_1st
    end
    return (var_fitness_set, var_fitness_set_additive, change_in_fitness_set,
            Fav_1st_2nd_set, Fav_1st_2nd_set_additive)
end

function get_Drift_of_Fitness(time_unique, i_t_offset, Fav_1st_2nd_set, v_actual)
    n_time_step = length(time_unique)
    Drift_actual = zeros(n_time_step)
    t_offset = time_unique[i_t_offset]
    F_offset = Fav_1st_2nd_set[1, i_t_offset]
    for i_t in 1:n_time_step
        t    = time_unique[i_t]
        F1   = Fav_1st_2nd_set[1, i_t]
        F2   = Fav_1st_2nd_set[2, i_t]
        t_eff  = t - t_offset
        vt     = v_actual * t_eff + F_offset
        Drift_actual[i_t] = F2 - 2 * F1 * vt + vt^2
    end
    return Drift_actual
end

# ── Command-line arguments ────────────────────────────────────────────────────
# Usage: julia run_inference_expanded_grid_Jun2026.jl \
#            i_N i_mu i_sig i_rec data_base_dir out_base_dir i_ens_start i_ens_end time_upto
#
# Changes from May2026 version:
#   1. Incremental C_mpl accumulation — avoids storing the full Ξ matrix
#      (O(N×T×rank_x) memory/compute) by updating C_mpl step-by-step.
#   2. Ground-truth files written once per job, not at every checkpoint.
#   3. Per-checkpoint L×L matrix outputs written as gzip-compressed .txt.gz.
i_N            = parse(Int, ARGS[1])
i_mu           = parse(Int, ARGS[2])
i_σ            = parse(Int, ARGS[3])
i_rec          = parse(Int, ARGS[4])
data_base_dir  = ARGS[5]
out_base_dir   = ARGS[6]
i_ensemble_start = parse(Int, ARGS[7])
i_ensemble_end   = parse(Int, ARGS[8])
time_upto        = parse(Int, ARGS[9])

# ── Parameter grids (must match WF job_code_expanded_grid_Jun2026.jl) ─────────
list_N   = [100, 500, 1000, 5000, 10000]
list_mu  = [1e-4, 5e-4, 1e-3, 5e-3]
list_σ   = collect(range(0.001, 0.1, length=16))
list_rec = vcat([0.0, 0.0001], collect(range(0.0, 0.01, length=16))[2:end])

N        = list_N[i_N]
mu_sim   = list_mu[i_mu]
σ        = list_σ[i_σ]
rec_sim  = list_rec[i_rec]

σ_str    = @sprintf("%.5f", σ)
rec_str  = @sprintf("%.5f", rec_sim)
mu_str   = @sprintf("%.5f", mu_sim)
rho      = 0.1

@printf("N=%d  mu=%.5f  σ=%.5f  rec=%.5f  reps=%d-%d\n",
        N, mu_sim, σ, rec_sim, i_ensemble_start, i_ensemble_end)

# ── Fixed model constants ─────────────────────────────────────────────────────
const L        = 50
const LLhalf   = Int(L * (L - 1) / 2)
const rank_x   = L + LLhalf

# ── Inference hyperparameters ─────────────────────────────────────────────────
const γ             = 1.0
const γ_KNS         = 1.0
const γ_KNSGC       = 1e-2
const i_t_offset    = 15
const freq_threshold = 0.001

time_estimate_epistasis = [30, 60, 90, 150, 210, 300, 450, 600, 810, 990]

# ── Path / key construction ───────────────────────────────────────────────────
N_mu_dir  = @sprintf("N-%d_mu-%.5f", N, mu_sim)
dir_key   = @sprintf("epistasis-strength-%s_recombination_strengh-%s_density-%.2f",
                     σ_str, rec_str, rho)

coupling_file_key = @sprintf("epistasis-strength-%s_recombination_strengh-%s", σ_str, rec_str)
log_file_key      = @sprintf("mu-%s_epistasis-strength-%s_recombination_strengh-%s_density-%.2f",
                              mu_str, σ_str, rec_str, rho)

data_dir     = joinpath(data_base_dir, N_mu_dir, dir_key) * "/"
coupling_dir = joinpath(data_base_dir, "coupling_localfield") * "/"
out_dir_base = joinpath(out_base_dir,  N_mu_dir, dir_key)

# ── Read ground-truth parameters (log / coupling / local-field) ───────────────
function get_S_J_GT_expanded(data_dir, coupling_dir, coupling_key, log_key)
    log_fname = data_dir * "log_" * log_key * "_id-1.txt"
    @show log_fname
    hyper_para = readdlm(log_fname)
    mu_v, rec_v, L_v, N_v, T_v = 0.0, 0.0, 0, 0, 0
    for n in 1:size(hyper_para, 1)
        x = split(string(hyper_para[n, 1]), ":")
        if x[1] == "mu"  mu_v  = parse(Float64, x[2]) end
        if x[1] == "rec" rec_v = parse(Float64, x[2]) end
        if x[1] == "L"   L_v   = parse(Int,     x[2]) end
        if x[1] == "N"   N_v   = parse(Int,     x[2]) end
        if x[1] == "T"   T_v   = parse(Int,     x[2]) end
    end
    @printf("mu:%.1e  rec:%.1e  L:%d  N:%d  T:%d\n", mu_v, rec_v, L_v, N_v, T_v)
    S_GT = vec(readdlm(coupling_dir * "local-field_" * coupling_key * ".txt"))
    J_GT = readdlm(coupling_dir * "coupling_"    * coupling_key * ".txt")
    return (mu_v, rec_v, L_v, N_v, T_v, S_GT, J_GT)
end

(mu, rec, L_log, N_log, T_max, S_GT, J_GT) = get_S_J_GT_expanded(
    data_dir, coupling_dir, coupling_file_key, log_file_key)

# ── Genomic distance matrix ───────────────────────────────────────────────────
pilied_col = repeat(collect(1:L), 1, L)
pilied_row = repeat(collect(1:L)', L, 1)
dist_i_j   = abs.(pilied_col .- pilied_row)

# ── Create output subdirectories ──────────────────────────────────────────────
for subdir in ["", "/detectability", "/x2", "/participation_rate",
               "/QLE_order_parameter", "/effective_population_size",
               "/pairwise_heterogeneity", "/fitness_wave",
               "/averaged_LD", "/averaged_Z", "/changes_in_Z", "/changes_in_LD",
               "/Ground_Trueth", "/LD", "/KNS", "/KNSGC", "/UFE", "/MPL", "/SL",
               "/diversity_time_course"]
    mkpath(out_dir_base * subdir)
end

# ── I/O helpers ───────────────────────────────────────────────────────────────
function read_trajectory_gz(fname)
    open(`gzip -dc $fname`) do io
        readdlm(io)
    end
end

# Accumulate plain-text files for a single batch gzip call at job end.
# Spawning gzip per-file (1600 times) dominates runtime on network filesystems;
# one batch call reduces that overhead to a single process spawn.
const files_to_gz = String[]

function writedlm_gz(fname_gz, data)
    fname_txt = fname_gz[1:end-3]   # strip .gz, write plain text
    writedlm(fname_txt, data)
    push!(files_to_gz, fname_txt)
end

# ── Write ground-truth files once per job ─────────────────────────────────────
# J_GT and S_GT are the same for all replicates and checkpoints, so write once.
let gt_j = out_dir_base * "/Ground_Trueth/epis_sigma-"     * σ_str * "_rec-" * rec_str * ".txt"
    gt_s = out_dir_base * "/Ground_Trueth/additive_sigma-" * σ_str * "_rec-" * rec_str * ".txt"
    isfile(gt_j) || writedlm(gt_j, J_GT)
    isfile(gt_s) || writedlm(gt_s, S_GT)
end

# ── Summary files (one per σ,rec; all reps/checkpoints appended) ──────────────
fname_key_summary = @sprintf("_sigma-%s_rec-%s", σ_str, rec_str)
fout_epis = open(out_dir_base * "/inferred_epistasis_summary" * fname_key_summary * ".txt", "w")
fout_sel  = open(out_dir_base * "/inferred_selection_summary"  * fname_key_summary * ".txt", "w")

@show i_ensemble_start, i_ensemble_end

# ══════════════════════════════════════════════════════════════════════════════
# Per-replicate loop
# ══════════════════════════════════════════════════════════════════════════════
for id_ensemble in i_ensemble_start:i_ensemble_end

    file_key_traj = @sprintf(
        "_mu-%s_epistasis-strength-%s_recombination_strengh-%s_density-%.2f",
        mu_str, σ_str, rec_str, rho)
    fname_traj = @sprintf(
        "%sallele-traject-time-10_id-%d_N-%d%s_wt_freq-selection.txt.gz",
        data_dir, id_ensemble, N, file_key_traj)

    if !isfile(fname_traj)
        @warn @sprintf("Trajectory missing for id-%d: %s", id_ensemble, fname_traj)
        continue
    end

    # ── Load trajectory ───────────────────────────────────────────────────────
    data_raw = try
        read_trajectory_gz(fname_traj)
    catch e
        @warn @sprintf("Failed to read id-%d: %s", id_ensemble, string(e))
        continue
    end
    read_upto = count(data_raw[:, 1] .<= time_upto)
    data      = data_raw[1:read_upto, :]
    time_list = Int.(unique(data[:, 1]))
    time_unique = copy(time_list)

    # ── Fitness wave analysis (every recorded step, dt=10) ────────────────────
    (x1_t, x2_t) = get_x1_x2_wave(L, data, time_unique)
    (var_fitness_set, var_fitness_set_additive,
     change_in_fitness_set, Fav_1st_2nd_set,
     Fav_1st_2nd_set_additive) = get_fitness_evolution(L, data, time_unique, S_GT, J_GT)

    v_F      = mean(var_fitness_set[i_t_offset:end])
    v_F_A    = mean(var_fitness_set_additive[i_t_offset:end])
    v_actual = mean(change_in_fitness_set[i_t_offset:end])
    Drift_actual = get_Drift_of_Fitness(time_unique, i_t_offset, Fav_1st_2nd_set, v_actual)
    Diff     = mean(Drift_actual[i_t_offset:end])

    fname_key_id = @sprintf("_sigma-%s_rec-%s_id-%d", σ_str, rec_str, id_ensemble)

    fout_fw = open(out_dir_base * "/fitness_wave/Fitness_evo" * fname_key_id * ".txt", "w")
    for i_t in 1:length(time_unique)
        vf   = i_t < i_t_offset ? 0.0 : var_fitness_set[i_t]
        vfa  = i_t < i_t_offset ? 0.0 : var_fitness_set_additive[i_t]
        drft = i_t < i_t_offset ? 0.0 : Drift_actual[i_t]
        dcf  = i_t < i_t_offset ? 0.0 : change_in_fitness_set[i_t]
        println(fout_fw, @sprintf("%d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e",
            time_unique[i_t],
            Fav_1st_2nd_set[1, i_t], Fav_1st_2nd_set[2, i_t],
            Fav_1st_2nd_set_additive[1, i_t], Fav_1st_2nd_set_additive[2, i_t],
            vf, vfa, drft, dcf))
    end
    close(fout_fw)
    open(out_dir_base * "/fitness_wave/Averaged_Fitness" * fname_key_id * ".txt", "w") do f
        println(f, @sprintf("%.6e %.6e %.6e %.6e", Diff, v_F, v_F_A, v_actual))
    end

    # ── Genetic diversity time course ─────────────────────────────────────────
    open(out_dir_base * "/diversity_time_course/diversity" * fname_key_id * ".txt", "w") do fout_div
        for i_t in 1:length(time_unique)
            x1_now   = x1_t[:, i_t]
            het      = x1_now .* (1 .- x1_now)
            H_mean   = mean(het)
            H_std    = std(het)
            H_max    = maximum(het)
            n_poly   = count(0.1 .<= x1_now .<= 0.9)
            n_high   = count(x1_now .> 0.9)
            n_low    = count(x1_now .< 0.1)
            n_het_80 = count(het .> 0.2)
            n_het_40 = count(het .> 0.1)
            println(fout_div, @sprintf("%d %.5e %.5e %d %d %d %.5e %d %d",
                time_unique[i_t], H_mean, H_std, n_poly, n_high, n_low, H_max,
                n_het_80, n_het_40))
        end
    end

    # ── Subsample to t % 30 == 0 for statistical accumulation ────────────────
    idx_sub       = data[:, 1] .% 30 .== 0
    data_sub      = data[idx_sub, :]
    time_list_sub = time_list[time_list .% 30 .== 0]
    n_time_sub    = length(time_list_sub)

    fout_pr  = open(out_dir_base * "/participation_rate/participation_rate"               * fname_key_id * ".txt", "w")
    fout_ph  = open(out_dir_base * "/pairwise_heterogeneity/pairwise_heterogeneity"       * fname_key_id * ".txt", "w")
    fout_az  = open(out_dir_base * "/averaged_Z/averaged_Z"                               * fname_key_id * ".txt", "w")
    fout_ald = open(out_dir_base * "/averaged_LD/averaged_LD"                             * fname_key_id * ".txt", "w")
    fout_dz  = open(out_dir_base * "/changes_in_Z/changes_in_Z"                          * fname_key_id * ".txt", "w")
    fout_dld = open(out_dir_base * "/changes_in_LD/changes_in_LD"                        * fname_key_id * ".txt", "w")
    fout_ne  = open(out_dir_base * "/effective_population_size/effective_population_size" * fname_key_id * ".txt", "w")
    fout_qle = open(out_dir_base * "/QLE_order_parameter/QLE_order_parameter"             * fname_key_id * ".txt", "w")

    # ── Initial subsampled time step (id_t = 1) ───────────────────────────────
    t = time_list_sub[1]
    (n_t, sample_t) = get_sample_at_t(data_sub, t)
    Ntemp = sum(n_t)
    println(fout_pr, @sprintf("%d %.5f", t, sum((n_t / Ntemp) .^ 2)))

    (x_1w2_set, x_1w2, d_mu_tot, d_rec_tot) = get_x_1w2_set_and_d_mu_random_compression(
        sample_t, n_t, L, rank_x)
    println(fout_ph, @sprintf("%d %.5f", t, sum(x_1w2[(L+1):end] .^ 2) / LLhalf))

    x1 = copy(x_1w2[1:L])
    x2 = flat_J_to_mat_J(x_1w2, L)
    (Z_temp, Z_av_temp)  = get_averaged_Z(L, x1, x2, freq_threshold)
    (LD_temp, LD_av_temp, r2_tot, r2_neu, r2_pos, r2_neg) = get_averaged_LD(L, x1, x2, J_GT, freq_threshold)
    println(fout_az,  @sprintf("%d %.5e", t, abs(Z_av_temp)))
    println(fout_ald, @sprintf("%d %.5e %.5e %.5e %.5e %.5e", t, abs(LD_av_temp), r2_tot, r2_neu, r2_pos, r2_neg))
    println(fout_dz,  @sprintf("%d %.5e", t, 0.0))
    println(fout_dld, @sprintf("%d %.5e", t, 0.0))
    println(fout_ne,  @sprintf("%d %.5f %.5f", t, 1.0, 1.0))
    println(fout_qle, @sprintf("%d %.5e %.5e %.5e %.5e", t, 0.0, 0.0, 0.0, 0.0))
    Z_old = copy(Z_temp); LD_old = copy(LD_temp)

    x_1w2_init = copy(x_1w2)

    # ── Initialise incremental MPL covariance matrix ──────────────────────────
    # Replaces x_1w2_set_alltime: accumulate C_mpl = Ξ*Ξ' step-by-step as
    #   C_mpl += (Dt_k/2) * x_1w2_set' * x_1w2_set   (trapezoidal weight)
    #   C_mpl += (dt_k/6) * dv * dv'                  (correction term)
    # Memory is O(rank_x²) = 13 MB regardless of N or T.
    Dt_k_init   = time_list_sub[2] - time_list_sub[1]
    C_mpl_accum = (Dt_k_init / 2) * (x_1w2_set' * x_1w2_set)
    d_mu_tot   *= Dt_k_init / 2
    d_rec_tot  *= Dt_k_init / 2

    x_1w2_old     = copy(x_1w2)
    x2_accum_temp = flat_J_to_mat_J(x_1w2, L)
    x1_accum_temp = copy(x_1w2[1:L])

    # ── Loop over remaining subsampled time steps ─────────────────────────────
    for id_t in 2:n_time_sub
        t_old = time_list_sub[id_t - 1]
        t     = time_list_sub[id_t]

        (n_t, sample_t) = get_sample_at_t(data_sub, t)
        N_temp = sum(n_t)
        println(fout_pr, @sprintf("%d %.5f", t, sum((n_t / N_temp) .^ 2)))

        (x_1w2_set, x_1w2, d_mu, d_rec) = get_x_1w2_set_and_d_mu_random_compression(
            sample_t, n_t, L, rank_x)
        println(fout_ph, @sprintf("%d %.5f", t, sum(x_1w2[(L+1):end] .^ 2) / LLhalf))

        x1 = copy(x_1w2[1:L])
        x2 = flat_J_to_mat_J(x_1w2, L)
        x2[diagind(x2)] = x1 .* (1 .- x1)
        C  = x2 - x1 * x1'

        (Z_temp, Z_av_temp)  = get_averaged_Z(L, x1, x2, freq_threshold)
        (LD_temp, LD_av_temp, r2_tot, r2_neu, r2_pos, r2_neg) = get_averaged_LD(L, x1, x2, J_GT, freq_threshold)
        Z_nz  = Z_temp  .!= 0;  LD_nz = LD_temp .!= 0
        ΔZ_av  = count(Z_nz)  > 0 ? sqrt(sum(((Z_temp[Z_nz]  - Z_old[Z_nz])  ) .^ 2) / count(Z_nz))  : 0.0
        ΔLD_av = count(LD_nz) > 0 ? sqrt(sum(((LD_temp[LD_nz] - LD_old[LD_nz])) .^ 2) / count(LD_nz)) : 0.0
        println(fout_az,  @sprintf("%d %.5e", t, abs(Z_av_temp)))
        println(fout_ald, @sprintf("%d %.5e %.5e %.5e %.5e %.5e", t, abs(LD_av_temp), r2_tot, r2_neu, r2_pos, r2_neg))
        println(fout_dz,  @sprintf("%d %.5e", t, ΔZ_av))
        println(fout_dld, @sprintf("%d %.5e", t, ΔLD_av))
        Z_old = copy(Z_temp); LD_old = copy(LD_temp)

        x1_old = copy(x_1w2_old[1:L])
        println(fout_ne, @sprintf("%d %.5e %.5e", t, var(x1_old), 2 * var(x1 .- x1_old)))

        # ── Update incremental C_mpl (trapezoidal weight + correction) ────────
        Dt_k = id_t < n_time_sub ? time_list_sub[id_t+1] - time_list_sub[id_t-1] :
                                    time_list_sub[id_t]   - time_list_sub[id_t-1]
        C_mpl_accum += (Dt_k / 2) * (x_1w2_set' * x_1w2_set)
        d_mu_tot  += d_mu  * (Dt_k / 2)
        d_rec_tot += d_rec * (Dt_k / 2)

        dt_k = time_list_sub[id_t] - time_list_sub[id_t - 1]
        dv   = x_1w2 - x_1w2_old
        C_mpl_accum += (dt_k / 6) * (dv * dv')

        (denom_qle, numer_qle, denom_qle_sc, numer_qle_sc) = get_expected_changes_in_moment(
            LLhalf, L, S_GT, J_GT, x1, x2, x_1w2_set, mu, d_mu, rec, d_rec)
        println(fout_qle, @sprintf("%d %.5e %.5e %.5e %.5e", t, denom_qle, numer_qle, denom_qle_sc, numer_qle_sc))

        x_1w2_old = copy(x_1w2)
        x1_accum_temp += x1
        x2_accum_temp += x2

        # ── Inference at checkpoint times ─────────────────────────────────────
        if t ∈ time_estimate_epistasis
            x1_var   = x1 .* (1 .- x1)
            x1_accum = x1_accum_temp / id_t
            x2_accum = x2_accum_temp / id_t
            C_accum  = x2_accum - x1_accum * x1_accum'
            C_accum[diagind(C_accum)] .= x1_accum .* (1 .- x1_accum)

            fname_key_t = @sprintf("_sigma-%s_rec-%s_id-%d_time-%d_density-%.2f",
                                    σ_str, rec_str, id_ensemble, t, rho)

            # ── LD ─────────────────────────────────────────────────────────────
            LD11_temp  = C       ./ sqrt.(abs.(x1_var * x1_var')                                               .+ 1e-10)
            LD11_accum = C_accum ./ sqrt.(abs.(x1_accum .* (1 .- x1_accum)) .* (x1_accum .* (1 .- x1_accum))' .+ 1e-10)
            LD11_temp[LD11_temp   .>  1.0] .=  1.0;  LD11_temp[LD11_temp   .< -1.0] .= -1.0
            LD11_accum[LD11_accum .>  1.0] .=  1.0;  LD11_accum[LD11_accum .< -1.0] .= -1.0

            # Detectability: x_i * x_j (threshold in accuracy script: > 0.05² = 0.0025)
            detectable_mat       = x1 * x1'
            x1_var_accum         = x1_accum .* (1 .- x1_accum)
            detectable_mat_accum = x1_accum * x1_accum'

            # ── UFE ────────────────────────────────────────────────────────────
            UFE_temp  = get_UFE_2nd(L, x1,       x2,       freq_threshold)
            UFE_accum = get_UFE_2nd(L, x1_accum, x2_accum, freq_threshold)

            # ── KNS ────────────────────────────────────────────────────────────
            cij = dist_i_j
            J_nMF_temp  = try -inv(C       + γ_KNS * I) catch; -inv(C       + 10*γ_KNS * I) end
            J_nMF_accum = try -inv(C_accum + γ_KNS * I) catch; -inv(C_accum + 10*γ_KNS * I) end
            KNS_temp    = J_nMF_temp  .* (rec * cij .+ 4*mu)
            KNS_accum   = J_nMF_accum .* (rec * cij .+ 4*mu)

            # ── KNS+GC ─────────────────────────────────────────────────────────
            KNSGC_temp  = (C       ./ (x1_var       * x1_var'       .+ γ_KNSGC)) .* (4*mu .+ rec * cij)
            KNSGC_accum = (C_accum ./ (x1_var_accum * x1_var_accum' .+ γ_KNSGC)) .* (4*mu .+ rec * cij)

            # ── MPL and SL (C_mpl_accum is already up-to-date) ─────────────────
            my_num = (x_1w2 - x_1w2_init) - mu * d_mu_tot + rec * d_rec_tot

            s_MPL          = (C_mpl_accum + γ * I) \ my_num
            additive_s_MPL = copy(s_MPL[1:L])
            epis_MPL       = flat_J_to_mat_J(s_MPL, L)

            s_SL           = my_num ./ (diag(C_mpl_accum) .+ γ)
            additive_s_SL  = copy(s_SL[1:L])
            epis_SL        = flat_J_to_mat_J(s_SL, L)

            # ── Write per-checkpoint outputs as gzip-compressed .txt.gz ────────
            writedlm_gz(out_dir_base * "/detectability/detectability-temp"  * fname_key_t * ".txt.gz", detectable_mat)
            writedlm_gz(out_dir_base * "/detectability/detectability-accum" * fname_key_t * ".txt.gz", detectable_mat_accum)
            writedlm_gz(out_dir_base * "/x2/x2-temp"                        * fname_key_t * ".txt.gz", x2)
            writedlm_gz(out_dir_base * "/x2/x2-accum"                       * fname_key_t * ".txt.gz", x2_accum)
            writedlm_gz(out_dir_base * "/LD/LD11-temp"                       * fname_key_t * ".txt.gz", LD11_temp)
            writedlm_gz(out_dir_base * "/LD/LD11-accum"                      * fname_key_t * ".txt.gz", LD11_accum)
            writedlm_gz(out_dir_base * "/UFE/UFE-temp"                       * fname_key_t * ".txt.gz", UFE_temp)
            writedlm_gz(out_dir_base * "/UFE/UFE-accum"                      * fname_key_t * ".txt.gz", UFE_accum)
            writedlm_gz(out_dir_base * "/KNS/KNS-temp"                       * fname_key_t * ".txt.gz", KNS_temp)
            writedlm_gz(out_dir_base * "/KNS/KNS-accum"                      * fname_key_t * ".txt.gz", KNS_accum)
            writedlm_gz(out_dir_base * "/KNSGC/KNSGC-temp"                   * fname_key_t * ".txt.gz", KNSGC_temp)
            writedlm_gz(out_dir_base * "/KNSGC/KNSGC-accum"                  * fname_key_t * ".txt.gz", KNSGC_accum)
            writedlm_gz(out_dir_base * "/MPL/epis_MPL"                        * fname_key_t * ".txt.gz", epis_MPL)
            writedlm_gz(out_dir_base * "/SL/epis_SL"                          * fname_key_t * ".txt.gz", epis_SL)
            writedlm_gz(out_dir_base * "/MPL/additive_MPL"                    * fname_key_t * ".txt.gz", additive_s_MPL)
            writedlm_gz(out_dir_base * "/SL/additive_SL"                      * fname_key_t * ".txt.gz", additive_s_SL)

            for i in 1:L
                for j in (i+1):L
                    println(fout_epis, @sprintf(
                        "%d %d %d %d %d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e",
                        id_ensemble, t, i, j, abs(j - i),
                        J_GT[i,j], S_GT[i],
                        LD11_accum[i,j],  LD11_temp[i,j],
                        UFE_accum[i,j],   UFE_temp[i,j],
                        KNS_accum[i,j],   KNS_temp[i,j],
                        KNSGC_accum[i,j], KNSGC_temp[i,j],
                        epis_MPL[i,j], epis_SL[i,j],
                        additive_s_MPL[i], additive_s_SL[i],
                        detectable_mat[i,j], detectable_mat_accum[i,j]))
                end

                println(fout_sel, @sprintf("%d %d %d %.5e %.5e %.5e %.5e",
                    id_ensemble, t, i, S_GT[i], additive_s_MPL[i], additive_s_SL[i], x1[i]))
            end
        end # checkpoint
    end # time loop

    close(fout_pr); close(fout_ph)
    close(fout_az); close(fout_ald)
    close(fout_dz); close(fout_dld)
    close(fout_ne); close(fout_qle)

end # replicate loop

close(fout_epis)
close(fout_sel)

# Compress all matrix output files in one batch gzip call instead of 1600 separate spawns.
isempty(files_to_gz) || run(Cmd(vcat(["gzip", "-f"], files_to_gz)))
