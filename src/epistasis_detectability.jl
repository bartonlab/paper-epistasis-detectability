using Statistics, Random

function get_TP_FP(bool_vec, idx_sort)
    TP_pos = cumsum(bool_vec[idx_sort])
    FP_pos = cumsum(.!bool_vec[idx_sort])
    return (TP_pos / TP_pos[end], FP_pos / FP_pos[end])
end

function get_ROC_AUC(TP, FP)
    myauc = 0.0
    total_positive = maximum(TP)
    total_negative = maximum(FP)
    for n in 2:length(TP)
        y  = 0.5 * (TP[n] + TP[n-1]) / total_positive
        dx = (FP[n] - FP[n-1]) / total_negative
        myauc += y * dx
    end
    return myauc
end

function my_cor(x, y)
    return_cor = 0.0
    if length(x) == length(y)
        if length(unique(x)) > 1 && length(unique(y)) > 1
            return_cor = cor(x, y)
        end
    end
    return return_cor
end

function mean_errorbars(data; n_bootstrap=100, ci=0.90)
    valid_data = filter(x -> !isnan(x) && !ismissing(x), data)
    n_invalid = length(data) - length(valid_data)
    if n_invalid > 0
        println("Warning: $n_invalid NaN or missing values found among $(length(data)) values.")
    end
    means = Float64[]
    n = length(valid_data)
    for _ in 1:n_bootstrap
        resample = rand(Random.GLOBAL_RNG, valid_data, n)
        push!(means, mean(resample))
    end
    lower_percentile = 100 * (1 - ci) / 2
    upper_percentile = 100 * (1 + ci) / 2
    ci_low  = quantile(means, lower_percentile / 100)
    ci_high = quantile(means, upper_percentile / 100)
    mean_val = mean(data)
    return mean_val, mean_val - ci_low, ci_high - mean_val
end
