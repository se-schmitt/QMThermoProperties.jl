# translational
function ln_partition_function_trans(ϱ, T, Mw)
    _2πkB_NAh² = 2π*k_B/N_A/h²
    return log((_2πkB_NAh² * Mw * T)^(3/2) / (ϱ*N_A))
end

# rotational
function ln_partition_function_rot(T, σ, I, ::Val{false})
    return log(sqrt(π) / σ * (8*π^2*k_B/h² * T)^(3/2) * sqrt(prod(I * (dalton * 1e-20))))
end

function ln_partition_function_rot(T, σ, I, ::Val{true})
    return log(8*π^2*k_B/h² * T * I[1] * (dalton * 1e-20)  / σ)
end

# vibrational
function ln_partition_function_vib(T, ν)
    _h_kB = h / k_B
    exp_hν_2kBT = -_h_kB .* ν ./ T ./ 2
    ln_1_m_exp_hνᵢ_kBT = log.(1 .- exp.(-_h_kB .* ν ./ T))
    return sum(exp_hν_2kBT) - sum(ln_1_m_exp_hνᵢ_kBT)
end

