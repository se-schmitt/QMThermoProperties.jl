# translational
function ln_partition_function_trans(ϱ, T, Mw)
    _2πkB_NAh² = 2π*CL.k_B/CL.N_A/h²
    return log((_2πkB_NAh² * Mw * T)^(3/2) / (ϱ*CL.N_A))
end

# rotational
function ln_partition_function_rot(T, σ, I, ::Val{false})
    _8π²kB_h² = 8*π^2*CL.k_B/h²/kgm²2uÅ²
    return log(sqrt(π) / σ * (_8π²kB_h² * T)^(3/2) * sqrt(prod(I)))
end

function ln_partition_function_rot(T, σ, I, ::Val{false})
    _8π²kB_h² = 8*π^2*CL.k_B/h²/kgm²2uÅ²
    return log(_8π²kB_h² * T * I[1]  / σ)
end

# vibrational
function ln_partition_function_vib(T, ν)
    _h_kB = CL.h / CL.k_B
    exp_hν_2kBT = -_h_kB .* ν ./ T ./ 2
    ln_1_m_exp_hνᵢ_kBT = log.(1 .- exp.(-_h_kB .* ν ./ T))
    return sum(exp_hν_2kBT) - sum(ln_1_m_exp_hνᵢ_kBT)
end

