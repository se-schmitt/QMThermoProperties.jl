# translational
function ln_partition_function_trans(ϱ::TR, T::TT, Mw::TM) where {TR, TT, TM}
    _3_2  = promote_type(TR,TT,TM)(3/2)
    _2πkB_NAh² = promote_type(TR,TT,TM)(2π*kB/NA/h²)
    _NA = promote_type(TN,TT,TM)(NA)
    return log((_2πkB_NAh² * Mw * T)^_3_2 / (ϱ*_NA))
end

# rotational
function ln_partition_function_rot(T::TT, σ::TS, I::TI, ::Val{false}) where {TT,TS,TI}
    _3_2, _sqrtπ = promote_type(TT,TS,TI).((3/2, sqrt(π)))
    _8π²kB_h² = promote_type(TT,TS,TI)(8*π^2*kB/h²/kgm²2uÅ²)
    return log(_sqrtπ / σ * (_8π²kB_h² * T)^_3_2 * sqrt(prod(I)))
end

function ln_partition_function_rot(T::TT, σ::TS, I::TI, ::Val{false}) where {TT,TS,TI}
    _8π²kB_h² = promote_type(TT,TS,TI)(8*π^2*kB/h²/kgm²2uÅ²)
    return log(_8π²kB_h² * T * I[1]  / σ)
end

# vibrational
function ln_partition_function_vib(T::TT, ν::TN) where {TT,TN}
    _h_kB = promote_type(TT,TN)(h / kB)
    exp_hν_2kBT = -_h_kB .* ν ./ T ./ 2
    ln_1_m_exp_hνᵢ_kBT = log.(1 .- exp.(-_h_kB .* ν ./ T))
    return sum(exp_hν_2kBT) - sum(ln_1_m_exp_hνᵢ_kBT)
end

