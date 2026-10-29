# translational (thermal de Broglie wavelength cubed)
function thermal_wavelength³(T, Mw)
    return (h² / (2π * k_B * Mw/N_A * T))^(3/2)
end

# rotational
function ln_partition_function_rot(T, σ, I)
    θ = 8*π^2*k_B/h² * T * (dalton * 1e-20)
    length(I) == 1 && return log(θ * I[1] / σ)
    return log(sqrt(π * θ^3 * prod(I)) / σ)
end

# vibrational
function ln_partition_function_vib(T, ν)
    θ = h / (k_B * T)
    return -sum(νᵢ -> θ*νᵢ/2 + log1p(-exp(-θ*νᵢ)), ν)
end
