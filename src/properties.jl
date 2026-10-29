"""
    zero_point_energy(model::QMIdeal, z=SA[1.])
    zpe(model::QMIdeal, z=SA[1.])

Default units: `[J]`

Calculates the zero-point energy (ZPE), defined as:

```julia
ZPE = A(T → 0) = ∑zᵢN_A∑ₖhνᵢₖ/2
```
"""
function zpe(model, z=CL.SA[1.0])
    ϵ = sqrt(eps(Float64))
    return helmholtz_energy(model, 1., ϵ, z)
end
