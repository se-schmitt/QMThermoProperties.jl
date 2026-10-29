"""
    zpe(model, z=[1.0])

Calculates the zero-point energy (ZPE) in J of the ideal gas model `model` for the amounts `z` (in mol):

```math
zpe = A_{\rm ideal}(T \rightarrow 0)
```
"""
function zpe(model, z=CL.SA[1.0])
    ϵ = sqrt(eps(Float64))
    return helmholtz_energy(model, 1., ϵ, z)
end
