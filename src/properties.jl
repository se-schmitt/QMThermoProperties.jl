"""
    zpe(model)

Calculates the zero-point energy (ZPE) of the given model:

```math
zpe = a_{\rm ideal}(T \rightarrow 0)
```
"""
function zpe(model)
    ϵ = sqrt(eps(Float64))
    return helmholtz_energy(model, 1., ϵ)
end