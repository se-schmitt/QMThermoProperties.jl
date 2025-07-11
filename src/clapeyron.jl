abstract type QMIdealModel <: IdealModel end

struct QMIdealParam{T1,T2,T3,T4} <: EoSParam
    ν::SingleParam{Vector{T1}}
    σ::SingleParam{T2}
    I::SingleParam{Vector{T3}}
    is_linear::SingleParam{Bool}
    reference_state::ReferenceState
    Mw::SingleParam{T4}
end

struct QMIdeal{T1,T2,T3,T4} <: QMIdealModel 
    components::Vector{String}
    params::QMIdealParam{T1,T2,T3,T4}
    references::Vector{String}
end

"""
    QMIdeal <: IdealModel

    QMIdeal(components;
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
None

## Description
Ideal model based on data from quantum mechanical (QM) simulations from *Gond et al. (2025)*.

## Model Construction Examples
```
idealmodel = QMIdeal("water")
idealmodel = QMIdeal(["water","carbon dioxide"])
```

## References
1. Gond, D., ...
"""
QMIdeal

export QMIdeal 

CL.default_locations(::Type{QMIdeal}) = ["properties/identifiers.csv", "properties/molarmass.csv"]

function QMIdeal(components; userlocations=String[], reference_state=nothing, verbose=false)

    # Get parameters (Mw and smiles)
    params = getparams(components, CL.default_locations(QMIdeal); userlocations, 
        ignore_headers=["dipprnumber","inchikey","cas"])

    # Load json files


    _components = CL.format_components
    references = String[""]

    model = QMIdeal(_components, params, references)
    CL.set_reference_state!(model, reference_state; verbose)

    return model
end

function CL.a_ideal(model::QMIdeal, V, T, z)
    Σz = sum(z)
    res = zero(V+T+first(z))
    V⁻¹ = 1/V
    ϱ = ∑z * V⁻¹

    ν, σ, I, is_linear, Mw = getfield.(getfield.(Ref(model.params), (:ν,:σ,:I,:is_linear,:Mw)),:values)

    @inbounds for i ∈ CL.@comps
        ln_q_trans = ln_partition_function_trans(ϱ, T, Mw[i]/1000)
        ln_q_rot = ln_partition_function_rot(T, σ[i], I[i], Val(is_linear[i]))
        ln_q_vib = ln_partition_function_vib(T, ν[i])
        res -= z[i]*(ln_q_trans + ln_q_rot + ln_q_vib)
        res += CL.xlogx(z[i],V⁻¹)
    end

    return res
end