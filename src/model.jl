abstract type QMIdealModel <: IdealModel end

@kwdef struct QMIdealParam{TF,TI} <: EoSParam
    ν::SingleParam{Vector{TF}}
    σ::SingleParam{TI}
    I::SingleParam{Vector{TF}}
    is_linear::SingleParam{Bool}
    reference_state::ReferenceState
    Mw::SingleParam{TF}
end

struct QMIdeal{TF,TI} <: QMIdealModel 
    components::Vector{String}
    params::QMIdealParam{TF,TI}
    references::Vector{String}
end

"""
    QMIdeal <: IdealModel

    QMIdeal(components;
    userlocations = String[],
    reference_state = nothing,
    verbose = false,
    vib_correction=0.9546)

## Input parameters

- `components`: Components either defined 
    - by name (requires name to be in the `Clapeyron.jl` or a custom database) or
    - by SMILES via `MolID`
- `userlocation`: path of local database or `NamedTuple` defining the parameters
- `reference-state`: see [`Clapeyron.jl`](https://clapeyronthermo.github.io/Clapeyron.jl/stable/properties/basic/#Reference-States)
- `vib_correction`: correction factor for vibrational frequencies (see Ref. (1) for details)

## Description
Ideal model based on data from quantum mechanical (QM) simulations from *Gond et al. (2025)*.

## Model Construction Examples
```julia
using QMIdealModels

idealmodel = QMIdeal("water")
idealmodel = QMIdeal(["water","carbon dioxide"])

idealmodel = QMIdeal(MolID("CCCO"))
```

## References

1. D. Gond, J. Arweiler, T. Specht, H. Hasse, and F. Jirasek: CHAOS - A Consistent Large-Scale Database for Sigma-Profiles and Other Molecular Descriptors, (2025), DOI: https://doi.org/10.48550/arXiv.2511.19002.

"""
QMIdeal

export QMIdeal 

CL.default_locations(::Type{QMIdeal}) = ["properties/identifiers.csv", "properties/molarmass.csv"]

function QMIdeal(components; userlocations=String[], reference_state=nothing, verbose=false, vib_correction=0.9546)

    # Get parameters (Mw and smiles)
    params = getparams(components, CL.default_locations(QMIdeal); userlocations, 
        ignore_headers=["dipprnumber","inchikey","cas"])
    _components = CL.format_components(components)

    return QMIdeal(MolID.(_components, params["canonicalsmiles"]); reference_state, verbose, vib_correction)
end

QMIdeal(component::AbstractMolID; kwargs...) = QMIdeal([component]; kwargs...)
function QMIdeal(components::Vector{<:AbstractMolID}; reference_state=nothing, verbose=false, vib_correction=0.9546)
    # Load json files
    params = get_qm_params(components; vib_correction)

    references = String["Gond et al., to be published."]

    model = QMIdeal(getfield.(components,:name), params, references)
    CL.set_reference_state!(model, reference_state; verbose)

    return model
end

function CL.a_ideal(model::QMIdeal, V, T, z)
    ∑z = sum(z)
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