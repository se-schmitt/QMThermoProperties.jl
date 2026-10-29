abstract type QMIdealModel <: IdealModel end

@kwdef struct QMIdealParam{TF,TI} <: EoSParam
    ν::SingleParam{Vector{TF}}
    σ::SingleParam{TI}
    I::SingleParam{Vector{TF}}
    ω₀::SingleParam{TI}
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
    vib_correction=0.95461)

## Input parameters

- `components`: Component names
- `userlocations`: path of local database or `NamedTuple` defining the parameters
- `reference_state`: see [`Clapeyron.jl`](https://clapeyronthermo.github.io/Clapeyron.jl/stable/properties/basic/#Reference-States)
- `vib_correction`: correction factor for vibrational frequencies (see Ref. (1) for details)

## Description
Ideal model based on data from quantum mechanical (QM) simulations from *Gond et al. (2025)*.

## Model Construction Examples
```julia
using QMThermoProperties

idealmodel = QMIdeal("water")
idealmodel = QMIdeal(["water","carbon dioxide"])

idealmodel = QMIdeal(smiles"CCCO")
```

## References

1. D. Gond, J. Arweiler, T. Specht, H. Hasse, and F. Jirasek: CHAOS - A Consistent Large-Scale Database for Sigma-Profiles and Other Molecular Descriptors, (2025), DOI: https://doi.org/10.48550/arXiv.2511.19002.

"""
QMIdeal

CL.default_locations(::Type{QMIdeal}) = ["properties/identifiers.csv", "properties/molarmass.csv"]

function QMIdeal(components; 
    userlocations=String[], 
    reference_state=nothing, 
    verbose=false, 
    vib_correction=0.95461
)

    _components = CL.format_components(components)
    
    _params = CL.getparams(components,CL.default_locations(QMIdeal);
        userlocations,ignore_headers=["dipprnumber","inchikey","cas"], ignore_missing_singleparams=["canonicalsmiles", "Mw"])

    entries = chaos_entries(load_chaos(CHAOS_PATH_ID), _params, length(_components))

    _k = h / (8*π^2) * 1e-9 / (dalton * 1e-20)
    ν = SingleParam("ν", _components, [e["Frequencies"] .* (1e2*c₀*vib_correction) for e in entries])
    σ = SingleParam("σ", _components, [e["SymNumber"] for e in entries])
    I = SingleParam("I", _components, [_k ./ e["RotConst"] for e in entries])
    ω₀ = SingleParam("ω₀", _components, [e["Multiplicity"] for e in entries])
    Mw = SingleParam("Mw", _components, [_params["Mw"].ismissingvalues[i] ? get_mw(_params["SMILES"].values[i]) : _params["Mw"].values[i] for i in eachindex(_components)])

    params = QMIdealParam(; ν, I, σ, ω₀, reference_state=ReferenceState(), Mw)
    references = String["10.48550/arXiv.2511.19002"]

    model = QMIdeal(_components, params, references)
    CL.set_reference_state!(model, reference_state; verbose)

    return model
end

function CL.a_ideal(model::QMIdeal, V, T, z)
    (; ν, σ, I, ω₀, Mw) = model.params
    a = zero(V+T+first(z))
    for i ∈ CL.@comps
        ln_q_rot = ln_partition_function_rot(T, σ[i], I[i])
        ln_q_vib = ln_partition_function_vib(T, ν[i])
        ln_q_elec = log(ω₀[i])
        a += CL.xlogx(z[i], N_A*thermal_wavelength³(T, Mw[i]/1000)/V) - z[i]*(ln_q_rot + ln_q_vib + ln_q_elec)
    end
    return a/sum(z) - 1
end