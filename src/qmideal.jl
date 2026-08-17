abstract type QMIdealModel <: IdealModel end

@kwdef struct QMIdealParam{TF,TI} <: EoSParam
    ν::SingleParam{Vector{TF}}
    σ::SingleParam{TI}
    I::SingleParam{Vector{TF}}
    ω₀::SingleParam{TI}
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

- `components`: Component names
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

idealmodel = QMIdeal(smiles"CCCO")
```

## References

1. D. Gond, J. Arweiler, T. Specht, H. Hasse, and F. Jirasek: CHAOS - A Consistent Large-Scale Database for Sigma-Profiles and Other Molecular Descriptors, (2025), DOI: https://doi.org/10.48550/arXiv.2511.19002.

"""
QMIdeal

export QMIdeal 

CL.default_locations(::Type{QMIdeal}) = ["properties/identifiers.csv", "properties/molarmass.csv"]

function QMIdeal(components; 
    userlocations=String[], 
    reference_state=nothing, 
    verbose=false, 
    vib_correction=0.9546
)

    _components = CL.format_components(components)
    
    _params = CL.getparams(components,CL.default_locations(QMIdeal);
        userlocations,ignore_headers=["dipprnumber","inchikey","cas"], ignore_missing_singleparams=["canonicalsmiles", "Mw"])

    smiles = [
        _params["canonicalsmiles"].ismissingvalues[i] ?
        canonicalize.(_params["SMILES"].values[i]) :
        _params["canonicalsmiles"].values[i]
    for i in eachindex(_components)]

    # Load CHAOS db (if not already loaded)
    if isnothing(CHAOS_DB_ID)
        chaos = load_chaos_id()
        global CHAOS_DB_ID = chaos
    else
        chaos = CHAOS_DB_ID
    end
    
    # Check availability
    for s in smiles
        s in keys(chaos) || error("Molecule with SMILES '$(s)' not found in CHAOS database")
        !chaos[s]["not_converged"] || error("Molecule with SMILES '$(s)' marked not converged in the CHAOS database")
    end

    # Create QM parameters
    _chaos_entries = [chaos[s] for s in smiles]
    ν = SingleParam(
        "frequencies in s⁻¹",
        _components,
        [_entry["Frequencies"].*(1e2*c₀) for _entry in _chaos_entries] .* vib_correction
    )
    σ = SingleParam(
        "rotational symmetry number",
        _components,
        [_entry["SymNumber"] for _entry in _chaos_entries]
    )
    _k = h / (8*π^2) * 1e-9 / (dalton * 1e-20)
    I = SingleParam(
        "moments of inertia in u Å²",
        _components,
        [_k ./ _entry["RotConst"] for _entry in _chaos_entries]
    )
    ω₀ = SingleParam(
        "spin multiplicity",
        _components,
        [_entry["Multiplicity"] for _entry in _chaos_entries]
    )
    is_linear = SingleParam(
        "molecule is linear",
        _components,
        [length(i) == 1 for i in I]
    )
    _Mw = [
        _params["Mw"].ismissingvalues[i] ?
        get_mw(_params["SMILES"].values[i]) :
        _params["Mw"].values[i]
    for i in eachindex(_components)]
    Mw = SingleParam("molecular mass in g mol⁻¹", _components, _Mw)

    params = QMIdealParam(; ν, I, σ, ω₀, is_linear, reference_state=ReferenceState(), Mw)

    references = String["Gond et al., to be published."]

    model = QMIdeal(_components, params, references)
    CL.set_reference_state!(model, reference_state; verbose)

    return model
end

function CL.a_ideal(model::QMIdeal, V, T, z)
    ∑z = sum(z)
    x = z ./ ∑z
    a = zero(V+T+first(z))

    ν, σ, I, ω₀, is_linear, Mw = getfield.(getfield.(Ref(model.params), (:ν,:σ,:I,:ω₀,:is_linear,:Mw)),:values)
    g_el = [1]

    @inbounds for i ∈ CL.@comps
        ln_q_trans = ln_partition_function_trans(z[i]/V, T, Mw[i]/1000)
        ln_q_rot = ln_partition_function_rot(T, σ[i], I[i], Val(is_linear[i]))
        ln_q_vib = ln_partition_function_vib(T, ν[i])
        ln_q_elec = log(ω₀[i])

        a -= x[i]*(ln_q_trans + ln_q_rot + ln_q_vib + ln_q_elec + 1)
    end

    return a*∑z
end