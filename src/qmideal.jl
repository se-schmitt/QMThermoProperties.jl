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
    vib_correction = 0.95461)

## Input parameters

- `SMILES`: Single Parameter (`String`) - SMILES representation of the components (canonicalized using RDKit)
- `canonicalsmiles`: Single Parameter (`String`) (Optional) - canonical SMILES (RDKit), used instead of `SMILES` if given
- `Mw`: Single Parameter (`Float64`) (Optional) - Molecular Weight `[g·mol⁻¹]`, calculated from `SMILES` using RDKit if not given

## Model parameters

- `ν`: Single Parameter (`Vector{Float64}`) - vibrational frequencies `[s⁻¹]`
- `σ`: Single Parameter (`Int`) - rotational symmetry number
- `I`: Single Parameter (`Vector{Float64}`) - principal moments of inertia `[u·Å²]` (one value for linear molecules)
- `ω₀`: Single Parameter (`Int`) - spin multiplicity of the electronic ground state
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`

## Description

Ideal gas model based on quantum mechanical (QM) data from the CHAOS database [1].
The molecular partition function is evaluated in the rigid-rotor/harmonic-oscillator approximation:

```
A/(nRT) = ∑xᵢ(ln(ρxᵢN_AΛᵢ³) - ln(qᵢ,rot) - ln(qᵢ,vib) - ln(ω₀ᵢ)) - 1
Λᵢ = h/√(2πmᵢk_BT)
qᵢ,rot = √π/σᵢ (8π²k_BT/h²)^(3/2) √(Iᵢ,₁Iᵢ,₂Iᵢ,₃)    (nonlinear)
qᵢ,rot = 8π²k_BTIᵢ/(σᵢh²)                          (linear)
qᵢ,vib = ∏ₖ exp(-hνᵢₖ/2k_BT)/(1-exp(-hνᵢₖ/k_BT))
```

The harmonic frequencies from CHAOS are scaled by `vib_correction` (default from [1]).
The zero-point energy is included in the vibrational contribution (see `zpe`).
The molecules are identified in CHAOS by their canonical SMILES, which requires `RDKitMinimalLib.jl` to be loaded (unless `canonicalsmiles` is given).

## Model Construction Examples
```julia
using QMThermoProperties, Clapeyron, RDKitMinimalLib

# Using the default database
idealmodel = QMIdeal("water") #single input
idealmodel = QMIdeal(["water","carbon dioxide"]) #multiple components
idealmodel = QMIdeal(smiles"CCCO") #by SMILES (Clapeyron database)

# Passing parameters directly
idealmodel = QMIdeal("propanol"; userlocations = (; SMILES = ["CCCO"]))
idealmodel = QMIdeal("propanol"; userlocations = (; canonicalsmiles = ["CCCO"], Mw = [60.096]))

# Unscaled vibrational frequencies
idealmodel = QMIdeal("water"; vib_correction = 1.0)

# As ideal model in an equation of state
model = PR(["water","ethanol"]; idealmodel = QMIdeal)
```

## References

1. Gond, D., Arweiler, J., Specht, T., Hasse, H., & Jirasek, F. (2025). CHAOS – A consistent large-scale database for sigma-profiles and other molecular descriptors. arXiv:2511.19002. [doi:10.48550/arXiv.2511.19002](https://doi.org/10.48550/arXiv.2511.19002)
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