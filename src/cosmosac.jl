default_locations_cosmo_chaos() = ["properties/identifiers.csv"]

"""
    COSMOSAC02_CHAOS(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters

- `SMILES`: Single Parameter (`String`) - SMILES representation of the components (canonicalized using RDKit)
- `canonicalsmiles`: Single Parameter (`String`) (Optional) - canonical SMILES (RDKit), used instead of `SMILES` if given

## Input models

- `puremodel`: model to calculate pure pressure-dependent properties

## Description

Constructs a `COSMOSAC02` model (COSMO-SAC 2002) with the parameters (`Pi`, `V`, `A`) taken from the CHAOS database [1].
See [`COSMOSAC02`](https://clapeyronthermo.github.io/Clapeyron.jl/stable/eos/activity/#Clapeyron.COSMOSAC02) in the Clapeyron.jl documentation for details on the model and its references.

The molecules are identified in CHAOS by their canonical SMILES, which requires `RDKitMinimalLib.jl` to be loaded (unless `canonicalsmiles` is given).

## Model Construction Examples
```julia
using QMThermoProperties, Clapeyron, RDKitMinimalLib

model = COSMOSAC02_CHAOS(["toluene","ethanol"])
model = COSMOSAC02_CHAOS(["toluene","ethanol"]; userlocations = (; SMILES = ["CC1=CC=CC=C1","CCO"]))
```

## References

1. Gond, D., Arweiler, J., Specht, T., Hasse, H., & Jirasek, F. (2025). CHAOS – A consistent large-scale database for sigma-profiles and other molecular descriptors. arXiv:2511.19002. [doi:10.48550/arXiv.2511.19002](https://doi.org/10.48550/arXiv.2511.19002)
"""
function COSMOSAC02_CHAOS(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing
)

    _components = CL.format_components(components)
    
    _params = CL.getparams(components, default_locations_cosmo_chaos();
        userlocations,ignore_headers=["dipprnumber","inchikey","cas"], ignore_missing_singleparams=["canonicalsmiles"])

    entries = chaos_entries(load_chaos(CHAOS_PATH_COSMO02), _params, length(_components))

    _A = [e["CavArea"] * a0^2 / 1e-20 for e in entries]
    _V = [e["CavVolume"] * a0^3 / 1e-30 for e in entries]
    _Pi = [sigma_profile(e["Sigma_total"]) for e in entries]
    A = SingleParam("A",_components,_A)
    V = SingleParam("V",_components,_V)
    Pi = SingleParam("Pi",_components,_Pi)

    _puremodel = CL.init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = CL.COSMOSAC02Param(Pi,V,A)
    references = String["10.1021/j100007a062","10.1021/ie001047w","10.48550/arXiv.2511.19002"]
    model = COSMOSAC02(_components, packagedparams, _puremodel, 1e-12, references)
    set_reference_state!(model,reference_state,verbose = verbose)

    return model
end

"""
    COSMOSAC10_CHAOS(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters

- `SMILES`: Single Parameter (`String`) - SMILES representation of the components (canonicalized using RDKit)
- `canonicalsmiles`: Single Parameter (`String`) (Optional) - canonical SMILES (RDKit), used instead of `SMILES` if given

## Input models

- `puremodel`: model to calculate pure pressure-dependent properties

## Description

Constructs a `COSMOSAC10` model (COSMO-SAC 2010) with the parameters (`Pnhb`, `POH`, `POT`, `V`, `A`) taken from the CHAOS database [1].
See [`COSMOSAC10`](https://clapeyronthermo.github.io/Clapeyron.jl/stable/eos/activity/#Clapeyron.COSMOSAC10) in the Clapeyron.jl documentation for details on the model and its references.

The molecules are identified in CHAOS by their canonical SMILES, which requires `RDKitMinimalLib.jl` to be loaded (unless `canonicalsmiles` is given).

## Model Construction Examples
```julia
using QMThermoProperties, Clapeyron, RDKitMinimalLib

model = COSMOSAC10_CHAOS(["toluene","ethanol"])
model = COSMOSAC10_CHAOS(["toluene","ethanol"]; userlocations = (; SMILES = ["CC1=CC=CC=C1","CCO"]))
```

## References

1. Gond, D., Arweiler, J., Specht, T., Hasse, H., & Jirasek, F. (2025). CHAOS – A consistent large-scale database for sigma-profiles and other molecular descriptors. arXiv:2511.19002. [doi:10.48550/arXiv.2511.19002](https://doi.org/10.48550/arXiv.2511.19002)
"""
function COSMOSAC10_CHAOS(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing
)

    _components = CL.format_components(components)
    
    _params = CL.getparams(components, default_locations_cosmo_chaos();
        userlocations,ignore_headers=["dipprnumber","inchikey","cas"], ignore_missing_singleparams=["canonicalsmiles"])

    entries = chaos_entries(load_chaos(CHAOS_PATH_COSMO10), _params, length(_components))

    _A = [e["CavArea"] * a0^2 / 1e-20 for e in entries]
    _V = [e["CavVolume"] * a0^3 / 1e-30 for e in entries]
    _Pnhb = [sigma_profile(e["Sigma_NHB"], A*e["Fraction_PartSigmas"][1]) for (e,A) in zip(entries,_A)]
    _POH = [sigma_profile(e["Sigma_OH"], A*e["Fraction_PartSigmas"][2]) for (e,A) in zip(entries,_A)]
    _POT = [sigma_profile(e["Sigma_OT"], A*e["Fraction_PartSigmas"][3]) for (e,A) in zip(entries,_A)]
    A = SingleParam("A",_components,_A)
    V = SingleParam("V",_components,_V)
    Pnhb = SingleParam("Pnhb", _components, _Pnhb)
    POH = SingleParam("POH", _components, _POH)
    POT = SingleParam("POT", _components, _POT)

    _puremodel = CL.init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = CL.COSMOSAC10Param(Pnhb,POH,POT,V,A)
    references = String["10.1021/j100007a062","10.1021/ie001047w","10.1016/j.fluid.2010.06.011","10.48550/arXiv.2511.19002"]
    model = COSMOSAC10(_components, packagedparams, _puremodel, 1e-12, references)
    set_reference_state!(model,reference_state,verbose = verbose)

    return model
end