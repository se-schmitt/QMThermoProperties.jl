default_locations_cosmo_chaos() = ["properties/identifiers.csv", "properties/molarmass.csv"]

"""
    COSMOSAC02_CHAOS(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters:
- `Pi` :Single Parameter{String} 
- `V`: Single Parameter{Float64}
- `A`: Single Parameter{Float64}

## Description
An activity coefficient model using molecular solvation based on the COSMO-RS method. Sigma profiles are taken from the CHAOS database (see Gond et al. (2025)).
 
## References
1.  A. Klamt: Conductor-like Screening Model for Real Solvents: A New Approach to the Quantitative Calculation of Solvation Phenomena, J. Phys. Chem. 99 (1995) 2224–2235, DOI: https://doi.org/10.1021/j100007a062.
2.  S.-T. Lin and S. I. Sandler: A Priori Phase Equilibrium Prediction from a Segment Contribution Solvation Model, Ind. Eng. Chem. Res. 41 (2002) 899–913, DOI: https://doi.org/10.1021/ie001047w.
3.  D. Gond, J. Arweiler, T. Specht, H. Hasse, and F. Jirasek: CHAOS - A Consistent Large-Scale Database for Sigma-Profiles and Other Molecular Descriptors, (2025), DOI: https://doi.org/10.48550/arXiv.2511.19002.

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
        userlocations,ignore_headers=["dipprnumber","inchikey","cas"], ignore_missing_singleparams=["canonicalsmiles", "Mw"])

    smiles = [
        _params["canonicalsmiles"].ismissingvalues[i] ?
        canonicalize.(_params["SMILES"].values[i]) :
        _params["canonicalsmiles"].values[i]
    for i in eachindex(_components)]

    # Load CHAOS db (if not already loaded)
    if isnothing(CHAOS_DB)
        chaos = load_chaos_cosmo02()
        global CHAOS_DB_COSMO02 = chaos
    else
        chaos = CHAOS_DB_COSMO02
    end

    # Parse parameters
    _A = zeros(length(_components))
    _V = zeros(length(_components))
    _Pi = [zeros(51) for i in 1:length(_components)]
    
    for (i,s) in enumerate(smiles)
        _A[i] = chaos[s]["CavArea"] * a0^2 / 1e-20
        _V[i] = chaos[s]["CavVolume"] * a0^3 / 1e-30
        setindex!.(Ref(_Pi[i]), chaos[s]["Sigma_total"]["val"], chaos[s]["Sigma_total"]["idx"])
    end
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

## Input parameters:
- `Pnhb` :Single Parameter{String}
- `POH` :Single Parameter{String}
- `POT` :Single Parameter{String}
- `V`: Single Parameter{Float64}
- `A`: Single Parameter{Float64}

## Description
An activity coefficient model using molecular solvation based on the COSMO-RS method. Sigma profiles are now split by non-hydrogen bonding, hydrogen acceptor and hydrogen donor. Sigma profiles are taken from the CHAOS database (see Gond et al. (2025)).
 
## References
1.  A. Klamt: Conductor-like Screening Model for Real Solvents: A New Approach to the Quantitative Calculation of Solvation Phenomena, J. Phys. Chem. 99 (1995) 2224–2235, DOI: https://doi.org/10.1021/j100007a062.
2.  S.-T. Lin and S. I. Sandler: A Priori Phase Equilibrium Prediction from a Segment Contribution Solvation Model, Ind. Eng. Chem. Res. 41 (2002) 899–913, DOI: https://doi.org/10.1021/ie001047w.
3.  C.-M. Hsieh, S. I. Sandler, and S.-T. Lin: Improvements of COSMO-SAC for Vapor–Liquid and Liquid–Liquid Equilibrium Predictions, Fluid Phase Equilibria 297 (2010) 90–97, DOI: https://doi.org/10.1016/j.fluid.2010.06.011.
4.  D. Gond, J. Arweiler, T. Specht, H. Hasse, and F. Jirasek: CHAOS - A Consistent Large-Scale Database for Sigma-Profiles and Other Molecular Descriptors, (2025), DOI: https://doi.org/10.48550/arXiv.2511.19002.
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
        userlocations,ignore_headers=["dipprnumber","inchikey","cas"], ignore_missing_singleparams=["canonicalsmiles", "Mw"])

    smiles = [
        _params["canonicalsmiles"].ismissingvalues[i] ?
        canonicalize.(_params["SMILES"].values[i]) :
        _params["canonicalsmiles"].values[i]
    for i in eachindex(_components)]

    # Load CHAOS db (if not already loaded)
    if isnothing(CHAOS_DB)
        chaos = load_chaos_cosmo10()
        global CHAOS_DB_COSMO10 = chaos
    else
        chaos = CHAOS_DB_COSMO10
    end

    # Parse parameters
    _A = zeros(length(_components))
    _V = zeros(length(_components))
    _Pnhb = [zeros(51) for i in 1:length(components)]
    _POH = [zeros(51) for i in 1:length(components)]
    _POT = [zeros(51) for i in 1:length(components)]
    
    for (i,s) in enumerate(smiles)
        _A[i] = chaos[s]["CavArea"] * a0^2 / 1e-20
        _V[i] = chaos[s]["CavVolume"] * a0^3 / 1e-30
        _f = chaos[s]["Fraction_PartSigmas"]
        setindex!.(Ref(_Pnhb[i]), chaos[s]["Sigma_NHB"]["val"].*_A[i].*_f[1], chaos[s]["Sigma_NHB"]["idx"])
        setindex!.(Ref(_POH[i]), chaos[s]["Sigma_OH"]["val"].*_A[i].*_f[2], chaos[s]["Sigma_OH"]["idx"])
        setindex!.(Ref(_POT[i]), chaos[s]["Sigma_OT"]["val"].*_A[i].*_f[3], chaos[s]["Sigma_OT"]["idx"])
    end
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