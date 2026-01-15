abstract type AbstractMolID end

struct MolID{T1,T2} <: AbstractMolID
    name::T1
    smiles::T2
    MolID(name, _smiles; canonicalize=true) = begin
        if canonicalize
            mol = get_mol(_smiles)
            smiles = get_smiles(mol)
        else
            smiles = _smiles
        end
        return new{typeof(name),typeof(smiles)}(name, smiles)
    end
end
MolID(smiles; canonicalize=true) = MolID(smiles, smiles; canonicalize)
macro mol_str(str, canonicalize=true)
    MolID(str; canonicalize)
end

canonize_smiles(smiles) = get_smiles(get_mol(smiles))

function get_qm_params(comps::Vector{<:MolID}; T=Float64)
    if isnothing(CHAOS_DB)
        chaos = load_chaos()
        global CHAOS_DB = chaos 
    else
        chaos = CHAOS_DB
    end
    names = getfield.(comps, :name)

    smiles = getfield.(comps, :smiles)
    for s in smiles
        s in keys(chaos) || error("Molecule with SMILES '$(s)' not found in CHAOS database")
        !chaos[s]["not_converged"] || error("Molecule with SMILES '$(s)' marked not converged in the CHAOS database")
    end

    Mw = SingleParam(
        "molecular mass in kg mol⁻¹",
        names,
        [T.(chaos[c.smiles]["MolecularMass"]) for c in comps]
    )
    ν = SingleParam(
        "frequencies in s⁻¹",
        names,
        [T.(chaos[c.smiles]["Frequencies"].*(1e2*c₀)) for c in comps]
    )
    _k = T(CL.h / (8*π^2) * 1e11 / dalton)
    I = SingleParam(
        "moments of inertia in ??",     #TODO
        names,
        [_k ./ T.(chaos[c.smiles]["RotConst"]) for c in comps]
    )
    σ = SingleParam(
        "rotational symmetry number",
        names,
        [chaos[c.smiles]["SymNumber"] for c in comps]
    )
    is_linear = SingleParam(
        "molecule is linear",
        names,
        [length(i) == 1 for i in I]
    )

    return QMIdealParam(; Mw, ν, I, σ, is_linear, reference_state=ReferenceState())
end

function load_chaos()
    if !isfile(CHAOS_PATH)
        error("CHAOS database not found at $(CHAOS_PATH)! Download CHAOS database with `QMIdealModels.download_chaos()`")
    end
    
    return JSON.parsefile(CHAOS_PATH)
end

function create_chaos_json(path; do_save=true)
    filelist = filter(_s -> endswith(_s,".json"), readdir(path))

    chaos = Dict{String,Any}()
    for file in filelist
        add_mol_to_chaos!(chaos, joinpath(path, file))
    end

    if do_save
        open(CHAOS_PATH, "w") do io 
            JSON.json(io, chaos; pretty=false)
        end
    end

    return chaos
end

function add_mol_to_chaos!(chaos, filepath)
    name_section = [
        ("MolecularMass",   "general"),
        ("not_converged",   "general"),
        ("RotConst",        "structural"),
        ("SymNumber",       "structural"),
        ("Frequencies",     "vibrational"),
    ]

    mol = JSON.parsefile(filepath)
    smiles = canonize_smiles(mol["general"]["CanonicalSMILES"])
    chaos[smiles] = Dict(name => mol[section][name] for (name,section) in name_section)
    return nothing
end
