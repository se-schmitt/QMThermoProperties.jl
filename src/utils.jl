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

function get_qm_params(comps::Vector{<:MolID}; T=Float64)
    chaos = load_chaos()
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

function download_chaos()
    mkpath(DB_PATH)

    @info "Downloading CHAOS database from $(CHAOS_URL)..."
    
    # Download to temporary file and save only required parts
    #TODO update as soon as CHAOS is online
    # path_tmp = Downloads.download(CHAOS_URL, CHAOS_PATH) 
    path_tmp = joinpath(DB_PATH,"CHAOS_database_all.json")
    raw_chaos = JSON.parsefile(path_tmp)
    stripped_chaos = strip_chaos(raw_chaos)
    open(CHAOS_PATH, "w") do io 
        JSON.print(io, stripped_chaos)
    end
    # rm(path_tmp)
    
    @info "CHAOS database saved to $(CHAOS_PATH)."
end

function strip_chaos(raw)
    stripped = Dict{String,Any}()
    field_part = [
        ("MolecularMass",   "general"),
        ("not_converged",   "general"),
        ("RotConst",        "structural"),
        ("SymNumber",       "structural"),
        ("Frequencies",     "vibrational"),
    ]

    for (smiles,body) in raw
        stripped[smiles] = Dict(field => body[part][field] for (field,part) in field_part)
    end
    return stripped
end