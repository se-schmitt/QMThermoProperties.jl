canonicalize(smiles) = get_smiles(get_mol(smiles))

function load_chaos()
    if !isfile(CHAOS_PATH_ID)
        error("CHAOS database not found at $(CHAOS_PATH_ID)! Download CHAOS database with `QMIdealModels.download_chaos()`")
    end
    
    return JSON.parsefile(CHAOS_PATH_ID)
end

function create_chaos_json(path; do_save=true)
    JSON.VERSION ≤ v"1" && error("JSON ≥v1 required for this function!")
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
        ("Multiplicity",    "electronic"),
        ("Frequencies",     "vibrational"),
    ]

    mol = JSON.parsefile(filepath)
    smiles = canonize_smiles(mol["general"]["CanonicalSMILES"])
    chaos[smiles] = Dict(name => mol[section][name] for (name,section) in name_section)
    return nothing
end
