# canonicalization
function _canonicalize_fallback(smiles)
    @warn """
    SMILES are not canonicalized by `QMIdealModels.jl`!
    Either ensure that the used SMILES are canonicalized or load `RDKitMinimalLib.jl` (not working on Windows).
    """
    return smiles
end

const _canonicalize = Ref{Function}(_canonicalize_fallback)

canonicalize(smiles) = _canonicalize[](smiles)

# molecular weight
_get_mw_fallback(smiles) = error("""
    No molecular weight found for SMILES '$(smiles)'. Provide it via `userlocations` or load `RDKitMinimalLib.jl` (not working on Windows).
    """)

const _get_mw = Ref{Function}(_get_mw_fallback)

get_mw(smiles) = _get_mw[](smiles)

# Load databases
function load_chaos_id()
    if !isfile(CHAOS_PATH_ID)
        error("CHAOS database not found at $(CHAOS_PATH_ID)! Download CHAOS database with `QMIdealModels.download_chaos()`")
    end
    
    return JSON.parsefile(CHAOS_PATH_ID)
end

function load_chaos_cosmo02()
    if !isfile(CHAOS_PATH_COSMO02)
        error("CHAOS database not found at $(CHAOS_PATH_COSMO02)! Download CHAOS database with `QMIdealModels.download_chaos()`")
    end
    
    return JSON.parsefile(CHAOS_PATH_COSMO02)
end

function load_chaos_cosmo10()
    if !isfile(CHAOS_PATH_COSMO10)
        error("CHAOS database not found at $(CHAOS_PATH_COSMO10)! Download CHAOS database with `QMIdealModels.download_chaos()`")
    end
    
    return JSON.parsefile(CHAOS_PATH_COSMO10)
end

function create_chaos_json(path; do_save=true)
    JSON.VERSION ≤ v"1" && error("JSON ≥v1 required for this function!")
    filelist = filter(_s -> endswith(_s,".json"), readdir(path))

    chaos_id = Dict{String,Any}()
    chaos_cosmo02 = Dict{String,Any}()
    chaos_cosmo10 = Dict{String,Any}()
    for file in filelist
        add_mol_to_chaos_id!(chaos_id, joinpath(path, file))
        add_mol_to_chaos_cosmo02!(chaos_cosmo02, joinpath(path, file))
        add_mol_to_chaos_cosmo10!(chaos_cosmo10, joinpath(path, file))
    end

    if do_save
        for (_chaos, _path) in zip((chaos_id,chaos_cosmo02,chaos_cosmo10),(CHAOS_PATH_ID,CHAOS_PATH_COSMO02,CHAOS_PATH_COSMO10))
            open(_path, "w") do io 
                JSON.json(io, _chaos; pretty=false)
            end
        end
    end

    return chaos_id, chaos_cosmo02, chaos_cosmo10
end

function add_mol_to_chaos_id!(chaos, filepath)
    name_section = [
        ("MolecularMass", "general"),
        ("not_converged", "general"),
        ("RotConst", "structural"),
        ("SymNumber", "structural"),
        ("Multiplicity", "electronic"),
        ("Frequencies", "vibrational"),
    ]

    return _add_mol_to_chaos!(chaos, filepath, name_section)
end

function add_mol_to_chaos_cosmo02!(chaos, filepath)
    name_section = [
        ("MolecularMass", "general"),
        ("CavArea", "solvation"),
        ("CavVolume", "solvation"),
        ("Sigma_total", "solvation"),
    ]

    return _add_mol_to_chaos!(chaos, filepath, name_section)
end

function add_mol_to_chaos_cosmo10!(chaos, filepath)
    name_section = [
        ("MolecularMass", "general"),
        ("CavArea", "solvation"),
        ("CavVolume", "solvation"),
        ("Fraction_PartSigmas","solvation"),
        ("Sigma_NHB", "solvation"),
        ("Sigma_OH", "solvation"),
        ("Sigma_OT", "solvation"),
    ]

    return _add_mol_to_chaos!(chaos, filepath, name_section)
end

function _add_mol_to_chaos!(chaos, filepath, name_section)
    mol = JSON.parsefile(filepath)
    smiles = canonicalize(mol["general"]["CanonicalSMILES"])
    chaos[smiles] = Dict(name => _save_dict(mol, section, name) for (name,section) in name_section)
    return nothing
end

function _save_dict(mol, section, name)
    if startswith(name, "Sigma")
        _s_arr = sparse(mol[section][name])
        return Dict("idx" => _s_arr.nzind, "val" => _s_arr.nzval)
    else
        return mol[section][name]
    end
end
