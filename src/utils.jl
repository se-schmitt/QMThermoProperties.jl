# canonicalization
function _canonicalize_fallback(smiles)
    @warn """
    SMILES are not canonicalized by `QMThermoProperties.jl`!
    Either ensure that the used SMILES are canonicalized or load `RDKitMinimalLib.jl` (not working on Windows).
    """ maxlog=1
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

# Download and convert database
"""
    download_chaos(; force=false, raw=false)

Make the CHAOS database available to `QMThermoProperties.jl`. The databases are stored in a
scratch space and persist across sessions; call this once before building any model.

By default the already converted databases (~60 MB) are downloaded from
[the release assets]($(CHAOS_DB_URL)) and their checksum is verified.

With `raw=true` the full CHAOS database (~4.3 GB) is instead downloaded from Zenodo
($(CHAOS_URL)) and converted locally, keeping only the properties the models need. This
is how the converted databases above are produced; it needs ~16 GB of free disk space
(the extracted raw data is removed again afterwards) and takes considerably longer.
Use `bundle_chaos` to package the result for upload.

Pass `force=true` to replace an existing database.
"""
function download_chaos(; force=false, raw=false)
    if !force && all(isfile, (CHAOS_PATH_ID, CHAOS_PATH_COSMO02, CHAOS_PATH_COSMO10))
        @info "CHAOS database already available in $(DB_PATH). Use `force=true` to rebuild it."
        return DB_PATH
    end

    raw ? _build_chaos() : _fetch_chaos()

    # invalidate anything loaded from a previous database
    global CHAOS_DB_ID = nothing
    global CHAOS_DB_COSMO02 = nothing
    global CHAOS_DB_COSMO10 = nothing

    @info "CHAOS database ready in $(DB_PATH)."
    return DB_PATH
end

# Download the pre-converted databases
function _fetch_chaos(url=CHAOS_DB_URL, sha=CHAOS_DB_SHA256)
    mktempdir(DB_PATH) do tmp
        archive = joinpath(tmp, "databases.zip")
        @info "Downloading converted CHAOS databases from $(url) ..."
        _dl_last[] = -1.0
        Downloads.download(url, archive; progress=_download_progress)
        _verify_sha256(archive, sha)
        run(`$(p7zip()) x -y -bso0 -bsp0 -o$(DB_PATH) $(archive)`)
    end
    return nothing
end

# Rebuild the converted databases from the raw CHAOS data
function _build_chaos()
    archive = joinpath(DB_PATH, "CHAOS_raw.zip")
    if isfile(archive)
        @info "Reusing already downloaded archive $(archive)"
    else
        @info "Downloading raw CHAOS database from $(CHAOS_URL) ..."
        _dl_last[] = -1.0
        Downloads.download(CHAOS_URL, archive*".part"; progress=_download_progress)
        mv(archive*".part", archive)    # only a complete download is reusable
    end

    # the extracted raw data (~11 GB) is always cleaned up again
    mktempdir(DB_PATH) do raw
        @info "Extracting archive ..."
        run(`$(p7zip()) x -y -bso0 -bsp0 -o$(raw) $(archive)`)

        @info "Converting molecules ..."
        create_chaos_json(raw)
    end
    rm(archive)                         # a failed conversion keeps the download for a retry
    return nothing
end

# run after `download_chaos(raw=true)` (internal use only)
"""
    bundle_chaos(out=joinpath(DB_PATH, "CHAOS_databases.zip"))

Package the converted databases into a single zip file for upload as a release asset,
and report its SHA256 checksum (to be stored in `CHAOS_DB_SHA256`). Counterpart to
`download_chaos(; raw=true)`.
"""
function bundle_chaos(out=joinpath(DB_PATH, "CHAOS_databases.zip"))
    for path in (CHAOS_PATH_ID, CHAOS_PATH_COSMO02, CHAOS_PATH_COSMO10)
        isfile(path) || error("No CHAOS database at $(path)! Build it with `download_chaos(; raw=true)`.")
    end
    rm(out; force=true)
    # stored without paths, so that `_fetch_chaos` can unpack them straight into `DB_PATH`
    run(`$(p7zip()) a -tzip -mx=9 -bso0 -bsp0 $(out) $(CHAOS_PATH_ID) $(CHAOS_PATH_COSMO02) $(CHAOS_PATH_COSMO10)`)
    sha = bytes2hex(open(sha256, out))
    @info "Wrote $(out) ($(round(filesize(out)/1e6; digits=1)) MB)\n  SHA256: $(sha)"
    return out, sha
end

function _verify_sha256(path, expected)
    if isempty(expected)
        @warn "No checksum pinned in `CHAOS_DB_SHA256`; skipping verification of $(basename(path))."
        return nothing
    end
    got = bytes2hex(open(sha256, path))
    got == expected || error("""
        Checksum mismatch for $(basename(path)):
          expected $(expected)
          got      $(got)
        The download is corrupted or the release asset changed. Please retry.
        """)
    return nothing
end

function _download_progress(total, now)
    total == 0 && return nothing
    scale, unit = total < 1e9 ? (1e6, "MB") : (1e9, "GB")
    done = round(now/scale; digits=total < 1e9 ? 1 : 2)
    done == _dl_last[] && return nothing        # only report actual progress
    _dl_last[] = done
    print("\r  $(done) / $(round(total/scale; digits=total < 1e9 ? 1 : 2)) $(unit)")
    now == total && println()
    return nothing
end

const _dl_last = Ref(-1.0)

# Load databases
function _load_chaos(path)
    if !isfile(path)
        error("CHAOS database not found at $(path)! Download CHAOS database with `QMIdealModels.download_chaos()`")
    end

    return JSON.parsefile(path)
end

load_chaos_id() = _load_chaos(CHAOS_PATH_ID)
load_chaos_cosmo02() = _load_chaos(CHAOS_PATH_COSMO02)
load_chaos_cosmo10() = _load_chaos(CHAOS_PATH_COSMO10)

function create_chaos_json(path; do_save=true)
    # The database is keyed by canonicalized SMILES and looked up with `canonicalize` too, so
    # building it without `RDKitMinimalLib.jl` yields keys that no lookup can ever match.
    _canonicalize[] === _canonicalize_fallback && error("""
        Cannot build the CHAOS database without SMILES canonicalization.
        Load `RDKitMinimalLib.jl` (not working on Windows) before calling this function.
        """)

    # molecule files are named `<id>.json`; the archive also ships e.g. `DICT_SMILES_TO_ID.json`
    filelist = filter(_s -> !isnothing(match(r"^\d+\.json$", _s)), readdir(path))

    chaos_id = Dict{String,Any}()
    chaos_cosmo02 = Dict{String,Any}()
    chaos_cosmo10 = Dict{String,Any}()
    for (i,file) in enumerate(filelist)
        mol = JSON.parsefile(joinpath(path, file))     # parse once, feed all three databases
        if !haskey(mol, "general")
            @warn "Skipping '$(file)': not a CHAOS molecule file"
            continue
        end
        add_mol_to_chaos_id!(chaos_id, mol)
        add_mol_to_chaos_cosmo02!(chaos_cosmo02, mol)
        add_mol_to_chaos_cosmo10!(chaos_cosmo10, mol)
        i % 5000 == 0 && @info "  converted $(i)/$(length(filelist)) molecules"
    end

    if do_save
        for (_chaos, _path) in zip((chaos_id,chaos_cosmo02,chaos_cosmo10),(CHAOS_PATH_ID,CHAOS_PATH_COSMO02,CHAOS_PATH_COSMO10))
            open(_path, "w") do io
                JSON.print(io, _chaos)
            end
        end
    end

    return chaos_id, chaos_cosmo02, chaos_cosmo10
end

function add_mol_to_chaos_id!(chaos, mol)
    name_section = [
        ("MolecularMass", "general"),
        ("not_converged", "general"),
        ("RotConst", "structural"),
        ("SymNumber", "structural"),
        ("Multiplicity", "electronic"),
        ("Frequencies", "vibrational"),
    ]

    return _add_mol_to_chaos!(chaos, mol, name_section)
end

function add_mol_to_chaos_cosmo02!(chaos, mol)
    name_section = [
        ("MolecularMass", "general"),
        ("CavArea", "solvation"),
        ("CavVolume", "solvation"),
        ("Sigma_total", "solvation"),
    ]

    return _add_mol_to_chaos!(chaos, mol, name_section)
end

function add_mol_to_chaos_cosmo10!(chaos, mol)
    name_section = [
        ("MolecularMass", "general"),
        ("CavArea", "solvation"),
        ("CavVolume", "solvation"),
        ("Fraction_PartSigmas","solvation"),
        ("Sigma_NHB", "solvation"),
        ("Sigma_OH", "solvation"),
        ("Sigma_OT", "solvation"),
    ]

    return _add_mol_to_chaos!(chaos, mol, name_section)
end

# accepts either a path to a CHAOS molecule file or an already parsed one
_chaos_mol(mol::AbstractString) = JSON.parsefile(mol)
_chaos_mol(mol) = mol

function _add_mol_to_chaos!(chaos, mol, name_section)
    mol = _chaos_mol(mol)
    smiles = canonicalize(mol["general"]["CanonicalSMILES"])
    chaos[smiles] = Dict(name => _save_dict(mol, section, name) for (name,section) in name_section)
    return nothing
end

function _save_dict(mol, section, name)
    if startswith(name, "Sigma")
        _s_arr = sparse(_narrow(mol[section][name]))
        return Dict("idx" => _s_arr.nzind, "val" => _narrow(_s_arr.nzval))
    else
        return _narrow(mol[section][name])
    end
end

# JSON parses arrays as `Vector{Any}`; narrowing them keeps the databases small in memory
_narrow(x::AbstractVector) = identity.(x)
_narrow(x) = x
