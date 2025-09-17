module QMIdealModels

using JSON, Downloads, RDKitMinimalLib
using Clapeyron
const CL = Clapeyron

const h² = CL.h^2
const c₀ = 299792458.                   # m/s

const DB_PATH = joinpath(pkgdir(@__MODULE__),"database")
const CHAOS_PATH = joinpath(DB_PATH,"CHAOS_database.json")
const CHAOS_URL = "https://raw.githubusercontent.com/gond/CHAOS/main/qm.json"

# Conversions
const dalton = 1.66053906892e-27        # kg
const kgm²2uÅ² = 1/1e-20/dalton         # kgm² to uÅ²

function __init__()
    if !isfile(CHAOS_PATH)
        download_chaos()
    end
    return nothing
end

include("utils.jl")
include("model.jl")
include("partition_functions.jl")

export MolID, QMIdeal, @smiles_str

end
