module QMThermoProperties

using JSON, Downloads, SHA, Scratch, p7zip_jll
using Clapeyron
const CL = Clapeyron

const k_B = 1.380649e-23                # J/K
const N_A = 6.02214076e23               # mol⁻¹
const h = 6.62607015e-34                # Js
const h² = h^2
const c₀ = 299792458.                   # m/s
const a0 = 5.29177210545e-11            # Bohr radius

# Raw CHAOS database (~4.3 GB) — only needed to rebuild the converted databases
const CHAOS_URL = "https://zenodo.org/records/18701299/files/CHAOS_v2.zip?download=1"

# Converted databases (~60 MB), built from the raw data with `download_chaos(; raw=true)`
const CHAOS_DB_URL = "https://github.com/se-schmitt/QMThermoProperties.jl/releases/download/db-v2/CHAOS_v2_databases.zip"
const CHAOS_DB_SHA256 = "7b57eded570420f416126e62799bf1928416761ed52e26459aad98dd6f8e4c8d"

# Database location (scratch space, filled in `__init__`)
DB_PATH = ""
CHAOS_PATH_ID = ""
CHAOS_PATH_COSMO02 = ""
CHAOS_PATH_COSMO10 = ""

function __init__()
    global DB_PATH = @get_scratch!("CHAOS")
    global CHAOS_PATH_ID = joinpath(DB_PATH,"CHAOS_ID_database.json")
    global CHAOS_PATH_COSMO02 = joinpath(DB_PATH,"CHAOS_COSMO02_database.json")
    global CHAOS_PATH_COSMO10 = joinpath(DB_PATH,"CHAOS_COSMO10_database.json")
    return nothing
end

# Conversions
const dalton = 1.66053906892e-27        # kg

include("utils.jl")
include("qmideal.jl")
include("partition_functions.jl")
include("properties.jl")
include("cosmosac.jl")

export QMIdeal, COSMOSAC02_CHAOS, COSMOSAC10_CHAOS, zpe

end
