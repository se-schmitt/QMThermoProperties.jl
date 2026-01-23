module QMIdealModels

using JSON, Downloads, RDKitMinimalLib
using Clapeyron
const CL = Clapeyron

const k_B = 1.380649e-23                # J/K
const N_A = 6.02214076e23               # mol⁻¹
const R = k_B * N_A
const h = 6.62607015e-34                # Js
const h² = h^2
const c₀ = 299792458.                   # m/s

const DB_PATH = joinpath(pkgdir(@__MODULE__),"database")
const CHAOS_PATH = joinpath(DB_PATH,"CHAOS_database.json")
const CHAOS_URL = "https://raw.githubusercontent.com/gond/CHAOS/main/qm.json"
CHAOS_DB = nothing

# Conversions
const dalton = 1.66053906892e-27        # kg
const kgm²2uÅ² = 1/1e-20/dalton         # kgm² to uÅ²

include("utils.jl")
include("model.jl")
include("partition_functions.jl")

export QMIdeal, MolID, @mol_str

end
