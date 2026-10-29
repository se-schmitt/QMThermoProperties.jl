module RDKitMinimalLibExt

using QMThermoProperties
using RDKitMinimalLib: RDKitMinimalLib as RDK

function _mol(smiles)
    mol = RDK.get_mol(smiles)
    isnothing(mol) && error("RDKit could not parse SMILES '$(smiles)'")
    return mol
end

_canonicalize_rdk(smiles) = RDK.get_smiles(_mol(smiles))

_get_mw_rdk(smiles) = RDK.get_descriptors(_mol(smiles))["amw"]

function __init__()
    QMThermoProperties._canonicalize[] = _canonicalize_rdk
    QMThermoProperties._get_mw[] = _get_mw_rdk
    return nothing
end

end
