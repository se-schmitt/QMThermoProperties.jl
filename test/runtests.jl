using QMIdealModels, Clapeyron
using Test

@testset "QMIdealModels.jl" begin
    smiles_list = ["CCC","CC1=CC=CC=C1","C(C(F)(F)F)F","COCCl","O=[N+]([O-])C1CCCCC1"]
    names_list = ["water","butane","tert-butylbenzene","dichloroethane","p-dibromobenzene"]
    models = []
    @testset "Constructors" begin
        for smiles in smiles_list
            push!(models, QMIdeal(MolID(smiles)))
        end
        for name in names_list
            push!(models, QMIdeal(name))
        end
    end
end
