@testitem "COSMOSAC" begin
    using Clapeyron

    smiles = ["Cc1ccccc1", "CCO"]
    names = ["toluene", "ethanol"]
    mw = [92.14099, 46.069]

    model02 = COSMOSAC02_CHAOS(names; userlocations=(;SMILES=smiles, Mw=mw))
    model10 = COSMOSAC10_CHAOS(names; userlocations=(;SMILES=smiles, Mw=mw))

    T = 300.
    p = 1e5
    x = [.3,.7]
    γ02 = activity_coefficient(model02, p, T, x)
    γ10 = activity_coefficient(model10, p, T, x)

    @test γ02[1] ≈ 1.811009305776511 rtol=1e-8
    @test γ02[2] ≈ 1.078574899399766 rtol=1e-8
    @test γ10[1] ≈ 1.875489492607190 rtol=1e-8
    @test γ10[2] ≈ 1.084077874193272 rtol=1e-8

    @test COSMOSAC10_CHAOS("ethanol"; userlocations=(;SMILES=["CCO"])) isa COSMOSAC10
    @test_throws ErrorException COSMOSAC02_CHAOS("xenon"; userlocations=(;SMILES=["[Xe]"]))
end
