using QMIdealModels, Clapeyron
using Test

@testset "QMIdealModels.jl" begin
    v = 0.025
    T = 300.
    smiles_list = ["CCC","CC1=CC=CC=C1","C(C(F)(F)F)F","COCCl","O=[N+]([O-])C1CCCCC1"]
    names_list = ["water","butane","tert-butylbenzene","dichloroethane","p-dibromobenzene"]

    a_ref =  [203334.178, 250253.866, 45571.9682, 114686.038, 357968.677, 16712.5898, 269023.735, 454935.052, 81002.2349, 125314.802]
    p_ref =   99773.5514
    s_ref =  [231.106926, 297.790687, 279.541771, 264.612555, 344.884919, 149.871513, 263.683428, 363.469799, 263.591498, 324.808877]
    cᵥ_ref = [64.4620409, 101.902287, 77.9280581, 66.9442229, 135.480290, 25.2148861, 85.6266104, 170.741425, 64.8916907, 110.369237]

    models = []
    @testset "Constructors" begin
        for smiles in smiles_list
            model_i = QMIdeal(MolID(smiles))
            @test model_i isa QMIdeal
            push!(models, model_i)
        end
        for name in names_list
            model_i = QMIdeal(name)
            @test model_i isa QMIdeal
            push!(models, model_i)
        end
    end

    @testset "Values" begin
        for (i,model_i) in enumerate(models)
            @test Clapeyron.VT_helmholtz_energy(model_i, v, T) ≈ a_ref[i] rtol=1e-8
            @test Clapeyron.VT_pressure(model_i, v, T) ≈ p_ref rtol=1e-8
            @test Clapeyron.VT_entropy(model_i, v, T) ≈ s_ref[i] rtol=1e-8
            @test Clapeyron.VT_isochoric_heat_capacity(model_i, v, T) ≈ cᵥ_ref[i] rtol=1e-8
        end
    end

    @testset "Gaussian" begin
        p = 101325.             # Pa
        T = 298.15              # K

        # Gaussian reference values
        # 1
        smiles1 = "BrC1C=CC(Br)C1"
        model1 = QMIdeal(MolID(smiles1); vib_correction=1.)
        zpe1 = 259155.8         # J mol⁻¹
        heatcap1 = 24.299       # cal (mol K)⁻¹
        entropy1 = 88.745       # cal (mol K)⁻¹
        e_thermal1 = 66.277     # kcal mol⁻¹
        # 2
        smiles2 = "c1ccc2c(c1)ccc1cc3ccc4ccccc4c3cc12"
        model2 = QMIdeal(MolID(smiles2); vib_correction=1.)
        zpe2 = 763171.4         # J mol⁻¹
        heatcap2 = 63.122       # cal (mol K)⁻¹
        entropy2 = 119.721      # cal (mol K)⁻¹
        e_thermal2 = 191.614    # kcal mol⁻¹
        
        @test zpe1 ≈ helmholtz_energy(model1, p, 1e-6) rtol=1e-4
        @test heatcap1 ≈ isochoric_heat_capacity(model1, p, T) / 4.184 rtol=1e-4
        @test entropy1 ≈ entropy(model1, p, T) / 4.184 rtol=1e-4
        @test e_thermal1 ≈ internal_energy(model1, p, T) / 4184. rtol=1e-4

        @test zpe2 ≈ helmholtz_energy(model2, p, 1e-6) rtol=1e-4
        @test heatcap2 ≈ isochoric_heat_capacity(model2, p, T) / 4.184 rtol=1e-4
        @test entropy2 ≈ entropy(model2, p, T) / 4.184 rtol=1e-4
        @test e_thermal2 ≈ internal_energy(model2, p, T) / 4184. rtol=1e-4
    end
end
