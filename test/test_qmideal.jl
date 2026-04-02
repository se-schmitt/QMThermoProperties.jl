@testitem "QMIdeal" begin
    using Clapeyron

    v = 0.025
    T = 300.
    smiles_list = ["CCC","CC1=CC=CC=C1","C(C(F)(F)F)F","COCCl","O=[N+]([O-])C1CCCCC1"]
    names_list = ["water","butane","tert-butylbenzene","dichloroethane","p-dibromobenzene"]

    a_ref =  [191638.524, 238558.212, 33876.3144, 102990.384, 346273.023, 5016.93592, 257328.082, 443239.398, 69306.5811, 113619.148]
    p_ref =   99773.5514
    s_ref =  [270.092439, 336.776200, 318.527284, 303.598068, 383.870432, 188.857027, 302.668941, 402.455312, 302.577011, 363.794390]
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
        smiles =  ("BrC1C=CC(Br)C1", "O=C=O", "O=O",   "c1ccc2c(c1)ccc1cc3ccc4ccccc4c3cc12")
        zpe_ =    (259155.8,         31304.7, 10263.6, 763171.4,                           )
        heatcap = (24.299,           6.774,   5.003,   63.122,                             )
        entr =    (88.745,           50.983,  48.942,  119.721,                            )
        e_therm = (66.277,           9.115 ,  3.936,   191.614,                            )
        
        for i in eachindex(smiles)
            s = smiles[i]
            @testset "SMILES: $(length(s) ≤ 15 ? s : s[1:15]*"...")" begin
                model = QMIdeal(MolID(s); vib_correction=1)
                @test zpe_[i] ≈ zpe(model) rtol=5e-4
                @test heatcap[i] ≈ isochoric_heat_capacity(model, p, T) / 4.184 rtol=5e-4
                @test entr[i] ≈ entropy(model, p, T) / 4.184 rtol=5e-4
                @test e_therm[i] ≈ internal_energy(model, p, T) / 4184. rtol=5e-4
            end
        end
    end
end
