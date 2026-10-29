@testitem "QMIdeal" begin
    using Clapeyron, RDKitMinimalLib

    v = 0.025
    T = 300.
    smiles_list = ["CCC","CC1=CC=CC=C1","C(C(F)(F)F)F","COCCl","O=[N+]([O-])C1CCCCC1"]
    names_list = ["water","butane","tert-butylbenzene","dichloroethane","p-dibromobenzene"]

    a_ref =  [191638.383, 238558.507, 33876.8313, 102968.511, 346275.460, 5016.58135, 257328.962, 443242.147, 69270.5231, 113588.971]
    p_ref =   99773.5514
    s_ref =  [270.101755, 336.786009, 318.529264, 303.677052, 383.877003, 188.860098, 302.677268, 402.464052, 302.702114, 363.901603]
    cᵥ_ref = [64.4614137, 101.901123, 77.9274818, 66.9437141, 135.478792, 25.2148707, 85.6257697, 170.739593, 64.8911967, 110.368298]

    models = []
    @testset "Constructors" begin
        for (name,smiles) in zip(names_list,smiles_list)
            model_i = QMIdeal(name, userlocations=(;SMILES=[smiles]))
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

    @testset "Mixtures" begin
        model = QMIdeal(["water","carbon dioxide"])
        z = [0.3,0.7]
        @test Clapeyron.VT_helmholtz_energy(model, 2v, T, 2z) ≈ 2*Clapeyron.VT_helmholtz_energy(model, v, T, z) rtol=1e-12
        @test isfinite(Clapeyron.VT_helmholtz_energy(model, v, T, [0.,1.]))
        @test Clapeyron.VT_entropy(model, v, T, [0.,1.]) ≈ Clapeyron.VT_entropy(QMIdeal("carbon dioxide"), v, T) rtol=1e-12
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
                model = QMIdeal(s; userlocations=(;SMILES=[s]), vib_correction=1)
                @test zpe_[i] ≈ zpe(model) rtol=5e-4
                @test heatcap[i] ≈ isochoric_heat_capacity(model, p, T) / 4.184 rtol=5e-4
                @test entr[i] ≈ entropy(model, p, T) / 4.184 rtol=5e-4
                @test e_therm[i] ≈ internal_energy(model, p, T) / 4184. rtol=5e-4
            end
        end
    end
end
