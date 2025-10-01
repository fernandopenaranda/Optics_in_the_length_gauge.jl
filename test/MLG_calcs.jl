@testset "MLG_calcs" begin
    cp = Computation_presets(
        xbounds = [0, 2π/√3],
        ybounds = [-2π, 2π],
        ωlist = collect(0:1:2),
        broadening = 0.01,
        evals = 1e1
    )
    
    μ = 0
    h(q) = Presets.MLG_hamiltonian(μ, q) # system specific
    nabla_h(q) = Presets.MLG_nabla(q) 
    dos_presets = DOS_presets(h = h, computation = cp)
    jdos_presets = JDOS_presets(h, cp)
    sigma_ij_presets = σij_presets(:x, :x, h,
        nabla_h, ncp)

    @test dos(dos_presets)
    @test jdos(jdos_presets)
    @test linear_optical_conductivity(sigma_ij_presets)
end