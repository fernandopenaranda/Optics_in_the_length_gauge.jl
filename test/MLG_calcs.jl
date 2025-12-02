function does_not_throw(f::Function, args...)
    try
        f(args...)
        return true   # no error
    catch
        return false  # error occurred
    end
end

@testset "MLG_calcs" begin
    cpo = Optical_computation_presets(
        xbounds = [0, 2π/√3],
        ybounds = [-2π, 2π],
        ωlist = collect(0:1:2),
        broadening = 0.01,
        evals = 1e1)
    cpt = Transport_computation_presets(
        xbounds = [0, 2π/√3],
        ybounds = [-2π, 2π],
        evals = 1e1)

    μ = 0
    h(q) = Presets.MLG_hamiltonian(μ, q) # system specific
    nabla_h(q) = Presets.MLG_nabla(q) 
    # spectral tests
    dos_presets = DOS_presets(h = h, computation = cpo)
    jdos_presets = JDOS_presets(h, cpo)
    #optical tests
    sigma_ij_presets = σij_presets(:x, :x, h, nabla_h, cpo)
    #transport tests
    drude_xx_presets = Drude_presets(:x,:x,h,nabla_h,10,cpt)
    end
    @test does_not_throw(dos, dos_presets)
    @test does_not_throw(jdos, jdos_presets)
    @test does_not_throw(linear_optical_conductivity, sigma_ij_presets)
    @test does_not_throw(drude_conductivity, drude_xx_presets)
end