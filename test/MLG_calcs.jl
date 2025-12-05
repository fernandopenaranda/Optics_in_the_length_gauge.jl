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
    a0 = 2.46 #in Å
    h(q) = Presets.MLG_hamiltonian(μ, q) # system specific
    nabla_h(q) = Presets.MLG_nabla(q)
    dhx(q) = Presets.MLG_nabla(q)[1]
    dhy(q) = Presets.MLG_nabla(q)[2]
    dhxx(q) = Presets.MLG_2deriv(q,:x)
    rz_mat(q,ψs) = ψs'*ψs #no q dependence
    # spectral tests
    dos_presets = DOS_presets(a0, h, cpo)
    jdos_presets = JDOS_presets(a0, h, cpo)
    #optical tests
    sigma_ij_presets = σij_presets(a0, :x, :x, h, nabla_h, cpo)
    #transport tests
    drude_xx_presets = Drude_presets(a0, :x,:x,h,dhx,10,200e-15,cpt)
    ahe_xy_presets = AH_presets(a0, :x, :y, h, nabla_h, 2, cpt)
    planar_σijk_presets = Planar_σijk_presets(a0, :x,:x,:x, h,nabla_h, dhxx, rz_mat, 200e-15, 10, cpt, true, true, false, true)
    @test does_not_throw(dos, dos_presets)
    @test does_not_throw(jdos, jdos_presets)
    @test does_not_throw(linear_optical_conductivity, sigma_ij_presets)
    @test does_not_throw(drude_conductivity, drude_xx_presets)
    @test does_not_throw(σij_anomalous_hall, ahe_xy_presets)
    @test does_not_throw(linear_magneto_conductivity, planar_σijk_presets)
end