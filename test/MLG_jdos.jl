using Optics_in_the_length_gauge: jdos

@testset "MLG_jdos" begin
    M = 2pi/3 # modulus of the M point
    xbounds =  [0,2M]
    ybounds =  [-1,1] .* 2π/√3
    ωlist = collect(0:2:10)
    ν = 1
    broadening = 1
    evals = 10
    dirJ = :x
    dirE = :x
    h(q) = MLG_hamiltonian(0, ν, q)
    dh(q) = MLG_hamiltonian_derivatives(ν, q)
    @test dos(h, xbounds, ybounds, ωlist; η = broadening, evals = evals)
    @test jdos(h, xbounds, ybounds, ωlist; η = broadening, evals = evals)
    @test linear_optical_conductivity(dirJ, dirE, h, dh, xbounds, ybounds, ωlist; η = broadening, evals = evals)
end