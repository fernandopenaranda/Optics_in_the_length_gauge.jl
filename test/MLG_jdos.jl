using Optics_in_the_length_gauge: jdos

@testset "MLG_jdos" begin
    ham(q) = MLG_hamiltonian(0, 1, q)
    @test jdos(ham, [0,0], [1,1], [0.,1]; η = 0.1, evals = 100)
end


