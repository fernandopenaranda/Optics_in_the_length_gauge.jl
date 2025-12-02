module Optics_in_the_length_gauge
    using Arpack
    using LinearAlgebra
    using Cubature
    using ProgressMeter
    using Base.Threads
    using Distributed
    using Dierckx
    using PhysicalConstants
    using PhysicalConstants.CODATA2018
    using Unitful
    using SparseArrays
    using StaticArrays
    using Parameters
    const kB = (PhysicalConstants.CODATA2018.k_B |> u"eV/K").val
    const ħ = PhysicalConstants.CODATA2018.ħ
    const e = PhysicalConstants.CODATA2018.e
    const C = ((e^3 / ħ^2) |> u"μA/V^2/s").val
    const C_cd = ((e^2/ħ) |> u"μA/V").val
    const ħ_ev_s = (ħ |> u"eV*s").val
    const ang_to_m = 1e-10

    function does_not_throw(f::Function, args...)
        try
            f(args...)
            return true   # no error
        catch
            return false  # error occurred
        end
    end
    
    include("structs.jl")
    include("length_gauge_operators.jl")
    include("integration.jl")
    include("jdos.jl")
    include("linear_optical_conductivity.jl")
    include("drude_conductivity.jl")
    include("linear_magneto_transport.jl")
    #...
    export Optical_computation_presets, Transport_computation_presets, DOS_presets, JDOS_presets, σij_presets, Drude_presets, Planar_σijk_presets
    export does_not_throw, dos, jdos, linear_optical_conductivity, linear_magneto_conductivity, drude_conductivity
    
    # Export the presets submodule
    export Presets 
    include("presets/MLG_ham.jl") 
    using .MLGPresets  # Import the submodule
    module Presets # Create a nested submodule for the exported presets
        export MLG_hamiltonian, MLG_nabla, MLG_2deriv, K1, K2   # export functions in this submodule
        using ..MLGPresets: MLG_hamiltonian, MLG_nabla, MLG_2deriv, K1, K2
    end
end