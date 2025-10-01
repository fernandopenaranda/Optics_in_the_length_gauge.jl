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


    const k_B = (PhysicalConstants.CODATA2018.k_B |> u"eV/mK").val
    const ħ = PhysicalConstants.CODATA2018.ħ
    const e = PhysicalConstants.CODATA2018.e
    const C = ((e^3 / ħ^2) |> u"μA/V^2/s").val
    const C_cd = ((e^2/ħ) |> u"μA/V").val
    const ħ_ev_s = (ħ |> u"eV*s").val
    
    include("structs.jl")
    include("optics_operators.jl")
    include("integration.jl")
    include("jdos.jl")
    include("linear_conductivity.jl")

    
    # Presets
    # Export the presets submodule
    export Presets

    # Include the presets file(s)
    include("presets/MLG_ham.jl")

    # Import the submodule
    using .MLGPresets

    # Create a nested submodule for the exported presets
    module Presets
        export MLG_hamiltonian, MLG_nabla   # export functions in this submodule
        using ..MLGPresets: MLG_hamiltonian, MLG_nabla
    end



    export Computation_presets, DOS_presets, JDOS_presets, σij_presets
    export dos, jdos, linear_optical_conductivity
end