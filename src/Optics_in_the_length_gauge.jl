module Optics_in_the_length_gauge

    using Arpack
    using Cubature
    using ProgressMeter
    using Base.Threads
    using Distributed
    using Dierckx
    using PhysicalConstants
    using PhysicalConstants.CODATA2018
    using Unitful

    const k_B = (PhysicalConstants.CODATA2018.k_B |> u"eV/mK").val
    const ħ = PhysicalConstants.CODATA2018.ħ
    const e = PhysicalConstants.CODATA2018.e
    const C = ((e^3 / ħ^2) |> u"μA/V^2/s").val
    const C_cd = ((e^2/ħ) |> u"μA/V").val
    const ħ_ev_s = (ħ |> u"eV*s").val
    
    include("optics_operators.jl")
    include("jdos.jl")
    include("linear_conductivity.jl")
    
    export jdos, linear_optical_conductivity

end