""" `dos(h, xbounds, ybounds, ωlist; kws...)`
Computes the density of states, provided a Hamiltonian function with momentum dependence h(\vec k),
BZ bounds along the x (xbounds) and y (ybounds) direction.
- `a` and `b` are the in-plane (:x, :y) directions.
- `η` defines the energy broadening of the Lorentzian peaks.
- `evals` sets the mesh size in the BZ integration
"""

dos(params::DOS_presets) =
    dos(params.h, params.computation.xbounds, params.computation.ybounds, 
    params.computation.ωlist, η = params.computation.broadening, evals =
    params.computation.evals)

function dos(h, xbounds, ybounds, ωlist; η = 0.1, evals = 10000, kws...)
    dos = zeros(Float64, length(ωlist))
    half_dim = length(ωlist)÷2
    dos[1:half_dim] .= integral_dos(ωlist[1:half_dim], h, xbounds, ybounds, η, evals)
    dos[half_dim+1:length(ωlist)] .= integral_dos(ωlist[half_dim+1:length(ωlist)], h, xbounds, ybounds, η, evals)
    return ωlist, dos
end

function integral_dos(ωlist::Array, h, xbounds, ybounds, η, evals)
    integrand(q) = dos_ω(ωlist, h, q, η)
    bz_vol = (1/(2pi))^(length(xbounds)) 
    return bz_vol .* bz_integration_optical(integrand, xbounds, ybounds, ωlist, evals)
end

function dos_ω(ωlist::Array, h, q, η)
    ϵs, = eigen(Matrix(h(q)))
    return [sum((1/π * η) ./ ((ω .- ϵs).^2 .+ η^2) ) for ω in ωlist]  
end

""" `jdos(h, xbounds, ybounds, ωlist; kws...)`
Note that the fermi level is passed on h, thus the functions do not depend on μ.
Computes the joint density of states, provided a Hamiltonian function with momentum dependence h(\vec k),
BZ bounds along the x (xbounds) and y (ybounds) direction.
- `a` and `b` are the in-plane (:x, :y) directions.
- `η` defines the energy broadening of the Lorentzian peaks.
- `evals` sets the mesh size in the BZ integration
"""
jdos(params::JDOS_presets) =
    jdos(params.h, params.computation.xbounds, params.computation.ybounds, 
    params.computation.ωlist, η = params.computation.broadening, evals =
    params.computation.evals)

function jdos(h, xbounds, ybounds, ωlist; η = 0.1, evals = 10000, kws...)
    jdos = zeros(Float64, length(ωlist))
    half_dim = length(ωlist)÷2
    jdos[1:half_dim] .= integral_jdos(ωlist[1:half_dim], h, xbounds, ybounds, η, evals)
    jdos[half_dim+1:length(ωlist)] .= integral_jdos(ωlist[half_dim+1:length(ωlist)], h, xbounds, ybounds, η, evals)
    return ωlist, jdos
end

function integral_jdos(ωlist::Array, h, xbounds, ybounds, η, evals)
    integrand(q) = jdos_ω(ωlist, h, q, η)
    bz_vol = (1/(2pi))^(length(xbounds)) 
    return bz_vol .* bz_integration_optical(integrand, xbounds, ybounds, ωlist, evals)
end

function jdos_ω(ωlist::Array, h, q, η)
    ϵs, = eigen(Matrix(h(q)))
    return [sum_nondiag(-f(ϵs, 0, 0) .* lorentz(ϵs, ω, η)) for ω in ωlist]   
end