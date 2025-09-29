""" `jdos(h, xbounds, ybounds, ωlist; kws...)`
Computes the joint density of states, provided a Hamiltonian function with momentum dependence h(\vec k),
BZ bounds along the x (xbounds) and y (ybounds) direction.
- `a` and `b` are the in-plane (:x, :y) directions.
- `η` defines the energy broadening of the Lorentzian peaks.
- `evals` sets the mesh size in the BZ integration
"""
function jdos(h, xbounds, ybounds, ωlist; η = 0.1, evals = 10000, kws...)
    jdos = zeros(Float64, length(ωlist))
    half_dim = length(ωlist)÷2
    jdos[1:half_dim] .= integral_jdos(ωlist[1:half_dim], h, xbounds, ybounds, η, evals)
    jdos[half_dim+1:length(ωlist)] .= integral_jdos(ωlist[half_dim+1:length(ωlist)], h, xbounds, ybounds, η, evals)
    return ωlist, jdos
end

function integral_jdos(ωlist::Array, h, xbounds, ybounds, η, evals)
    integrand(q) = jdos_ω(ωlist, h, q, η)
    return bz_integration(integrand, xbounds, ybounds, ωlist, evals)
end

function jdos_ω(ωlist::Array, h, q, η)
    ϵs, = eigen(Matrix(h(q)))
    return [sum_nondiag(-f(ϵs, 0, 0) .* lorentz(ϵs, ω, η)) for ω in ωlist]   
end