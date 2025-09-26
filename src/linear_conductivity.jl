""" ωlist in meV
returns the interband linear optical conductivity
    \vec{dh} is \vec{\partial} H(\vec k)/\partial \vec{k} 

usage: 
    `ωlist, conds = σab_inter_linear(:x, :x, p, collect(0:1:20), evals = 4000)`
"""
function linear_optical_conductivity(a, b, h, xbounds, ybounds, ωlist; η = 0.1, evals = 10000, kws...)
    conds = zeros(Float64, length(ωlist))
    half_dim = length(ωlist)÷2
    conds[1:half_dim] .= integral_linear(ωlist[1:half_dim], a, b, h, dh, xbounds, ybounds, η, evals)
    conds[half_dim+1:length(ωlist)] .= integral_linear(ωlist[half_dim+1:length(ωlist)], a, b, h, dh, xbounds, ybounds, η, evals)
    return ωlist, conds
end

function integral_linear(ωlist::Array, a, b, h, dh, xbounds, ybounds, η, evals)
    integrand(q) = σab_linear_ω(ωlist, h(q), dh(q)[dir_to_ind(a)], dh(q)[dir_to_ind(b)], η)
    return bz_integration(integrand, xbounds, ybounds, ωlist, evals)
end

function σab_linear_ω(ωlist::Array, h, dha, dhb, η)
    ϵs, ψs = eigen(Matrix(h))
    mat = linear_integrand(ϵs, ψs, p, dha, dhb)
    return [π .* sum_nondiag(mat .* lorentz(ϵs, ω, η) .* -ω) for ω in ωlist] 
end

linear_integrand(ϵs, ψs, p, dha, dhb) = f(ϵs, p.μ, 0.0) .* (r(ϵs, ψs, dha) .* t(r(ϵs, ψs, dhb)))

""" 
adaptive integration for 2D
    bz_integration(f,xbounds, ybounds,  ωlist, evals) 
f integrand, xbounds, ybounds
"""
function bz_integration(f, xbounds, ybounds, ωlist, evals) 
    val, _ = Cubature.hcubature(length(ωlist), (x,v) -> v[:] = f(x), 
        [xbounds[1], ybounds[1]], [xbounds[2], ybounds[2]]; reltol = 1e-5, abstol=0, maxevals=evals);
    bz_surface = (1/(2pi*a0))^2 
    return bz_surface .* val
end

function dir_to_ind(dir::Symbol)::Int
    dir == :x ? 1 :
    dir == :y ? 2 :
    throw(ArgumentError("$dir direction is not implemented."))
end

flatten(v) = collect(Base.Iterators.flatten(v))