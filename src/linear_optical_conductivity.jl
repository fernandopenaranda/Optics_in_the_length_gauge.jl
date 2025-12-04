"""
    returns the real (absortive) part of the linear optical conductivity
    `σab_inter_linear(dirE, dirJ, h, dh, ωlist, η = broadening, evals = evals)``
    Input: 
    `dirJ::Symbol`, `dirE::Symbol`, directions of J and E (tensor component of the conductivity)`::Symbols -> {:x, :y}` unbounded 2D directions
    `h::Function` is a Hamiltonian function with dependence on momentum `q`
    `dh::Vector{Function}` is ∇h along the x and y directions (2D). An array of functions `dh = [dh_x(q), dh_y(q)]` of momentum `q`
    `broadening::Float64` energy scale 
    `evals::Int` number of steps in the adaptive integration
    `ωlist::Vector` frequency list of the incoming radiation
"""
linear_optical_conductivity(params::σij_presets) =
    linear_optical_conductivity(params.a0, params.dirJ, params.dirE, params.h, params.nabla_h,
    params.computation.xbounds, params.computation.ybounds, 
    params.computation.ωlist, η = params.computation.broadening, 
    evals = params.computation.evals)

function linear_optical_conductivity(a0, dirJ::Symbol, dirE::Symbol, h::Function, dh, xbounds, ybounds, 
    ωlist; η = 0.1, evals = 10000, kws...)
    conds = zeros(Float64, length(ωlist))
    half_dim = length(ωlist)÷2
    conds[1:half_dim] .= integral_linear(ωlist[1:half_dim], 
        dirJ, dirE, h, dh, xbounds, ybounds, η, evals, a0)
    conds[half_dim+1:length(ωlist)] .= integral_linear(ωlist[half_dim+1:length(ωlist)],
        dirJ, dirE, h, dh, xbounds, ybounds, η, evals, a0)
    return ωlist, π * conds # units of e^2/(16ħ)
end

function integral_linear(ωlist::Array, dirJ, dirE, h, dh, xbounds, ybounds, η, evals, a0)
    integrand(q) = real(σab_linear_ω(ωlist, h(q), dh(q)[dir_to_ind(dirJ)], dh(q)[dir_to_ind(dirE)], η))
    bz_vol = (1/(2pi*a0*ang_to_m))^(length(xbounds)) 
    return bz_vol .* bz_integration_optical(integrand, xbounds, ybounds, ωlist, evals)
end

function σab_linear_ω(ωlist::Array, h, dh_dirJ, dh_dirE, η)
    ϵs, ψs = eigen(Matrix(h))
    mat = linear_integrand(ϵs, ψs, dh_dirJ, dh_dirE)
    return [sum_nondiag(mat .* lorentz(ϵs, ω, η) .* -ω) for ω in ωlist] 
end

linear_integrand(ϵs, ψs, dha, dhb) = f(ϵs, 0.0, 0.0) .* (r(ϵs, ψs, dha) .* t(r(ϵs, ψs, dhb)))

function dir_to_ind(dir::Symbol)::Int
    dir == :x ? 1 :
    dir == :y ? 2 :
    throw(ArgumentError("$dir direction is not implemented."))
end

flatten(v) = collect(Base.Iterators.flatten(v))