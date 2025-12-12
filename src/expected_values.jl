"""
    `filling(h, μ, Δx, T; kws...)`
returns the band filling over the integration area set by Δx
    `filling(h, μ, Δx, Δy, T; kws...)`
returns the band filling over the integration area set by Δx, Δy
note that if the integration area is not the full BZ, as it happens
usually when dealing with continuum models, a postprocessing regularization
by a momentum cutoff will be required. Since this postprocessing is model
dependent no regularization path is included within `Optics_in_the_length_gauge`.
"""
filling(h, μ, Δx, T; kws...) = expected_value(h, μ, 1I, Δx, T; kws...)
filling(h, μ, Δx, Δy, T; kws...) = expected_value(h, μ, 1I, Δx, Δy, T; kws...)
"""
    `expected_value(h, μ, op, Δx, T; kws...)`
returns the expected value of a given operator (same dims as h)
over the integration area set by Δx
    ` expected_value(h, μ, op, Δx, Δy, T; evals = 100)`
returns the expected value of a given operator (same dims as h)
    over the integration area set by Δx , Δy
note that if the integration area is not the full BZ, as it happens
usually when dealing with continuum models, a postprocessing regularization
by a momentum cutoff will be required. Since this postprocessing is model
dependent no regularization path is included within `Optics_in_the_length_gauge`.
"""
expected_value(h, μ, op, Δx, T; kws...) = expected_value(h, μ, op, Δx, Δx, T; kws...)

function expected_value(h, μ, op, Δx, Δy, T; evals = 100, rel_tol= 1e-5, abs_tol = 0)
    integrand(q) = k_expected_value(h, μ, op, q, T)
    val = bz_integration_transport(integrand, Δx , Δy, evals; rel_tol = rel_tol, abs_tol = abs_tol) 
    return val /((Δx[2]-Δx[1])*(Δy[2]-Δy[1]))
end 

function k_expected_value(h, μ, op, q, T)
    check_same_size(h(q), op)
    ϵs, ψs = eigen(h(q))
    if  T == 0
        inds  = findall(ϵ-> ϵ < μ, real(ϵs))
        s = 0.0 
        for ind in inds
            s += ψs[:,ind]'* op * ψs[:,ind]
        end
    else 
        s = 0.0 
        for (i, ϵ) in ϵs
            s += ψs[:,i]'* op * ψs[:,i] * fn(real(ϵ), μ, T)
        end
    end

    return real(s)
end

check_same_size(A::AbstractMatrix, B::UniformScaling) = nothing
function check_same_size(A::AbstractMatrix, B::AbstractMatrix)
    if size(A) != size(B)
        error("Matrix dimensions do not match: size(A) = $(size(A)), size(B) = $(size(B))")
    end
end