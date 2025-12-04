
"""
computes the σ_ij in the presence in a system finite in systems with an intrinsic Ωz ≠ 0
with i, j in the plane (x,y) in units of e^2/h
σ^A_ij = ϵ_ij 2π e^2/h ∑n ∫ dk^2/(2π)^d Ωz,n(k) f_n(k)
where f is the occupation function (temperature dependent),
and ϵ_ij = {1 if i = x, j = y, -1if i = y, j = x, 0 otherwise }

"""
σij_anomalous_hall(params::Planar_σijk_presets) = 
    σij_anomalous_hall(params.a0, params.dirJ, params.dirE,
    params.h, params.dh, params.T, params.computation)


σij_anomalous_hall(params::AH_presets) = 
    σij_anomalous_hall(params.a0, params.dirJ, params.dirE,
    params.h, params.dh, params.T, params.computation)


function σij_anomalous_hall(a0, i, j, h, dh, T, cpt; rel_tol = 1e-5, abs_tol = 0) 
    warn_equalargs(i, j)
    integrand(q) = k_σij_anomalous_hall(i, j, h(q), dh(q), T)
    bz_vol = 1/(2pi*a0*ang_to_m)^length(cpt.xbounds)
    val = bz_integration_transport(integrand, cpt, 
        rel_tol = rel_tol, abs_tol = abs_tol)
    return 2π * bz_vol * val # in units of e^2/h
end

function k_σij_anomalous_hall(i,j,h,dh, T)
    ϵs, ψs = eigen(Matrix(h))     
    ri = r(ϵs, ψs, dh[which_ind(i)]) .* ang_to_m
    rj = r(ϵs, ψs, dh[which_ind(j)]) .* ang_to_m
    return sum(f(ϵs, 0, T) .* Ωz(i, ri, rj))
end
"""
    Ωz
computes the out-of-plane component of the Berry curvature given by
"""
Ωz(i, ri, rj) = 2ε(i) * imag(Σ_nondiag(ri, rj))

which_ind(i::Symbol) = i == :x ? 1 : 2
