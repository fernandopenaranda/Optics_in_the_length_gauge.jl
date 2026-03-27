"""
computes the σ_ij in the presence in a system finite in systems with an intrinsic Ωz ≠ 0
with i, j in the plane (x,y) in units of e^2/h
σ^A_ij = ϵ_ij 2π e^2/h ∑n ∫ dk^2/(2π)^d Ωz,n(k) f_n(k)
where f is the occupation function (temperature dependent),
and ϵ_ij = {1 if i = x, j = y, -1if i = y, j = x, 0 otherwise }
There are two methods written for practical convenience. Note although envisioned for quasi 2d,
    it is possible to use this function to compute sigma_xyz with z being periodic, see the more 
        general expression below
"""
σij_anomalous_hall(params::Planar_σijk_presets_orbital) = 
    σij_anomalous_hall(params.a0, params.dirJ, params.dirE,
    params.h, params.nabla_h, params.T, params.computation)


σij_anomalous_hall(params::AH_presets) = 
    σij_anomalous_hall(params.a0, params.dirJ, params.dirE,
    params.h, params.dh, params.T, params.computation)


function σij_anomalous_hall(a0, i, j, h, dh, T, cpt; rel_tol = 1e-5, abs_tol = 0) 
    warn_equalargs(i, j)
    integrand(q) = k_σij_anomalous_hall(i, j, h(q), dh(q), T)
    bz_vol = 1/(2pi*a0*ang_to_m)^length(cpt.xbounds)
    val = bz_integration_transport(integrand, cpt, 
        rel_tol = rel_tol, abs_tol = abs_tol)
    return -2π * bz_vol * val # in units of e^2/h
end

function k_σij_anomalous_hall(i,j,h,dh, T)
    ϵs, ψs = eigen(Matrix(h))     
    ri = r(ϵs, ψs, dh[which_ind(i)]) .* ang_to_m
    rj = r(ϵs, ψs, dh[which_ind(j)]) .* ang_to_m
    return sum([fn(ϵ, 0, T) for ϵ in ϵs] .* Ωz(i, ri, rj))
end

"""
    Ωz
computes the out-of-plane component of the Berry curvature given by
"""
Ωz(i, ri, rj) = -2ε(i) * imag(Σ_nondiag(ri, rj))
which_ind(i::Symbol) = i == :x ? 1 : 2
#_________________________________________________________________________________________

""" valid for 2d and 3d"""
function σij_anomalous_hall(p::AH_presets_3d) 
    warn_equalargs(p.dirJ, p.dirE)
    integrand(q) = k_σij_anomalous_hall_3d(p.dirJ, p.dirE, p.h(q), p.dh(q), p.T, p.gs)
    if length(p.gs) == 3
        VBZ = bz_volume(p.gs[1],p.gs[2],p.gs[3])
    else length(p.gs) == 2
        VBZ = bz_volume(p.gs[1],p.gs[2])
    end
    bz_vol = VBZ/(2pi*ang_to_m)^length(p.computation.xbounds) 
    val = bz_integration_transport_3d(integrand, p.computation, p.gs, rel_tol = 1e-5, abs_tol = 1e-7)
    return -2π * bz_vol * val * ang_to_m^2 # in units of e^2/h
end

function k_σij_anomalous_hall_3d(i,j,h,dh, T, gs)
    ϵs, ψs = eigen(Matrix(h))
    ωs = Ω(ϵs)
    if length(gs) == 3
        vels = [v(:x,ψs,dh), v(:y,ψs,dh), v(:z,ψs,dh)]
    else length(gs) == 2
        vels = [v(:x,ψs,dh), v(:y,ψs,dh)]
    end
    # return sum([fn(ϵ, 0, T) for ϵ in ϵs] .*Ωab(i,j, ωs,vels))
    return sum([fn(ϵ, 0, T) for ϵ in ϵs] .*Ωab(i,j, ωs,vels))

end


"""
expression for the Berry curvature assuming periodic boundary conditions along all axes.
`Ωz` in quantum-anomalous_hall is the z component of this formula that is, the only component,
valid in quasi 2d systems with z-bounded direction.
Ω^i_nn = i ħ^2 ∑_{m≠n} ϵ_ijk v^j_nm v^k_mn /ϵ_nm^2 
"""
function Ωi(i, ωs, vels) 
    s = zeros(ComplexF64, size(ωs,1))
    coords = [:x,:y,:z]
    for j in coords
        for k in coords
            s .+= levi_civita(i,j,k) .* Ωab(j,k, ωs ,vels) * 1/2
             # the 1/2 cancels the two below that upon the contraction of the two indices
        end
    end
    return s
end


""" note that this formula (2d and 3d) works because the two epsilons are contracted no need for levicivita tensors
(-2) comes out of the contraction ϵijk ϵkab = δia δjb - δib δja """
function  Ωab(j,k, ωs ,vels)
    vj = vels[symb_to_ind(j)]
    vk = vels[symb_to_ind(k)]
    return -2 * imag.(diag((vj .- diagm(diag(vj)))./(ωs.^2) *
        (vk .- diagm(diag(vk)))))
end