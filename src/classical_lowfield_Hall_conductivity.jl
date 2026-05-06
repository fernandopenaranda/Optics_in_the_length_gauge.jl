#=Classical (Lorentz) contribution to $\sigma_{ijk} = \frac{-e^3}{\hbar} 
\int_k f_0' \tau v_i \epsilon_{lmk} v_m (\partial_l \tau v_j + \tau \partial_l v_j)
 -e^3 \int_k f_0'' \tau^2 v_i \epsilon_{lmk} v_m v_l v_j $ with $E_j$ and $B_k$
 with j_i = \sigma_ijk Ej Bk. This tensor component will be employed 
 to compute the second order response in B_k Bα, being α the spin index.
 Remarks the tensor is already antysymmetric in i and j. 
    So it is not a longitudinal but a Hall response
 =#

 classical_contribution_sigmaijk(p::Quantum_correction_σijk_antisym) = 
    classical_contribution_sigmaijk(p.dirJ, p.dirE, p.dirB, p.h, p.nabla_h, 
    p.nabla_nabla_h, p.gs, p.τ, p.T, p.computation, p.fermi_surface)
    
 classical_contribution_sigmaijk(p::Classical_σijk_antisym) = 
    classical_contribution_sigmaijk(p.dirJ, p.dirE, p.dirB, p.h, p.nabla_h, 
     p.nabla_nabla_h, p.gs, p.τ, p.T, p.computation, p.fermi_surface)


function classical_contribution_sigmaijk(dirJ, dirE, dirB, h, dh, ddh, Gs, 
    τ, T, cpt, fermi_surface, rel_tol = 1e-5, abs_tol = 0)
    checkperiodicity(dirJ, dirE, dirB, cpt)
    # checkdims(cpt.xbounds)
    checkantisym(dirJ,dirE,dirB)
    integrand(q) = integrand_classical_contribution_sigmaijk_q(dirJ, dirE, dirB, h, dh, ddh,T, q, fermi_surface)
    integrator(observable) = bz_integration_transport_3d(observable, cpt, Gs, rel_tol = rel_tol, abs_tol = abs_tol)
    bz_vol = bz_volume(Gs)/(2pi)^length(cpt.xbounds)             # No a0 in the denominator
    val = bz_vol * integrator(integrand)  # This has units of eV^2s^2*L we need this to have units of L/TE in 3d, so three hbars dividing # yielding the e^3/hbar^4
    return  val * ecube_hbarfour/2* ang_to_m^(abs(length(cpt.xbounds)-4)) * τ^2 * fs_to_s^2 #units of S/m/T in 3d, S/T in 2d, tau must be in seconds
end

integrand_classical_contribution_sigmaijk_q(p::Classical_σijk_antisym, q) = 
    integrand_classical_contribution_sigmaijk_q(p.dirJ, p.dirE, p.dirB, p.h, p.nabla_h, p.nabla_nabla_h, p.T, q, p.fermi_surface)

function integrand_classical_contribution_sigmaijk_q(dirJ, dirE, dirB, h, dh, ddh,T, q, fermi_surface)
    s = zeros(ComplexF64, 4)
    icc_eval!(s, dirJ, dirE, dirB, h, dh, ddh,T, q, fermi_surface)
    return real(sum(s))
end

function icc_eval!(s, dirJ, dirE, dirB, h, dh, ddh,T, q, fermi_surface)
    ϵs, ψs = eigen(Matrix(h(q)))
    ϵ = kB*T
    ωs = Ω(ϵs) .+ 0im 
    ωs[real(ωs) .< 1e-4] .+= im*ϵ
    df = d_f(ϵs, 0, T)
    vels, vvels = set_vels(q, h, dh, ddh, ψs)
    if fermi_surface == true
        s .+= -df
    else 
        # s .+= df .* classical_contribution_q(dirJ, dirE, dirB,ϵs,vels,vvels)
        s .+= classical_contribution_q(dirJ, dirE, dirB,ϵs,vels,vvels)

    end
end

spin_pol_integrand_classical_contribution_sigmaijk_q(p::Classical_σijk_antisym, q) = 
spin_pol_integrand_classical_contribution_sigmaijk_q(p.dirJ, p.dirE, p.dirB, p.h, p.nabla_h, p.nabla_nabla_h, p.T, q, p.fermi_surface)

function spin_pol_integrand_classical_contribution_sigmaijk_q(dirJ, dirE, dirB, h, dh, ddh,T, q, fermi_surface)
    s = zeros(ComplexF64, 4)
    icc_eval!(s, dirJ, dirE, dirB, h, dh, ddh,T, q, fermi_surface)
    ϵs, ψs = eigen(Matrix(h(q)))
    spinprojection = [real(ψs[:, i]' * diagm([1, 1, -1, -1]) * ψs[:, i]) for i in 1:4]
    return real(sum(s .* spinprojection))
end

function set_vels(q, h, dh, ddh, ψs)
    if length(q) == 3
        dhs = [dh(q)[1],dh(q)[2],dh(q)[3]]
        ddhs = [[ddh(q)[1][1],ddh(q)[1][2],ddh(q)[1][3]], 
        [ddh(q)[2][1],ddh(q)[2][2],ddh(q)[2][3]],
        [ddh(q)[3][1],ddh(q)[3][2],ddh(q)[3][3]]]
        vels = [v(:x,ψs,dhs), v(:y,ψs,dhs), v(:z,ψs,dhs)]  #units [E*L]
        vvels = d_3dvs(ψs, ddhs)
    else length(q) == 2
        ϵs, ψs = eigen(Matrix(h(q)))   
        dhs = [dh(q)[1],dh(q)[2]]
        auxmat = zeros(size(dhs[1],1), size(dhs[1],2))
        ddhs = [[ddh(q)[1][1],ddh(q)[1][2], auxmat], 
                [ddh(q)[2][1],ddh(q)[2][2], auxmat], 
                [auxmat,auxmat,auxmat]]
        vels = [v(:x,ψs,dhs), v(:y,ψs,dhs), auxmat]  #units [E*L]
        vvels = d_3dvs(ψs, ddhs)
    end
    return vels, vvels
end

function classical_contribution_q(a,b,c,ϵs,vs,vvs)
    s = zeros(ComplexF64, size(ϵs,1))
    for (l,m,sign) in allowed_components(c)
        s .+= sign * ccq(symb_to_ind(a),symb_to_ind(b),symb_to_ind(l),symb_to_ind(m), vs,vvs)
    end
    return s
end

ccq(a,b,l,m,vs,vvs) = diag(vs[m]) .* (diag(vvs[b][m]) .* diag(vs[a]) .- diag(vvs[a][l]) .* diag(vs[b]))

function checkperiodicity(dirJ, dirE, dirB, cpt)
    if length(cpt.xbounds) < 3
        @warn "Molecular limit expressions not yet implemented. All z velocities are zero 2d limit. Not suited for quasi-2d materials"
    end
end