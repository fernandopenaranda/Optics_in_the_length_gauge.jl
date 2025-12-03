"""
    σij. Drude conductivity (i = j symmetric). B independent
    semiclassical expression for the B-independent conductivity σ_ii 
    in the presence of a constant electric in the plane. 
    Output units: [e^2/h]
    Input units: [eV]
    τ is the scattering time in seconds
"""
drude_conductivity(p::Drude_presets) = drude_conductivity(p.h, p.dhi, p.T, p.τ, 
        p.computation.xbounds, p.computation.ybounds, p.computation.evals)

function drude_conductivity(h, dh, T, τ, xbounds, ybounds, evals; rel_tol = 1e-5, abs_tol = 0)
    C = 2π * τ * ħ_ev_s/ ang_to_m^2 #wait do we need a correction in here
    integrand(q) = real(k_in_plane_bindependent_conductivity(h(q), dh(q), T))
    val = bz_integration_transport(integrand, xbounds, ybounds, evals; rel_tol = rel_tol, abs_tol = abs_tol)
    return C * val
end
"""
    k_linear_magnetorresistance
integrand of σ_xx with all the prefactors
"""
function k_in_plane_bindependent_conductivity(h, dhi, T)  
    ϵs, ψs = eigen(Matrix(h))
    vi = vel(ψs, dhi) * ang_to_m/ ħ_ev_s
    sum(d_f(ϵs, 0, T) .* real(diag(vi)).^2)
end

"""
    k_linear_magnetorresistance
integrand of σ_ijk with all the prefactors
"""
function k_in_plane_bindependent_conductivity_integrand(ϵs, ψs, dhi, T)
    vi = vel(ψs, dhi) * ang_to_m/ħ_ev_s
    return sum(d_f(ϵs, 0, T) .* real(diag(vi)).^2)
end