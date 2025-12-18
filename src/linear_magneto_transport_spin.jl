#=
Computes the spin contribution to the linear magnetoconductivity coming from the spin d.o.f
its build upon a p::Planar_σijk_presets_orbital struct together with information of the spin 
matrix representations in the degrees of freddom of p.h
mm = omm + smm, with smm = -g μB ⟨[σ_x,σ_y,σ_z]⟩, g ≈ 2

_________________________________________________________________________________________
Remarks: current method is restricted to i = j = k = :x
=#
"""
    Hamiltonian in eV
    Units: [σxx]/[B]
    returns the SMM symmetric σxxx linear magneto conductivity conductivity
    `linear_magneto_conductivity(params::Planar_σijk_presets_orbital)``
    Input params see structs.jl:
    `dirJ::Symbol`, `dirE::Symbol`, directions of J and E (tensor component of the
     conductivity)`::Symbols -> {:x, :y}` unbounded 2D directions
    `h::Function` is a Hamiltonian function with dependence on momentum `q`
    `dh::Vector{Function}` is ∇h along the x and y directions (2D). An array of 
    functions `dh = [dh_x(q), dh_y(q)]` of momentum `q`
    `τ::Float64` scattering time
    `T::Float64` temperature
    `evals::Int` number of steps in the adaptive integration
    `rz` rz operator(q)
    details on the Hamiltonian (the chemical potential and dielectric field...) 
    are set in the parameters of the Hamiltonian stored in `Planar_σijk_presets_orbital`
"""
linear_magneto_conductivity_spin(params::Planar_σijk_presets_spin) =
    linear_magneto_conductivity_spin(params.s_op, params.p.a0, params.p.dirB, params.p.h, 
        params.p.nabla_nabla_h, params.p.τ, params.p.T, params.p.computation)

function linear_magneto_conductivity_spin(s_op, a0, k, h,  ddh, τ, T, cpt;
        rel_tol = 1e-5, abs_tol = 0)
    integrand(q) = k_linear_magneto_conductivity_spin(s_op[dir_to_ind(k)], h, ddh, q; 
        T = T, τ = τ) 
    integrator(observable) = 
        bz_integration_transport(observable, cpt, rel_tol = rel_tol, abs_tol = abs_tol)
    bz_vol = 1/(2pi*a0*ang_to_m)^length(cpt.xbounds)
    val = bz_vol * integrator(integrand)
    return val 
end 

# generalize to σyyy too
function k_linear_magneto_conductivity_spin(s_op, h, ddhi, q; T = 2, τ = 1e-15)
    ϵs, ψs = eigen(Matrix(h(q)))                                                                                                        
    return  2π * τ * k_linear_mr_integrand_spin(s_op, ϵs, ψs, ddhi(q), T)
end

"the term with vij is only valid for sigma xxx, generalized if required"
function k_linear_mr_integrand_spin(s_op, ϵs, ψs, dhxx, T)
        vxx = vel(ψs, dhxx) * ang_to_m^2/ ħ_ev_s # units meters, eV, seconds
        return real(sum(d_f(ϵs, 0, T) .* (- mr_vij_spin(s_op, ψs, vxx))))    
        #only valid in the xx direction
end

mr_vij_spin(s_op, ψs, vij) = real(SMM(s_op, ψs) .* diag(vij))
""" spin_magnetic moment """
SMM(s_op, ψs; g= 2) = - g * μB * ψs' * s_op * ψs

