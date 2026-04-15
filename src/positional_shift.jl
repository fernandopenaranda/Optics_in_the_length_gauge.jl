# """
# classical (Lorentz) contribution to σijk^(A,1)
# """
# function classical_contribution()
# end

# function integrand_classical_contribution()
# end
# Important, all the 2D paths assume the z direction as bounded and x y to lay on the plane.
"""
quantum contribution to σijk^(A,1)
Implementation of the correction to the classical σ_ijk^(A,1)
including the effects of the positional shift (modification of the Berry connection
to second order in the EB fields) and the magnetic moment (MM) in the dispersion.
Note A implies ij are antisymmetric
The current implementaion assumes a bulk (3D) system, generalization yet to come. 
The magnetic moment has been written in general to accept both spin and 
orbital contributions. The keyword `which_mm` by default `:orbital` selects the 
orbital contribution only.

NOTE: This formalism assumes non-degenerate bands which is not true at the Weyl points

Units the integrand times the differential is in length units, because the velocities are 
in E * L, W = E * L^2, M = E * L^2
# parametrization k -> us; k = u1b1+u2b2+u3b3
    # note q is passed in u_i ∈ [-0.5,0.5], with k = sum_i u_i * b_i
the result is given in units of (S/m/T) 
Note that the units are ecube_hbarsquared because the velocity is expressed as [E][L] 
in the code insted of [L]/[T].

"""
quantum_contribution(p::Quantum_correction_σijk_antisym) = 
    quantum_contribution(p.a0, p.dirJ, p.dirE, p.dirB, p.h, p.nabla_h, 
        p.nabla_nabla_h, p.gs, p.τ, p.T, p.computation, p.which_mm, 
        p.Ω_MM_switch, p.PS_switch, p.PS_orbital_switch, p.QM_switch, p.fermi_surface, p.epsilon)

function quantum_contribution(a0, dirJ, dirE, dirB, h, dh, ddh, Gs, τ, T, cpt, which_mm, Ω_MM_switch, 
    PS_switch, PS_orbital_switch, QM_switch, fermi_surface, epsilon, rel_tol = 1e-5, abs_tol = 0)
    checkdims(cpt.xbounds)
    checkantisym(dirJ,dirE,dirB)
    VBZ = bz_volume(Gs)
    integrand(q) = integrand_quantum_contribution_q(dirJ, dirE, dirB, h, dh, ddh,
         T, q, Ω_MM_switch,  PS_switch, PS_orbital_switch, QM_switch, fermi_surface, epsilon, which_mm = which_mm)
    integrator(observable) = bz_integration_transport_3d(observable, cpt, Gs,
         rel_tol = rel_tol, abs_tol = abs_tol)
    bz_vol = VBZ/(2pi)^length(cpt.xbounds) # No a0 in the denominator
    val = bz_vol * integrator(integrand)#  Integral -) Å. Prefactor e^3/ħ^2 
    
    return  val * ecube_hbarsquared * ang_to_m #units of S/m/T
end

"there is also a contribution coming from the spin magnetic moment (path disabled for the moment)
f' * ∑_α (v_a F_{bc}^α - v_b F_{ac}^α + ϵ_{abd}Omega_d M_c^α).
Important not to confuse"

integrand_quantum_contribution_q(p, q) = integrand_quantum_contribution_q(p.dirJ, p.dirE, p.dirB, p.h, 
    p.nabla_h, p.nabla_nabla_h, p.T, collect(q), p.Ω_MM_switch, p.PS_switch, p.PS_orbital_switch,
    p.QM_switch, p.fermi_surface, p.epsilon, which_mm = p.which_mm)

function integrand_quantum_contribution_q(a, b, c, h, dh, ddh, T, q, Ω_MM_switch, 
    PS_switch, PS_orbital_switch, QM_switch, fermi_surface, epsilon; which_mm = :orbital)
    if length(q) == 3
        ϵs, ψs = eigen(Matrix(h(q)))   
        dhs = [dh(q)[1],dh(q)[2],dh(q)[3]]
        ddhs = [[ddh(q)[1][1],ddh(q)[1][2],ddh(q)[1][3]], 
                [ddh(q)[2][1],ddh(q)[2][2],ddh(q)[2][3]],
                [ddh(q)[3][1],ddh(q)[3][2],ddh(q)[3][3]]]
    else length(q) == 2
        ϵs, ψs = eigen(Matrix(h(q)))   
        dhs = [dh(q)[1],dh(q)[2]]
        ddhs = [[ddh(q)[1][1],ddh(q)[1][2]], 
                [ddh(q)[2][1],ddh(q)[2][2]]]
    end
    integrand_quantum_contribution(a, b, c, ϵs, ψs, dhs, ddhs, T, Ω_MM_switch, 
        PS_switch, PS_orbital_switch, QM_switch, which_mm, fermi_surface, epsilon)
end

function integrand_quantum_contribution(a, b, c, ϵs, ψs, dh, ddh, T, Ω_MM_switch, 
    PS_switch, PS_orbital_switch, QM_switch, which_mm, fermi_surface, ϵ = 1e-5)
    ϵ = kB*T
    ωs = Ω(ϵs) .+ 0im 
    ωs[real(ωs) .< 1e-4] .+= im*ϵ
    s = 0.0im
    if length(dh) == 3
        vels = [v(:x,ψs,dh), v(:y,ψs,dh), v(:z,ψs,dh)]  #units [E*L]
        vvels = d_3dvs(ψs, ddh)
    else length(dh) == 2 # this path assumes z to be bounded
        vels = [v(:x,ψs,dh), v(:y,ψs,dh)]  #units [E*L]
        vvels = d_2dvs(ψs, ddh)
    end
    df = d_f(ϵs, 0, T)
    if fermi_surface == false
        if PS_switch == true
            s += sum(df .* positional_shift(a, b, c, ωs, vels, vvels, PS_orbital_switch, QM_switch, which_mm = which_mm))
        end
        if Ω_MM_switch == true
            s += sum(df .* berry_OMM(a,b,c, ωs, vels, which_mm = which_mm))
        end
        return real(s)
        # return real(sum(d_f(ϵs, 0, T) .* (ifelse(PS_switch == true, 1, 0) .* 
        #     positional_shift(a, b, c, ωs, vels, vvels, QM_switch, which_mm = which_mm) .+ 
        #     ifelse(Ω_MM_switch == true, 1, 0) .* berry_OMM(a,b,c, ωs, vels, which_mm = which_mm))))
    else 
        # throw(ArgumentError(""))
        return real(sum(-d_f(ϵs, 0, T))) # DOS
    end
end

d_2dvs(ψs, ddh) = [[dv(:x, :x, ψs, ddh), dv(:x, :y, ψs, ddh)],
    [dv(:y, :x, ψs, ddh), dv(:y, :y, ψs, ddh)]] 

d_3dvs(ψs, ddh) = [[dv(:x, :x, ψs, ddh), dv(:x, :y, ψs, ddh), dv(:x, :z, ψs, ddh)],
    [dv(:y, :x, ψs, ddh), dv(:y, :y, ψs, ddh), dv(:y, :z, ψs, ddh)],
    [dv(:z, :x, ψs, ddh), dv(:z, :y, ψs, ddh), dv(:z, :z, ψs, ddh)]] 

"""contribution to the quantum correction coming from the product of the Berry curvature and 
the magnetic moment (orbital + spin). Spin deactivated.
Note that M in this expression is diagonal. That is it is the OMM of band n
The two expressed functions are equivalent. """
function berry_OMM(a,b,c, ωs, vels; which_mm = :orbital)
    s = zeros(ComplexF64, size(ωs,1))
    s .+= Ωab(a, b, ωs, vels) .* diag(interband_MM(c, ωs, vels; which_mm))
    return s
end

""" contribution to the quantum correction to σ_ijk^(A,1) coming from the positional 
shift contribution F """
function positional_shift(a, b, c, ωs, vels, vvels, PS_orbital_switch, QM_switch; which_mm = :orbital) 
    va = vels[symb_to_ind(a)]
    vb = vels[symb_to_ind(b)]
    return diag(va) .* F(c, b, ωs, vels, vvels, PS_orbital_switch, QM_switch, which_mm = which_mm) .- 
        diag(vb) .* F(c, a, ωs, vels, vvels, PS_orbital_switch, QM_switch, which_mm = which_mm) 
 end

function F(a, b, ωs, vels, vvels, PS_orbital_switch, QM_switch; which_mm = :orbital)   
    vb = vels[symb_to_ind(b)]
    nd_vb = copy(vb) - diagm(diag(vb))
    dim =  size(ωs,1)
    s = zeros(ComplexF64, dim)
    if PS_orbital_switch == true
        MM = interband_MM(a, ωs, vels, which_mm = which_mm)
        for n in 1:dim
            for m in 1:dim
                if n ≠ m
                    s[n] += -2*(-im * MM[n,m] * nd_vb[m,n] /(ωs[n,m]*ωs[m,n]))
                end
            end
        end
    end
    if QM_switch == true
        s .+= - 1/2 .* real(contracted_sum_qm(a, b, ωs, vels, vvels))
    end
    return s
end

"""
levi civita contracted. ωs are energy differences not to confuse with ϵs
"""
function contracted_sum_qm(a, b, ωs,vels, vvels; which_mm = :orbital)
    s = zeros(ComplexF64, size(ωs,1))
    for (c,d,sign) in allowed_components(a)
        s += sign .* qm_int(b, c, d, ωs,vels, vvels)
    end
    return s
end

function qm_int(a, c, d, ωs, vels, vvels)
    ωs_safe = copy(ωs) #.+ 1e-2*im
    va = vels[symb_to_ind(a)]
    vc = vels[symb_to_ind(c)]
    vd = vels[symb_to_ind(d)]
    Δd = diag(vd) .-  transpose(diag(vd))
    nd_va = copy(va) - diagm(diag(va))
    nd_vc = copy(vc) - diagm(diag(vc))
    nd_vd = copy(vd) - diagm(diag(vd))
    return diag((Δd .* nd_va ./ ωs_safe .^ 3) * nd_vc) + 
        diag((nd_vd ./ ωs_safe .^ 2)  * vvels[symb_to_ind(a)][symb_to_ind(c)])
end

"""
anomalous spin and orbital polarizability. Default only P see the definition of MM(which=:orbital)
dh = [∂x h, ∂y h, ∂z h]
ddh =  [[∂x∂x h, ∂x∂y h, ∂x∂z h], [∂y∂x h, ∂y∂y h, ∂y∂z h], [∂z∂x h, ∂z∂y h, ∂z∂z h]]
alternative return (slower) but equivalent:
    # return 2*real.(diag(va * (interband_MM(b, ωs, vels, which_mm = which_mm) ./ (ωs .^2))))  
        + 0 .* 1/2 * real(contracted_sum_qm(a, b, ωs, vels, vvels))
"""
function F(p::Quantum_correction_σijk_antisym, q,a,b)
    ϵs, ψs = eigen(Matrix(p.h(q)))   
    dhs = [p.nabla_h(q)[1],p.nabla_h(q)[2],p.nabla_h(q)[3]]
    ddhs = [[p.nabla_nabla_h(q)[1][1],p.nabla_nabla_h(q)[1][2],p.nabla_nabla_h(q)[1][3]], 
        [p.nabla_nabla_h(q)[2][1],p.nabla_nabla_h(q)[2][2],p.nabla_nabla_h(q)[2][3]],
        [p.nabla_nabla_h(q)[3][1],p.nabla_nabla_h(q)[3][2],p.nabla_nabla_h(q)[3][3]]]
    ϵ = kB*p.T
    ωs_epsilon = Ω(ϵs) .+ 0im 
    ωs_epsilon[real(ωs_epsilon) .< 1e-5] .+=  im*ϵ # this is to avoid divergences at band crossings.
    vels = [v(:x,ψs,dhs), v(:y,ψs,dhs), v(:z,ψs,dhs)]  #units [E*L]
    vvels = d_3dvs(ψs, ddhs) #units [E*L^2]
    return real.(F(a, b, ωs_epsilon, vels, vvels, p.PS_orbital_switch, p.QM_switch; which_mm = :orbital)) #.* real.(diag(vels[1])) # the commented is only for rapid access to vi fij do not consider it seriously
end

""" interband magnetic moment with orbital and spin parts """
function interband_MM(a,ωs, vels; which_mm = :orbital)
    if which_mm == :orbital
        return interband_OMM(a, ωs, vels) 
    elseif which_mm == :spin 
        return interband_SMM(a, ωs, vels)
    elseif which_mm == :both
        return interband_OMM(a, ωs, vels) + interband_SMM(a, ωs, vels)
    else
        return  throw(ArgumentError("which_contribution not in (:orbital,:spin,:both))"))
    end
end

interband_SMM(a, ωs, vels) = 0I #ignored for the moment

function levi_civita(i::Symbol, j::Symbol, k::Symbol)
    perm = (i, j, k)
    if length(unique(perm)) < 3
        return 0
    elseif perm in ((:x,:y,:z), (:y,:z,:x), (:z,:x,:y))
        return 1
    else
        return -1
    end
end
#_________________________________________________________________________________________

v(a::Symbol,ψs, dh) = v(symb_to_ind(a),ψs, dh)
v(a::Int,ψs, dh) = vel(ψs, dh[a])

dv(a::Symbol, b::Symbol, ψs, ddh) = dv(symb_to_ind(a), symb_to_ind(b), ψs, ddh)
dv(a::Int,b::Int, ψs, ddh) = vel(ψs, ddh[a][b])

function symb_to_ind(a::Symbol)
    if a == :x
        return 1
    elseif a == :y
        return 2
    elseif a == :z
        return 3
    else 
        return 0
    end
end

function checkdims(q) 
    if length(q) < 3
        throw(ArgumentError("q must be three dimensional- molecular limit not implemented"))
    else nothing end
end

function checkantisym(a,b,c)
    if a == b
        throw(ArgumentError("a and b cannot be the same coordinate, σ_abc should be 
        antisymmetric in a and b"))
    else nothing end
end

function transform_k(us, Gs)
    if length(Gs) == 3
        return us[1]*Gs[1] + us[2]*Gs[2] + us[3]*Gs[3]
    else
        return us[1]*Gs[1] + us[2]*Gs[2]
    end
end
    
""" levi civita tensor selection pairs"""
function allowed_components(a::Symbol)
    pairs = Dict(
    :x => [(:y,:z, 1), (:z,:y,-1)],
    :y => [(:z,:x, 1), (:x,:z,-1)],
    :z => [(:x,:y, 1), (:y,:x,-1)])
    return pairs[a]
end


#_________________________________________________________________________________________

# function contracted_sum_qm(a::Symbol, b::Symbol, c::Symbol)
#     s = 0.0
#     cords = [:x, :y, :z]
#     d = setdiff(cords, (b,c))[1]  # only one element remains
#     if (b,c,d) in [(:x,:y,:z), (:y,:z,:x), (:z,:x,:y)]
#         sign = 1
#     elseif (b,c,d) in [(:x,:z,:y), (:z,:y,:x), (:y,:x,:z)]
#         sign = -1
#     end
#     s += sign * partial_quantum_metric(c, d, a)
# end

# function F_2(p::Quantum_correction_σijk_antisym, q,a,b)
#     ϵs, ψs = eigen(Matrix(p.h(q)))   
#     dhs = [p.nabla_h(q)[1],p.nabla_h(q)[2],p.nabla_h(q)[3]]
#     ddhs = [[p.nabla_nabla_h(q)[1][1],p.nabla_nabla_h(q)[1][2],p.nabla_nabla_h(q)[1][3]], 
#         [p.nabla_nabla_h(q)[2][1],p.nabla_nabla_h(q)[2][2],p.nabla_nabla_h(q)[2][3]],
#         [p.nabla_nabla_h(q)[3][1],p.nabla_nabla_h(q)[3][2],p.nabla_nabla_h(q)[3][3]]]
#     ϵ = kB*p.T
#     ωs_safe = Ω(ϵs)  .+ 0.0im
#     ωs_safe[real(ωs_safe) .< 1e-4] .+= im*ϵ
#     vels = [v(:x,ψs,dhs), v(:y,ψs,dhs), v(:z,ψs,dhs)]  #units [E*L]
#     vvels = d_3dvs(ψs, ddhs) #units [E*L^2]
#     return real.(F_2(a, b, ωs_safe, vels, vvels, p.QM_switch; which_mm = :orbital))
# end

 # note that the indices of F have been altered

# function F_2(a, b, ωs, vels, vvels, QM_switch; which_mm = :orbital)   
#     s = zeros(ComplexF64, size(ωs,1))
#     for (c,d,sign) in allowed_components(a)
#         s += -sign * aux_F(b, c, d, ωs, vels, vvels, QM_switch)
#     end
#     return s
# end

# function aux_F(b, c, d, ωs,vels, vvels, QM_switch)
#     if QM_switch == true
#         switch = 1
#     else switch = 0 end
#     vb = vels[symb_to_ind(b)]
#     vc = vels[symb_to_ind(c)]
#     vd = vels[symb_to_ind(d)]
#     term1 = zeros(ComplexF64,size(vb,1))
#     term2 = zeros(ComplexF64,size(vb,1))
#     term3 = zeros(ComplexF64,size(vb,1))

#     for n in 1:size(vb,1)
#         for m in 1:size(vb,1)
#             if n ≠ m
#                 term1[n] +=  (2vc[n,n] + vc[m,m] + switch*(-vc[n,n] + vc[m,m]))/2 * (vd[n,m] * vb[m,n])/ ωs[n,m]^3
#                 term2[n] +=  switch*vd[n,m] * vvels[symb_to_ind(b)][symb_to_ind(c)][m,n] / (2*ωs[n,m]^2)
#                 for l in 1:size(vb,1)
#                     if l ≠ n && l ≠ m
#                         term3[n] += 0*vc[n,l] * vd[l,m] * vb[m,n]/ (ωs[l,m]* ωs[n,m]^2)
#                     end
#                 end
#             else nothing end
#         end
#     end
#     return term1 + term2 + term3
# end