# """
# classical (Lorentz) contribution to σijk^(A,1)
# """
# function classical_contribution()
# end

# function integrand_classical_contribution()
# end
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

the result is given in units of (S/mT) 
val * 2e^2/h * pi * e/ħ = val * g0 * pi*e/ħ
"""
quantum_contribution(p::Quantum_correction_σijk_antisym) = 
    quantum_contribution(p.a0, p.dirJ, p.dirE, p.dirB, p.h, p.nabla_h, p.nabla_nabla_h, p.gs, p.τ, 
        p.T, p.computation, p.which_mm, p.Ω_MM_switch, p.PS_switch, p.QM_switch, p.fermi_surface, p.epsilon)

function quantum_contribution(a0, dirJ, dirE, dirB, h, dh, ddh, Gs, τ, T, cpt, which_mm, Ω_MM_switch, 
    PS_switch, QM_switch, fermi_surface, epsilon, rel_tol = 1e-5, abs_tol = 0)
    checkdims(cpt.xbounds)
    checkantisym(dirJ,dirE,dirB)
    VBZ = bz_volume(Gs[1],Gs[2],Gs[3]) * cube_volume(cpt.xbounds, cpt.ybounds)
    # parametrization k -> us; k = u1b1+u2b2+u3b3
    # note q is passed in u_i ∈ [-0.5,0.5], with k = sum_i u_i * b_i
    integrand(q) = integrand_quantum_contribution_q(dirJ, dirE, dirB, h, dh, ddh, T, transform_k(q, Gs) , Ω_MM_switch, 
        PS_switch, QM_switch, fermi_surface, epsilon, which_mm = which_mm)
    integrator(observable) = bz_integration_transport_3d(observable, cpt, rel_tol = rel_tol, abs_tol = abs_tol)
    bz_vol = VBZ/(2pi)^length(cpt.xbounds) #if q had no units you would need to divide by a0
    val = bz_vol * integrator(integrand)# val has units of Å. the result integrand is in Angstroms times e^3/ħ^2 
    return  val * g0 * π * e_o_ħ_m * ang_to_m  #units of S/(mT)
end
"there is also a contribution coming from the spin magnetic moment (path disabled for the moment)
f' * ∑_α (v_a F_{bc}^α - v_b F_{ac}^α + ϵ_{abd}Omega_d M_c^α).
Important not to confuse"
function integrand_quantum_contribution_q(a, b, c, h, dh, ddh, T, q::Array, Ω_MM_switch, 
    PS_switch, QM_switch, fermi_surface, epsilon; which_mm = :orbital)
    ϵs, ψs = eigen(Matrix(h(q)))   
    dhs = [dh(q)[1],dh(q)[2],dh(q)[3]]
    ddhs = [[ddh(q)[1][1],ddh(q)[1][2],ddh(q)[1][3]], 
            [ddh(q)[2][1],ddh(q)[2][2],ddh(q)[2][3]],
            [ddh(q)[3][1],ddh(q)[3][2],ddh(q)[3][3]]]
    integrand_quantum_contribution(a, b, c, ϵs, ψs, dhs, ddhs, T, Ω_MM_switch, PS_switch, QM_switch, which_mm, fermi_surface, epsilon)
end

function integrand_quantum_contribution(a, b, c, ϵs, ψs, dh, ddh, T, Ω_MM_switch, 
    PS_switch, QM_switch, which_mm, fermi_surface, ϵ = 1e-5)
    ωs = Ω(ϵs)
    ωs .+= (findmax(abs.(ωs))[1]) * ϵ # this is to avoid divergences at band crossings.
    vels = [v(:x,ψs,dh), v(:y,ψs,dh), v(:z,ψs,dh)]  #units [E*L]
    vvels = [[dv(:x, :x, ψs, ddh), dv(:x, :y, ψs, ddh), dv(:x, :z, ψs, ddh)],
             [dv(:y, :x, ψs, ddh), dv(:y, :y, ψs, ddh), dv(:y, :z, ψs, ddh)],
             [dv(:z, :x, ψs, ddh), dv(:z, :y, ψs, ddh), dv(:z, :z, ψs, ddh)]] #units [E*L^2]
    if fermi_surface == false
        return real(sum(d_f(ϵs, 0, T) .* (ifelse(PS_switch == true, 1, 0) .* 
            positional_shift(a, b, c, ωs, vels, vvels, QM_switch, which_mm = which_mm) .+ 
            ifelse(Ω_MM_switch == true, 1, 0) .* berry_OMM(a,b,c, ωs, vels, which_mm = which_mm))))
    else 
        return real(sum(d_f(ϵs, 0, T))) # DOS
    end
end

"""contribution to the quantum correction coming from the product of the Berry curvature and 
the magnetic moment (orbital + spin). Spin deactivated.
Note that M in this expression is diagonal. That is it is the OMM of band n"""
function berry_OMM(a,b,c, ωs, vels; which_mm = :orbital)
    s = zeros(ComplexF64, size(ωs,1))
    for d in [:x,:y,:z]
        s .+= ϵ(a,b,d) .* Ωi(d, ωs, vels) .* diag(interband_MM(c, ωs, vels, which_mm = which_mm))
    end
    return s
end
"""
expression for the Berry curvature assuming periodic boundary conditions along all axes.
`Ωz` in quantum-anomalous_hall is the z component of this formula that is, the only component,
valid in quasi 2d systems with z-bounded direction.
Ω^i_nn = i ħ^2 ∑_{m≠n} ϵ_ijk v^j_nm v^k_mn /ϵ_nm 
"""
function Ωi(i, ωs, vels) 
    s = zeros(ComplexF64, size(ωs,1))
    coords = [:x,:y,:z]
    ωs .+= diagm(ones(size(ωs,1))) 
    for j in coords
        for k in coords
            s .+= ϵ(i,j,k) .* Ωi_aux(j,k, ωs ,vels)
        end
    end
    return 1im .* s
end

function  Ωi_aux(j,k, ωs ,vels)
    vj = vels[symb_to_ind(j)]
    vk = vels[symb_to_ind(k)]
    return diag((vj .- diagm(diag(vj))./(ωs.^2)) *(vk .- diagm(diag(vk)))) 
end

""" contribution to the quantum correction to σ_ijk^(A,1) coming from the positional 
shift contribution F """
function positional_shift(a, b, c, ωs, vels, vvels, QM_switch; which_mm = :orbital) 
    va = vels[symb_to_ind(a)]
    vb = vels[symb_to_ind(b)]
    return diag(va) .* F(b,c, ωs, vels, vvels, QM_switch, which_mm = which_mm) .- 
        diag(vb) .* F(a, c, ωs, vels, vvels, QM_switch, which_mm = which_mm) 
end

"""
anomalous spin and orbital polarizability. Default only P see the definition of MM(which=:orbital)
dh = [∂x h, ∂y h, ∂z h]
ddh =  [[∂x∂x h, ∂x∂y h, ∂x∂z h], [∂y∂x h, ∂y∂y h, ∂y∂z h], [∂z∂x h, ∂z∂y h, ∂z∂z h]]
alternative return (slower) but equivalent:
    # return 2*real.(diag(va * (interband_MM(b, ωs, vels, which_mm = which_mm) ./ (ωs .^2))))  
        + 0 .* 1/2 * real(contracted_sum_qm(a, b, ωs, vels, vvels))
"""
function F(a, b, ωs, vels, vvels, QM_switch; which_mm = :orbital)   
    va = vels[symb_to_ind(a)]
    va .-= diagm(diag(va))
    return 2*real.(vec(sum(va .* transpose(interband_MM(b, ωs, vels, which_mm = which_mm) ./ (ωs .^2)),
         dims=2))) + ifelse(QM_switch == true, 1, 0) .* 1/2 * 
        real(contracted_sum_qm(a, b, ωs, vels, vvels)) # vectorized version
end

"""
levi civita contracted. ωs are energy differences not to confuse with ϵs
"""
function contracted_sum_qm(a, b, ωs,vels, vvels; which_mm = :orbital)
    s = zeros(ComplexF64, size(ωs,1))
    cords = [:x,:y,:z]
    for c in cords
        for d in cords
            s += ϵ(b,c,d) .* qm_int(a, c, d, ωs,vels, vvels)  
        end
    end
    return real(s)
end

function qm_int(a, c, d, ωs, vels, vvels)
    va = vels[symb_to_ind(a)]
    vc = vels[symb_to_ind(c)]
    vd = vels[symb_to_ind(d)]
    Δd = diag(vd) .-  transpose(diag(vd))
    va .-= diagm(diag(va))
    vc .-= diagm(diag(vc))
    vd .-= diagm(diag(vd)) # because va[n,n] = 0 i can sum over all m
    return diag((Δd .* va ./ ωs .^ 3) * vc) + diag((vd ./ ωs .^ 2)  * vvels[symb_to_ind(a)][symb_to_ind(c)])
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

function interband_OMM(a, ωs, vels)
    cords = [:x,:y,:z]
    s = zeros(ComplexF64, size(ωs,1), size(ωs,1))
    for b in cords
        for c in cords
            s .+= ϵ(a, b, c) .* omm_int(vels[symb_to_ind(b)], vels[symb_to_ind(c)], ωs)
        end
    end
    return  s # the -i is the result of writting the off-diagonal terms of the connections as velocities
end

function omm_int(vb, vc, ωs)
    vb_diag = diag(vb) # store the diagonal elements
    vc_diag = diag(vc)
    vb .-= diagm(vb_diag) # clean the diagonal entries so I don't need to exclude from the sum
    vc .-= diagm(vc_diag)
    ωs .+= diagm(ones(size(vb,1))) # I add dummy entries to the diagonal to avoid 0/0
    M =  vb * (vc ./ ωs) +  (vb_diag .+ transpose(vb_diag)) .* (vc ./ ωs)
    return M .* (-1im/2)
end
#_________________________________________________________________________________________

function contracted_sum_qm(a::Symbol, b::Symbol, c::Symbol)
    s = 0.0
    for d in [:x,:y,:z]
        s += ϵ(b, c, d) * partial_quantum_metric(c, d, a)
    end
    return s
end

function ϵ(i::Symbol, j::Symbol, k::Symbol)
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

transform_k(us, Gs) = transform_k(us, Gs[1], Gs[2], Gs[3])
function transform_k(us,b1,b2,b3)
    return us[1]*b1 + us[2]*b2 + us[3]*b3
end
