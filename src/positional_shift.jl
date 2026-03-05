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
"""
quantum_contribution(p::Quantum_correction_σijk_antisym) = 
quantum_contribution(p.a0, p.dirJ, p.dirE, p.dirB, p.h, p.dh, p.ddh, p.τ, p.T, p.cpt, 
p.q, p.which_mm)

function quantum_contribution(a0, dirJ, dirE, dirB, h, dh, ddh, τ, T, cpt, q, which_mm)
    checkdims(q)
    checkantisym(a,b,c)
    integrand(q) = integrand_quantum_contribution(dirJ, dirE, dirB, h, dh, ddh, T, q, which_mm = which_mm)
    integrator(observable) = bz_integration_transport(observable, cpt, rel_tol = rel_tol, abs_tol = abs_tol)
    bz_vol = 1/(2pi*a0*ang_to_m)^length(q)
    val = bz_vol * integrator(integrand)
    return val
end
"there is also a contribution coming from the spin magnetic moment (path disabled for the moment)
f' * ∑_α (v_a F_{bc}^α - v_b F_{ac}^α + ϵ_{abd}Omega_d M_c^α).
Important not to confuse"
function integrand_quantum_contribution(a, b, c, h, dh, ddh, T, q; which_mm = :orbital)
    ϵs, ψs = eigen(Matrix(h(q)))   
    dhs = [dh(q)[1],dh(q)[2],dh(q)[3]]
    ddhs = [[ddh(q)[1][1],ddh(q)[1][2],ddh(q)[1][3]], 
            [ddh(q)[2][1],ddh(q)[2][2],ddh(q)[2][3]],
            [ddh(q)[3][1],ddh(q)[3][2],ddh(q)[3][3]]]
    integrand_quantum_contribution(a, b, c, ϵs, ψs, dhs, ddhs, T, which_mm = which_mm)
end

function integrand_quantum_contribution(a, b, c, ϵs, ψs, dh, ddh, T; which_mm = :orbital)
    ωs = Ω(ϵs)
    vels = [v(:x,ψs,dh), v(:y,ψs,dh), v(:z,ψs,dh)] # factor needed
    vvels = [[dv(:x, :x, ψs, ddh), dv(:x, :y, ψs, ddh), dv(:x, :z, ψs, ddh)],
             [dv(:y, :x, ψs, ddh), dv(:y, :y, ψs, ddh), dv(:y, :z, ψs, ddh)],
             [dv(:z, :x, ψs, ddh), dv(:z, :y, ψs, ddh), dv(:z, :z, ψs, ddh)]]   #* ang_to_m^2/ ħ_ev_s
    return real(sum(d_f(ϵs, 0, T) .* (positional_shift(a, b, c, ωs, vels, vvels, which_mm = which_mm) .+ berry_OMM(a,b,c, ωs, vels, which_mm = which_mm))))
end
"""contribution to the quantum correction coming from the product of the Berry curvature and 
the magnetic moment (orbital + spin). Spin deactivated."""
function berry_OMM(a,b,c, ωs, vels; which_mm = :orbital)
    s = 0.0im
    for d in [:x,:y,:z]
        s += ϵ(a,b,d)*Berry_curvature(d)*interband_MM(c, ωs, vels, which_mm = which_mm)
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
    s = 0.0im
    coords = [:x,:y,:z]
    ωs .+= diagm(ones(size(ωs,1))) 

    for j in coords
        for k in coords
            s .+= ϵ(i,j,k) .* Ωi_aux(j,k, ωs ,vels)
        end
    end
    return 1im .* s #hbar^2
end
function  Ωi_aux(j,k, ωs ,vels)
    vj = vels[symb_to_ind(j)]
    vk = vels[symb_to_ind(k)]
    return diag((vj .- diagm(diag(vj))./(ωs.^2)) *(vk .- diagm(diag(vk)))) 
end

""" contribution to the quantum correction to σ_ijk^(A,1) coming from the positional 
shift contribution F """
function positional_shift(a, b, c, ωs, vels, vvels; which_mm = :orbital) 
    va = vels[symb_to_ind(a)]
    vb = vels[symb_to_ind(b)]
    return diag(va) .* F(b,c, ωs, vels, vvels, which_mm = which_mm) .-
             diag(vb) .* F(a, c, ωs, vels, vvels, which_mm = which_mm) 
end

"""
anomalous spin and orbital polarizability. Default only P see the definition of MM(which=:orbital)
dh = [∂x h, ∂y h, ∂z h]
ddh =  [[∂x∂x h, ∂x∂y h, ∂x∂z h], [∂y∂x h, ∂y∂y h, ∂y∂z h], [∂z∂x h, ∂z∂y h, ∂z∂z h]]
"""
function F(a, b, ωs, vels, vvels; which_mm = :orbital)   
    va = vels[symb_to_ind(a)]
    va .-= diagm(diag(va))
    2*imag.(diag(va * (interband_MM(b, ωs, vels, which_mm = which_mm) ./ (ϵs .^2)))) #using re(-i*x) == im(x)
    + 1/2 * real(contracted_sum_qm(a, b, ωs, vels, vvels))
end

"""
levi civita contracted. ωs are energy differences not to confuse with ϵs
"""
function contracted_sum_qm(a, b, ωs,vels, vvels; which_mm = :orbital)
    s = 0.0im
    cords = [:x,:y,:z]
    for c in cords
        for d in cords
            s += epsilon(b,c,d) .* qm_int(a, c, d, ωs,vels, vvels)  
        end
    end
    return real(sum)
end

function qm_int(a, c, d, ωs, vels, vvels)
    va = vels[symb_to_ind(a)]
    vc = vels[symb_to_ind(c)]
    vd = vels[symb_to_ind(d)]
    Δd = diag(vd) .-  diag(vd)'
    ωs .+= diagm(ones(size(va,1))) # avoid 0/0
    va .-= diagm(diag(va))
    vc .-= diagm(diag(vc))
    vd .-= diagm(diag(vd)) # because va[n,n] = 0 i can sum over all m
    return diag((Δd .* va ./ ωs .^ 3) * vc) + diag((vd ./ ωs .^ 2)  * vvels[a][c])

end

""" interband magnetic moment with orbital and spin parts """
function interband_MM(a,ωs, vels; which_mm = :orbital)
    if which_mm == orbital
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
    for b in cords
        for c in cords
            s += ϵ(a, b, c) .* omm_int(vels[symb_to_ind(b)], vels[symb_to_ind(c)], ωs)
        end
    end
    return  -1im * s/2 # the -i is the result of writting the off-diagonal terms of the connections as velocities
end

function omm_int(vb, vc, ωs)
    vb_diag = diag(vb) # store the diagonal elements
    vc_diag = diag(vb)
    vb .-= diagm(vb_diag) # clean the diagonal entries so I don't need to exclude from the sum
    vc .-= diagm(vc_diag)
    ωs .+= diagm(ones(size(vb,1))) # I add dummy entries to the diagonal to avoid 0/0
    M = vb * (vc ./ ωs) + (vb_diag .+ vb'_diag) .* (vc ./ ωs)
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
v(a::Int,ψs, dh) = vel(ψs, dh[a]) #* ang_to_m/ ħ_ev_s

dv(a::Symbol, b::Symbol, ψs, ddh) = dv(symb_to_ind(a), symb_to_ind(b), ψs, ddh)
dv(a::Int,b::Int, ψs, ddh) = vel(ψs, ddh[a][b]) #* ang_to_m/ ħ_ev_s

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