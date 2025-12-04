"""
    Hamiltonian in eV
    Units: [σxx]/[B]
    returns the symmetric σxxx linear magneto conductivity conductivity
    `linear_magneto_conductivity(params::planar_σijk_presets)``
    Input params see structs.jl:
    `dirJ::Symbol`, `dirE::Symbol`, directions of J and E (tensor component of the conductivity)`::Symbols -> {:x, :y}` unbounded 2D directions
    `h::Function` is a Hamiltonian function with dependence on momentum `q`
    `dh::Vector{Function}` is ∇h along the x and y directions (2D). An array of functions `dh = [dh_x(q), dh_y(q)]` of momentum `q`
    `τ::Float64` scattering time
    `T::Float64` temperature
    `evals::Int` number of steps in the adaptive integration
    `rz` rz operator(q)
    details on the Hamiltonian (the chemical potential and dielectric field...) are set in the parameters of the Hamiltonian stored in `planar_σijk_presets`
"""
linear_magneto_conductivity(params::Planar_σijk_presets) =
    linear_magneto_conductivity(
        params.dirJ, params.dirE, params.dirB, params.h, params.nabla_h, 
        params.nabla_nabla_h, params.rz, params.τ, params.T,
        params.berry_contribution, params.omm_contribution, params.fermi_surface, params.with_shift,
        params.computation.xbounds, params.computation.ybounds, params.computation.evals)

function linear_magneto_conductivity(i,j,k, h, dh, ddh, rz, τ, T, Ω_contr, omm_contr, #N
    fermi_surface, with_shift, xbounds, ybounds, evals; rel_tol = 1e-5, abs_tol = 0)
    integrand(q) = k_linear_magneto_conductivity(i, j, k, h, dh, ddh, rz, q; T = T, τ = τ, 
        Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface, with_shift = with_shift) 
    val = bz_integration_transport(integrand, xbounds, ybounds, evals, rel_tol = rel_tol, abs_tol = abs_tol)
    bz_vol = (1/(2pi*ang_to_m))^(length(xbounds))                                                              
    return bz_vol * val
end

function k_linear_magneto_conductivity(i::Symbol, j::Symbol, k::Symbol, h, dh, ddhi, rz::Function, q; 
    T = 2, τ = 1e-15, Ω_contr = true, omm_contr = true, fermi_surface = false, with_shift = true)
    ϵs, ψs = eigen(Matrix(h(q)))                                                                          
    C = 2π * τ                                       
    σxxx = C * k_linear_mr_integrand(i, j, k, ϵs, ψs, rz(q, ψs), dh(q)[1], dh(q)[2], ddhi(q), 0, T,           
        Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface)                           # generalize to σyyy too   
    if with_shift == false                                                               
        return σxxx
    else
        σxxx_shift = C * k_linear_mr_integrand_shift(i, j, k, ϵs, ψs, rz(q, ψs), dh(q)[1], dh(q)[2], ddhi(q), μ, T)
        return σxxx + σxxx_shift  
    end                      
end

"the term with vij is only valid for sigma xxx"
function k_linear_mr_integrand(i, j, k, ϵs, ψs, rzmat, dhx, dhy, dhxx, μ, T;
         Ω_contr = true, omm_contr = true, fermi_surface = false)                       # units meters, meV, seconds
    omega = Ω(ϵs)
    Δx = Δ(ψs, dhx) * ang_to_m
    Δy = Δ(ψs, dhy) * ang_to_m
    rx = r(ϵs, ψs, dhx) * ang_to_m
    ry = r(ϵs, ψs, dhy) * ang_to_m
    vx = vel(ψs, dhx) * ang_to_m/ ħ_ev_s
    vy = vel(ψs, dhy) * ang_to_m/ ħ_ev_s
    vxx = vel(ψs, dhxx) * ang_to_m^2/ ħ_ev_s
    rzmat *= ang_to_m
    Ω_switch = ifelse(Ω_contr == true, 1, 0)
    omm_switch = ifelse(omm_contr == true, 1, 0)
    if fermi_surface == true
        return sum(d_f(ϵs, 0, T))
    else
        return real(sum(d_f(ϵs, 0, T) .*
            (omm_switch .* mr_omm(i, j, omega, rx, ry, vx, vy, Δx, Δy, rzmat) + 
            Ω_switch .* mr_Ω(i, j, k, rzmat, rx, ry, vx, vy) + 
            - omm_switch .* mr_vij(i, vy, rzmat, vxx))))                                              #only valid in the xx direction
    end
end

"""shift correction due to the magnetic field effect on the bandstructure"""
function k_linear_mr_integrand_shift(i, j, k, ϵs, ψs, rzmat, dhx, dhy, dhxx, μ, T)
    ry = r(ϵs, ψs, dhy) * ang_to_m
    rzmat *= ang_to_m
    vxx = vel(ψs, dhxx) * ang_to_m^2/ ħ_ev_s
    vy = vel(ψs, dhy) * ang_to_m/ ħ_ev_s
    if fermi_surface == true
        return sum(d_f(ϵs, 0, T))
    else
        δμ_shift(i, ϵs, T, vy, ry, rzmat) * vij_shift(ϵs, T, vxx)/sum(d_f(ϵs, 0, T))
    end
end

""""
correction due to switching into the canonical ensemble
"""
vij_shift(ϵs, T, vij) = sum(d_f(ϵs, 0, T) .* real(diag(vij))) #check this
mr_vij(i, vj,rz, vij) = real(OMM(i, vj, rz) .* diag(vij))

δμ_shift(i, ϵs, T, vj, rj, rz) = sum(d_f(ϵs, 0, T) .* (OMM(i, vj, rz)) + 
   (Ωin(i, rj, rz)) .* f(ϵs, 0, T)/ħ_ev_s) #units 1/e m^2

# ----------------------------------------------------------------------------------------
#           Magnetorresistance integrand terms: m_Ω, m_OMM, and n-fixing correction
#                          before convoluting with Fermi derivatives                              
# ----------------------------------------------------------------------------------------
"""
Orbital magnetic moment contribution to the 
linear magnetorresistance (planar case)                          UNITS:  e
# .-vijmn * OMM() this is 0 read comment above.
"""
function mr_omm(i, j, omega, rx, ry, vx, vy, Δx, Δy, rz) 
    vi = which_mat(i, vx, vy)
    vj = which_mat(j, vx, vy)
    return 0.5 .* (real(diag(vi)) .* d_OMM(j, i, omega, rx, ry, Δx, Δy, rz) .+ 
    real(diag(vj)) .* d_OMM(i, j, omega, rx, ry, Δx, Δy, rz)) 
end
""" 
Berry curvature contribution to the 
linear magnetorresistance (planar case)                         UNITS:  e
note: k ∈ {x,y}
"""
function mr_Ω(i, j, k, rz, rx, ry, vx, vy)
    err_knotin(k)
    vi = which_mat(i, vx, vy)
    vj = which_mat(j, vx, vy)
    r_not_k = which_mat(k, ry, rx)
    # note that above ry and rx are interchanged the reason is that 
    #if k == :x then we need notk == :y for Ωin
    return real((diag(vi) .* diag(vj) .* Ωin(k, r_not_k, rz) .- 
            (δ_kron(j,k) .* diag(vi) .+ δ_kron(i,k) .* diag(vj)) .* 
            (diag(vx) .* Ωin(:x, ry, rz) .+ diag(vy) .* Ωin(:y, rx, rz)))) # this is real diag(vx) is in R
end

function mr_corr() #not sure if it vanishes due to vij. Lets not include it for now
    ε()
end
# ----------------------------------------------------------------------------------------
#              Berry curvature, orbital magnetic moment and its derivative
# ----------------------------------------------------------------------------------------
"""
    Ωin(ϵs, ψs, dh, k)
computes the in-plane Berry curvature
Ω^α = 2ħ ϵ^{ij} Re[Σ_n'≠n -iv_{nn'}^j * z_{n'n} / ω_{nn'}]      
"""
Ωin(i, rj, rz) = 2ε(i) * imag(Σ_nondiag(rj, rz))
"""
    Planar orbital magnetic momentum: OMM                          UNITS e
    one hbar comes from v and the other from the definition of OMM
"""
OMM(i, vj, rz) = 2ε(i) * real(Σ_nondiag(vj,rz))
"""
    k-derivative of the planar, orbital magnetic moment.            UNITS e
"""
function d_OMM(j, k, omega, rx, ry, Δx, Δy, rz) 
    Δj = which_mat(j, Δx, Δy)
    Δk = which_mat(k, Δy, Δx)
    rj = which_mat(j, rx, ry)
    rk = which_mat(k, ry, rx) # note that here rx and ry have to be reversed
    return -ε(j) * 
        imag(
            Σ_nondiag(Δj .* rk, rz) +   
            Σ_nondiag(omega .* rk, z_covariant(rz, rj)) + 
            Σ_nondiag(omega .* r_covariant(omega, rk, Δk, rj, Δj), rz))/ħ_ev_s
end # units L^3/s

# ----------------------------------------------------------------------------------------
#                               Auxiliary functions
# ----------------------------------------------------------------------------------------
""" 
    off diagonal product of two matrices returns a vector
    Let A and B be two matrices. It computes 
    [A_{nq}B_{qn} for n in nlist] with n neq q.
"""
Σ_nondiag(A, B) = diag(A * nondiag(B))
nondiag(A) = A .- diagm(diag(A))

function offdiag_sum(A, B)
    @assert size(A) == size(B) "A and B must be the same size"
    s = 0.0
    for n in axes(A, 1), m in axes(A, 2)
        if n != m
            s += A[n, m] * B[m, n]
        end
    end
    return s
end

""" Levi-civita tensor """
function ε(a::Symbol)
    if a==:x
        return 1 
    elseif a==:y
        return -1 
    end
end

function ε(a::Symbol, b::Symbol)
    if a == b
        return 0  
    elseif a==:x && b == :y
        return 1  
    elseif a==:y && b == :z
        return 1  
    elseif a==:y && b == :x
        return -1 
    elseif a==:z && b == :y
        return -1 
    end
end

function ε(a::Symbol, b::Symbol, c::Symbol)
    if a == b || a==c || b==c
        return 0  
    elseif a==:x && b == :y && c == :z
        return 1  
    elseif b==:x && c == :y && a == :z
        return 1  
    elseif c==:x && a == :y && b == :z
        return 1  
    elseif a==:z && b == :y && c == :x
        return -1  
    elseif a==:y && b == :x && c == :z
        return -1
    elseif a==:x && b == :z && c == :y
        return -1
    end
end

""" kronecker delta tensor """
function δ_kron(a::Symbol, b::Symbol)
    if a == b
        return 1
    else return 0
    end
end

""" matrix selector with direction """
which_mat(i, matx, maty) = ifelse(i == :x, matx, maty)

function err_knotin(k::Symbol)
    if k == :z
        throw(ArgumentError(
            "Error: k = :z is not allowed; k ∈ {x, y}"))
    end
end

function err_dirs(i,j,k)
    if in(:z, [i,j,k])
        throw(ArgumentError(
            ":z is not allowed; i, j, k ∈ {x, y}"))
    end
end