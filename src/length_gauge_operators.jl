#_________________________________________________________________________________________
# position operators
#_________________________________________________________________________________________

r(ϵs, ψs, dh) = -1im .* vel(ψs, dh) ./ Ω(ϵs)                                    # Units: Å  
# rz(ψs, dim) = ψs' * (zop(dim) .* ψs)                                          # Units: Å
function zop(dim)                                                                           
    dim2 = div(dim, 2)  
    return  3.3/2 .* vcat(ones(dim2), -1 .* ones(dim2))       # 3.3 Å = interlayer distance                                            
end

function r_covariant(omega, ra, Δa, rb, Δb)                                       # Units: Å^2
    raa_cov = zeros(ComplexF64, size(ra,1), size(ra,1))#similar(ra)
    dim = size(omega, 1) 
    for n in 1:dim
        for m in 1:dim
            if m != n && n > m 
                aux = 0
                for p in 1:dim
                    if p != n && p != m
                        aux += 1im * (omega[n,p]*ra[n,p]*rb[p,m] - 
                                        omega[p,m]*rb[n,p]*ra[p,m])            # Units: eV * Å * Å 
                    else
                        aux += 0
                    end
                end
                raa_cov[n,m] = - ((ra[n,m] * Δb[n,m] + rb[n,m] * Δa[n,m])
                               + aux)/omega[n,m]                            # Units: (Å * Å * eV) / eV
            else 
                nothing
            end
        end
    end
    return raa_cov
end

function z_covariant(zmat, rmat)
    dim = size(zmat, 1)
    mat = zeros(ComplexF64, size(rmat,1),size(rmat,1))
    @inbounds begin
        for n in 1:dim
            for m in 1:dim
                if m != n && n > m 
                    for p in 1:dim
                        if  p == n 
                            val = zmat[n,p] * rmat[p,m] 
                        elseif p == m 
                            val = -rmat[n,p] * zmat[p, m]
                        else
                            val = zmat[n,p] * rmat[p,m] - rmat[n,p] * zmat[p, m]
                        end
                        mat[n,m] += val    
                    end
                else nothing end
            end
        end
    end
    return -1im .* mat
end

#_________________________________________________________________________________________
# observables
#_________________________________________________________________________________________

lorentz(ϵs, ω, η) = 1/π .* η ./ (δ(ϵs, ω).^2 .+ η^2) 
δ(ϵs, ω) = ω .- Ω(ϵs)  # the minus sign comes from delta(e_f-e_i-omega); ωnm = -ωmn
# this definition is weak against degeneracies but it is highly efficient (benchmark)
Ω(ϵs) = ϵs .- t(ϵs) +  1e-7 * Diagonal(ϵs) .+ 1e-35# the last term is just to avoid infinities 
#it'll cancel out later
# Δ(ψs, dh) = vel(ψs, dh) - t(vel(ψs, dh))
vel(ψs, dh) = ψs' * dh * ψs # It really is hbar v  # Units: [H] Å  

function omega_r_Δ(ϵs, ψs, dh)
    vmat = vel(ψs, dh)
    mat = similar(vmat)
    omega = Ω(ϵs)
    for i in 1:size(vmat, 1)
        for j in 1:size(vmat, 1)
            mat[i, j] = vmat[i, i] - vmat[j, j]
        end
    end
    return omega, -1im .* vmat ./ omega, mat 
end

function Δ(ψs, dh)
    vmat = vel(ψs, dh)
    mat = similar(vmat)
    for i in 1:size(vmat, 1)
        for j in 1:size(vmat, 1)
            mat[i, j] = vmat[i, i] - vmat[j, j]
        end
    end
    return mat
end

"""
fn - fm where f are the Fermi Dirac distribution in the limit of T = 0.
I set p.
"""
function f!(mat, ϵs, μ, T, tol = 1e10) 
    replace!(mat, NaN => tol, NaN + NaN*im => tol *(1+im), 1im * NaN => tol*1im )
    mat .*= f(ϵs, μ, T)
end

#_________________________________________________________________________________________
# FERMI FUNCTION AND ITS DERIVATIVES
#_________________________________________________________________________________________

f(ϵs, μ, T) = [fn(ϵs[i], μ, T) - fn(ϵs[j], μ, T) for i in 1:length(ϵs), j in 1:length(ϵs)] #nothe this is a matrix difference of fermi functions
fn(ϵn,μ::Float64) = ifelse(ϵn < μ, 1.0, 0.0)
function fn(ϵn, μ, T)
    if T == 0
        return ifelse(ϵn < μ, 1.0, 0.0)
    else
        # return 1/(exp((ϵn - μ)/(kB * T)) + 1) # For RHG
        return 1/(exp((ϵn - μ)/(1e3*kB * T)) + 1)

    end
end
# first derivative with respect to E
d_f(ϵs, μ, T) = [d_fn(ϵn, μ, T) for ϵn in ϵs]


function d_fn(ϵn, μ, T)
        -1/(1e3*kB*T) * fn(ϵn, μ, T) * (1-fn(ϵn, μ, T)) # TBG
end 

# function d_fn(ϵn, μ, T) # this is for RHG
#     # if T === 0
#     η = π * kB*T
#     return (-1/π * η) / ((μ - ϵn) ^2 + η^2) # lorentzian fastest convergency
#     # else T ≠ 0
#         # -1/(kB*T) * fn(ϵn, μ, T) * (1-fn(ϵn, μ, T)) # not used at the moment sech method
#     # end
# end 
# second derivative with respect to E
d_d_f(ϵs, μ, T) = [d_d_fn(ϵn, μ, T) for ϵn in ϵs]
d_d_fn(ϵn, μ, T) = - 1/(kB*T)^2 * d_d_fn((ϵn-μ)/(kB*T))
d_d_fn(x) = exp(x) * (exp(x)-1)/(exp(x)+1)^3
#-------------------------------------------------------------------------------------------
# Aux operations
#-------------------------------------------------------------------------------------------
function whichBPGE(part, r_cov, r, r_cov2, r2)
    if part == :REAL
        return symmetrized_imagsum(r_cov, r, r_cov2, r2) 
    else
        return antisymmetrized_realsum(r_cov, r, r_cov2, r2)
    end
end

function whichBPGE(part, r_cov, r)
    if part == :REAL
        return symmetrized_imagsum(r_cov, r)
    else 
        return antysymmetrized_realsum(r_cov, r)
    end
end

symmetrized_imagsum(mat1, mat2)  =  imag(mat1 .* t(mat2) .+ mat2 .* t(mat1))
symmetrized_imagsum(mat1, mat2, mat3, mat4) = imag((mat1 .* t(mat2)) .+ (mat3 .* t(mat4))) 

function antisymmetrized_realsum(mat1, mat2) 
    auxmat = similar(mat1)
    dim = size(mat1, 1)
    for n in 1:dim
        for m in 1:dim
            auxmat[n,m] = real(mat1[n,m] * mat2[m,n] - mat2[n,m] * mat1[m,n])
        end
    end
    return auxmat
end

antisymmetrized_realsum(mat1, mat2, mat3, mat4) = real(mat1 .* t(mat2) .- mat3 .* t(mat4))

t(mat) = transpose(mat)
sum_nondiag(mat) = sum(mat) - sum(Diagonal(mat))
sum_lowerdiag(mat::Array) =  sum_lowerdiag(mat, size(mat,1))
sum_lowerdiag(mat, dim) = sum(mat .* (tril(ones(dim, dim), 0) - Diagonal(ones(dim,dim))))

function return_kwargs(kws, name::Symbol)
    kwargs = NamedTuple(kws)
    # println(kwargs)
    if haskey(kwargs, name)
        return kwargs[name]
    else  nothing
    end
end


#_________________________________________________________________________________________
# Extra
#_________________________________________________________________________________________
"""warn functions"""
function warn_equalargs(a, b)
    if a == b
        @warn "Arguments are equal; stopping."
        error("Arguments must not be equal.")
    end
end