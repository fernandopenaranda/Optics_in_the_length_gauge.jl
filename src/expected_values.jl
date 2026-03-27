"""
    `filling(h, μ, Δx, T; kws...)`
returns the band filling over the integration area set by Δx
    `filling(h, μ, Δx, Δy, T; kws...)`
returns the band filling over the integration area set by Δx, Δy
note that if the integration area is not the full BZ, as it happens
usually when dealing with continuum models, a postprocessing regularization
by a momentum cutoff will be required. Since this postprocessing is model
dependent no regularization path is included within `Optics_in_the_length_gauge`.
"""
filling(h, μ, Δx, T; kws...) = expected_value(h, μ, 1I, Δx, T; kws...)
filling(h, μ, Δx, Δy, T; kws...) = expected_value(h, μ, 1I, Δx, Δy, T; kws...)

"""
    `expected_value(h, μ, op, Δx, T; kws...)`
returns the expected value of a given operator (same dims as h)
over the integration area set by Δx
    ` expected_value(h, μ, op, Δx, Δy, T; evals = 100)`
returns the expected value of a given operator (same dims as h)
    over the integration area set by Δx , Δy
note that if the integration area is not the full BZ, as it happens
usually when dealing with continuum models, a postprocessing regularization
by a momentum cutoff will be required. Since this postprocessing is model
dependent no regularization path is included within `Optics_in_the_length_gauge`.
"""
expected_value(h, μ, op, Δx, T; kws...) = expected_value(h, μ, op, Δx, Δx, T; kws...)

function expected_value(h, μ, op, Δx, Δy, T; evals = 100, rel_tol= 1e-5, abs_tol = 0)
    integrand(q) = k_expected_value(h, μ, op, q, T)
    val = bz_integration_transport(integrand, Δx, Δy, evals; rel_tol = rel_tol, abs_tol = abs_tol) 
    return val /((Δx[2]-Δx[1])*(Δy[2]-Δy[1]))
end 

function k_expected_value(h, μ, op, q, T)
    check_same_size(h(q), op)
    ϵs, ψs = eigen(h(q))
    if  T == 0
        inds  = findall(ϵ-> ϵ < μ, real(ϵs))
        s = 0.0 
        for ind in inds
            s += ψs[:,ind]'* op * ψs[:,ind]
        end
    else 
        s = 0.0 
        for (i, ϵ) in enumerate(ϵs)
            s += ψs[:,i]'* op * ψs[:,i] * fn(real(ϵ), μ, T)
        end
    end

    return real(s)
end

check_same_size(A::AbstractMatrix, B::UniformScaling) = nothing
function check_same_size(A::AbstractMatrix, B::AbstractMatrix)
    if size(A) != size(B)
        error("Matrix dimensions do not match: size(A) = $(size(A)), size(B) = $(size(B))")
    end
end

"""
creates a k mesh with the symmetries of the problem, usefull for visualization in k-space
evaluation of the observable along the ab plane with ab ≠ c ∈ {x,y,z} at u_c = u0 ∈ [-0.5,0.5]
a, b are the plane coordinates quantity is the function in Optics_in_the_length_gauge. 
and p the presets compatible with such function
-CairoMakie function for visualization.
    function plot_kresolved(kx, ky, Z, u0; label = "f(x, y)")
        fig = Figure(size = (600, 500))
        ax = Axis(fig[1, 1], xlabel="kx [1/Å]", ylabel="ky [1/Å]", title = "uz =  [adimensional]")
        cmap = :grays  # CairoMakie built-in colormap, light-to-dark blues
        hm = surface!(ax, kx, ky, zeros(size(Z)),
        color=Z,
        shading=NoShading,
        colormap=cmap,
        colorrange=(-maximum(abs,Z), maximum(abs,Z)))
    Colorbar(fig[1, 2], hm, label=label)
    ylims!(ax, -2π-π/2, 2π+π/2)
    end


"""
function k_mesh_eval(quantity, p, R1, R2, R3, c_symb::Symbol; u0 = 0,  botbounds = [-0.5,-0.5], topbounds = [0.5,0.5], kpoints = 100)
    Gs = dualbasis([R1,R2,R3])
    N = floor(Int, kpoints^(1/3))
    uas = range(botbounds[1], topbounds[1], length=N)
    ubs = range(botbounds[2], topbounds[2], length=N)
    c = Optics_in_the_length_gauge.symb_to_ind(c_symb)
    u_param(ua,ub, c) = ifelse(c == 1, [u0,ua,ub], ifelse(c == 2, [ua,u0,ub], [ua,ub,u0]))
    aux_f(ua,ub,c) = Optics_in_the_length_gauge.transform_k(SVector(ua,ub,c), Gs)
    a, b = setdiff(1:3, [c])
    ka = [transform_k(u_param(ua,ub,c), Gs)[a] for ua in uas, ub in ubs]
    kb = [transform_k(u_param(ua,ub,c), Gs)[b] for ua in uas, ub in ubs]
    f(ua,ub) = quantity(p, transform_k(u_param(ua,ub,c),Gs))
    Z = [f(ua,ub) for ua in uas, ub in ubs]
    return (a,b), ka, kb, Z
end

function  k_mesh_eval(p, R1, R2;  botbounds = [-0.5,-0.5], topbounds = [0.5,0.5], kpoints = 100)
    Gs = dualbasis([R1,R2,R3])
    N = floor(Int, kpoints^(1/2))
    xs = range(botbounds[1], topbounds[1], length=N)
    ys = range(botbounds[2], topbounds[2], length=N)
    kx = [transform_k(SVector(x,y), Gs)[1] for y in ys, x in xs]
    ky = [transform_k(SVector(x,y), Gs)[2] for y in ys, x in xs]
    f(x,y) = integrand_quantum_contribution_q(comp_pres, [x,y,zval])
    Z = [quantity(p, [kx[i,j], ky[i,j]]) for i in 1:N, j in 1:N]
    return kx, ky, Z
end