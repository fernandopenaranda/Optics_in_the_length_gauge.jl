""" 
adaptive integration for 2D
    bz_integration_optical(f,xbounds, ybounds,  ωlist, evals) 
f integrand, xbounds, ybounds
it assumes a rectangular BZ. with frequency ω dependency
old---
xbounds and y bounds have to be adimensional
up to a bz_surface = (1/(2pi))^d with d = 2, 3.
which is included when calling bz_integration_functs...
new---
xbounds and ybounds are in [-0.5,0.5]. Assuming periodicity in the BZ
k = sum u_i * b_i, where b_i are the reciprocal lattice vector.
with this parametrization of k in terms of u_i, we must correct the integral by
multiplying with the Jacobian of the transformation given by `bz_volume`. 
Note that b_i may or may not be adimensional
"""
bz_integration_optical(f, p::Optical_computation_presets; kws...) = 
    bz_integration_optical(f, p.xbounds, p.ybounds, p.ωlist, p.evals; kws...) 

function bz_integration_optical(f, xbounds, ybounds, ωlist, evals; rel_tol = 1e-5, abs_tol = 0) 
    val, _ = Cubature.hcubature(length(ωlist), (x,v) -> v[:] = f(x), 
        [xbounds[1], ybounds[1]], [xbounds[2], ybounds[2]]; reltol = rel_tol, abstol= abs_tol, maxevals=Int(evals));
    return val
end
""" 
adaptive integration for 2D
    bz_integration_optical(f,xbounds, ybounds,  ωlist, evals) 
f integrand, xbounds, ybounds
it assumes a rectangular BZ.

xbounds and y bounds have to be adimensional
up to a bz_surface = (1/(2pi))^d with d = 2, 3
"""
bz_integration_transport(f, p::Transport_computation_presets; kws...) = 
    bz_integration_transport(f, p.xbounds, p.ybounds, p.evals; kws...) 

function bz_integration_transport(f, xbounds, ybounds, evals; rel_tol = 1e-5, abs_tol = 0) 
    val, _ = Cubature.hcubature(f, [xbounds[1], ybounds[1]], [xbounds[2], ybounds[2]];
        reltol = rel_tol, abstol= abs_tol, maxevals= Int(round(evals, digits = 0)))
    return val
end

function bz_integration_transport_3d(f, p::Transport_computation_3d_presets; kws...)
    if p.integration_method == :hcubature
        println("Adaptive integration")
        bz_integration_transport_3d_hcubature(f, p.xbounds, p.ybounds, p.evals; kws...) 
    elseif p.integration_method == :montecarlo
        println("Montecarlo integration")
        bz_integration_transport_3d_montecarlo(f, p.xbounds, p.ybounds, p.evals)
    elseif p.integration_method == :uniform_grid
        println("Uniform integration")
        bz_integration_transport_3d_uniformgrid(f, p.xbounds, p.ybounds, p.evals)
    else 0.0 end
end
function bz_integration_transport_3d_hcubature(f, xbounds, ybounds, evals; rel_tol=1e-5, abs_tol=0)
    val, _ = Cubature.hcubature(
        f,
        xbounds,
        ybounds;
        reltol = rel_tol,
        abstol = abs_tol,
        maxevals = Int(round(evals))
    )
    return val
end

function bz_integration_transport_3d_montecarlo(f, xbounds, ybounds, evals)
    dim = length(xbounds)
    pts = QuasiMonteCarlo.sample(evals, dim, SobolSample())
    acc = 0.0
    for i in 1:evals
        u = pts[:,i]
        x = xbounds + u .* (ybounds - xbounds)
        acc +=  f(x)
    end
    return cube_volume(xbounds, ybounds) * acc / evals
end

function bz_integration_transport_3d_uniformgrid(f, xbounds, ybounds, evals)
    N = floor(Int, evals^(1/3))
    it = 0
    if N == 0
        throw(ArgumentError("N = 0, increase number of evals"))
    end
    xs = range(xbounds[1], ybounds[1], length=N)
    ys = range(xbounds[2], ybounds[2], length=N)
    zs = range(xbounds[3], ybounds[3], length=N)
    acc = 0.0
    for x in xs, y in ys, z in zs
        it += 1
        acc += f(SVector(x,y,z))
    end
    return cube_volume(xbounds, ybounds) * acc / it
end



bz_volume(b1,b2,b3) = abs(dot(b1, cross(b2,b3)))

""" volume of the cube in parameter u_i space. By default it is 1 for [-0.5, 0.5] integration domains"""
cube_volume(xbounds, ybounds) = prod(ybounds .- xbounds)