""" 
adaptive integration for 2D
    bz_integration_optical(f,xbounds, ybounds,  ωlist, evals) 
f integrand, xbounds, ybounds
it assumes a rectangular BZ. with frequency ω dependency

xbounds and y bounds have to be adimensional
up to a bz_surface = (1/(2pi))^d with d = 2, 3
"""
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
function bz_integration_transport(f, xbounds, ybounds, evals; rel_tol = 1e-5, abs_tol = 0) 
    val, _ = Cubature.hcubature(f, [xbounds[1], ybounds[1]], [xbounds[2], ybounds[2]];
        reltol = rel_tol, abstol= abs_tol, maxevals=Int(evals));
    return val
end