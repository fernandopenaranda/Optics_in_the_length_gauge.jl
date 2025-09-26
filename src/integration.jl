
""" 
adaptive integration for 2D
    bz_integration(f,xbounds, ybounds,  ωlist, evals) 
f integrand, xbounds, ybounds
"""
function bz_integration(f, xbounds, ybounds, ωlist, evals) 
    val, _ = Cubature.hcubature(length(ωlist), (x,v) -> v[:] = f(x), 
        [xbounds[1], ybounds[1]], [xbounds[2], ybounds[2]]; reltol = 1e-5, abstol=0, maxevals=evals);
    bz_surface = (1/(2pi*a0))^2 
    return bz_surface .* val
end
