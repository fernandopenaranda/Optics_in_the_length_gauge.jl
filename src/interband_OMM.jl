"""
computes the interband_orbital magnetic moment. 
IT IS WRITTEN IN [Ev * Å^2] without the e/ħ prefactor
q has length units of Å^-1.
Valid for 3d any direction 
Valid for 2d only for the bounded dir that is fixed to be z.
In plane OMM follows a different path.
"""
function interband_OMM(p::OMM_presets, q)
    ϵs, ψs = eigen(Matrix(p.h(q)))
    dims = length(q)
    if dims == 3
        dhs = [p.nabla_h(q)[1],p.nabla_h(q)[2],p.nabla_h(q)[3]]
        vels = [v(:x,ψs,dhs), v(:y,ψs,dhs), v(:z,ψs,dhs)]  #units [E*L]
    else 
        if p.dir != :z
            throw(ArgumentError("Bounded direction in 2d must be z"))
        end
        dhs = [p.nabla_h(q)[1],p.nabla_h(q)[2]]
        vels = [v(:x,ψs,dhs), v(:y,ψs,dhs)]  #units [E*L]
    end
    ωs = Ω(ϵs) 
    return  interband_OMM(p.dir, ωs, vels)
end


function interband_OMM(a, ωs, vels)
    s = zeros(ComplexF64, size(ωs,1), size(ωs,1))
    if length(vels) == 3
        cords = [:x,:y,:z]
    else  
        cords = [:x,:y] # it assumes the bounded direction is z and a and b are in the plane
    end
    for b in cords
        for c in cords
            s .+= levi_civita(a, b, c) .* omm_int(vels[symb_to_ind(b)], vels[symb_to_ind(c)], ωs)
        end
    end
    return  s
end

function omm_int(vb, vc, ωs)
    ωs_safe = ωs #.+ Diagonal(fill(Inf, size(ωs,1)))
    vb_diag = diag(vb) # store the diagonal elements
    non_diagvb = copy(vb) - diagm(diag(vb))
    non_diagvc = copy(vc)- diagm(diag(vc))
    M = non_diagvb * (non_diagvc  ./ ωs_safe ) +  (vb_diag .+ vb_diag') .* (non_diagvc  ./ ωs_safe )
    return  M .* (-1im/2)
end


