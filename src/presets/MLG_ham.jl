#_________________________________________________________________________________________
# 
#                             MLG HAMILTONIAN
#
#_________________________________________________________________________________________

function MLG_hamiltonian(μ, ν, q)
    mat = spzeros(ComplexF64, 2, 2) 
    st = _sublattice_coupling(q, ν)
    mat[1:2,1:2] = [-μ st; conj(st) -μ]
    return mat
end

"intralayer coupling of graphene proyected into the layer and sublat
space:  `t * 1im * Σ_i δ^l_i'* (q + n1*G1 + n2*G2 - κ^l) exp(1im*δ_i^l*K_l)`. 
`l`` refers to the layer index and `δ^l_i` are the real space vectors connecting first 
neighbours in the rotated layer"
function _sublattice_coupling(q, ν)
    t = -2.46575
    b1 = 2π * [1, 1/√3]
    b2 = 2π * [-1, 1/√3]
    if ν == 1
        K = κ(b1, b2)
    else
        K = κ(b2, b1)
    end
    δ1 = [0, 1/√3]
    δ2 = [-1/2, -1/2/√3]
    δ3 =  [1/2, -1/2/√3]
    δs = [δ1, δ2, δ3]
    qvecs =  q -  K 
    return t * 1im * sum([δs[j]' * qvecs * exp(1im*δs[j]'*  K) for j in 1:3])
end

κ(G1, G2) =  1/3 * (2*G1+G2) 
m(G1, G2) = 1/2 * (G1+G2)


#_________________________________________________________________________________________
# 
#               PLOT BANDS
#
#_________________________________________________________________________________________

"""given a k path it computes the bands of the `bistritzer_hamiltonian()
along the linecut specified by `kpath` (i.e., K1-Γ-M-K2 line)"""
function bands_MLG(μ, ν; kws...)
    kmesh = kpath()
    sys_dim = 2
    es = zeros(Float64, sys_dim, size(kmesh,1))
    for i in 1:size(kmesh,1)
        es[:,i] = real.(MLG_eigs(μ, ν, kmesh[i,:])[1])
    end
    return (es, )
end

function MLG_eigs(μ, ν, q; kws...)
    l = eigen(Matrix(MLG_hamiltonian(μ, ν, q)))
    return l.values, l.vectors
end

"""given a mesh density with `pointsk` it computes the line K1, Γ, M, K2, where K1 and K2 
are the non-equiv point of the MBZ"""
function kpath(pointsk = 100)
    b1 = 2π * [1, 1/√3]
    b2 = 2π * [-1, 1/√3]
    κ(G1, G2) =  1/3 * (2*G1+G2) 
    m(G1, G2) = 1/2 * (G1+G2)
    κ1 = κ(b1, b2)
    κ2 = κ(b2, b1)
    meshdim = Int(5*pointsk/2 + 1)
    kpoints = zeros(Float64, meshdim, 2)
    for q = 1:meshdim
        if q <= pointsk
            qx = κ1[1] - κ1[1]*(q - 1)/pointsk
            qy = κ1[2] - κ1[2]*(q - 1)/pointsk
        elseif q > pointsk && q <= 2*pointsk
            qx = (κ2[1] + κ1[1])/2*(q - 1 - pointsk)/pointsk
            qy = (κ2[2] + κ1[2])/2*(q - 1 - pointsk)/pointsk
        else
            qx = (κ2[1] + κ1[1])/2 + (κ2[1] - κ1[1])*(q - 1 - 2*pointsk)/pointsk
            qy = (κ2[2] + κ1[2])/2 + (κ2[2] - κ1[2])*(q - 1 - 2*pointsk)/pointsk
        end
        kpoints[q,:] = [qx, qy]
    end
    return kpoints
end
    
######## Plots

function plotMLGbands(μ)
    mat, = bands_MLG(μ, 1)
    mat2, = bands_MLG(μ, -1)
    f = plotbands(mat, mat2)
    return f
end

plotbands(mat, mat2; kw...) = plotbands!(Figure(), mat, mat2; kw...)
function plotbands!(f, mat, mat2; dots = false, color = missing, ylimits = missing, xlimits = missing)
    ax = Axis(f[1, 1]; xlabel = "k", ylabel = "E [eV]")
    xarr = collect(1:size(mat,2))
    pointsk = 2/5 * (length(xarr)-1)
    if dots == false
        for i in 1:size(mat, 1)
            lines!(ax, xarr , mat[i,:], color = ifelse(isa(color,Missing), :lightgray, color))
            lines!(ax, xarr , mat2[i,:], color = ifelse(isa(color,Missing), :orange, color))
        end
    else
        for i in 1:size(mat, 1)
            scatter!(ax, collect(1:size(mat,2)) , mat[i,:], markersize = 5, markeralpha = 0.8)
            scatter!(ax, collect(1:size(mat,2)) , mat2[i,:], markersize = 5, markeralpha = 0.8)
        end
    end   
    ax.xticks = ([1, pointsk+1, 2*pointsk+1, 5*pointsk/2 + 1], ["K1", "Γ", "M", "K2"])
    if isa(ylimits,Missing)
        nothing
    else
        ylims!(ax, ylimits[1], ylimits[2])
    end
    if isa(xlimits,Missing)
        nothing
    else
        xlims!(ax, ylimits[1], ylimits[2])
    end
    return f
end
    