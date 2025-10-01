
const t1 = -2.46575 # hopping amplitude
const δ1 = [1/2, √3/2]
const δ2 = [1/2, -√3/2]
const δ3 =  [-1,0]
const δs = [δ1, δ2, δ3] #real space vectors connecting A and B sites
const K1 = 2π/3 * [1, 1/√3] # K point
const K2 = 2π/3 * [1,-1/√3] # K' point

"""
 Monolayer graphene Hamiltonian μ is the chemical potential in eV, the two valleys are considered, q::Array,
 is the momentum (in adimensional units). Spin degeneracy is not included.
"""
function MLG_hamiltonian(μ, q)
    st = _sublattice_coupling(q)
    return [-μ st; conj(st) -μ]
end

function MLG_hamiltonian_derivatives(q)
    st = _d_sublattice_coupling(q)
    return [0 st[1]; conj(st[1]) 0], [0 st[2]; conj(st[2]) 0]
end

_sublattice_coupling(q) = t1 * sum([exp(1im*δs[j]'*  q)  for j in 1:3])
"""
[∂H/∂kx, ∂H/∂ky] for MLG hamiltonian along the :x and .:y directions.
"""
_d_sublattice_coupling(q) = 1im * t1 .* [sum([δs[j][1] * exp(1im*δs[j]'*  q) for j in 1:3]), 
    sum([δs[j][2] * exp(1im*δs[j]'*  q) for j in 1:3])]

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
    [es[:,i] = real.(MLG_eigs(μ, ν, kmesh[i,:])[1]) for i in 1:size(kmesh,1)]
    return (es, )
end

function MLG_eigs(μ, ν, q; kws...)
    l = eigen(Matrix(MLG_hamiltonian(μ, q)))
    return l.values, l.vectors
end

"""given a mesh density with `pointsk` it computes the line K1, Γ, M, K2, where K1 and K2 
are the non-equiv point of the MBZ"""
function kpath(pointsk = 100)   
    meshdim = Int(5*pointsk/2 + 1)
    kpoints = zeros(Float64, meshdim, 2)
    for q = 1:meshdim
        if q <= pointsk
            qx = K1[1] - K1[1]*(q - 1)/pointsk
            qy = K1[2] - K1[2]*(q - 1)/pointsk
        elseif q > pointsk && q <= 2*pointsk
            qx = (K2[1] + K1[1])/2*(q - 1 - pointsk)/pointsk
            qy = (K2[2] + K1[2])/2*(q - 1 - pointsk)/pointsk
        else
            qx = (K2[1] + K1[1])/2 + (K2[1] - K1[1])*(q - 1 - 2*pointsk)/pointsk
            qy = (K2[2] + K1[2])/2 + (K2[2] - K1[2])*(q - 1 - 2*pointsk)/pointsk
        end
        kpoints[q,:] = [qx, qy]
    end
    return kpoints
end
    
######## Plots

function plotMLGbands(μ)
    mat, = bands_MLG(μ, 1)
    f = plotbands(mat)
    return f
end

plotbands(mat; kw...) = plotbands!(Figure(), mat; kw...)
function plotbands!(f, mat; dots = false, color = missing, ylimits = missing, xlimits = missing)
    ax = Axis(f[1, 1]; xlabel = "k", ylabel = "E [eV]")
    xarr = collect(1:size(mat,2))
    pointsk = 2/5 * (length(xarr)-1)
    if dots == false
        for i in 1:size(mat, 1)
            lines!(ax, xarr , mat[i,:], color = ifelse(isa(color,Missing), :orange, color))
        end
    else
        for i in 1:size(mat, 1)
            scatter!(ax, collect(1:size(mat,2)) , mat[i,:], markersize = 5, markeralpha = 0.8)
         
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
    