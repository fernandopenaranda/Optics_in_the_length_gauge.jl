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
    