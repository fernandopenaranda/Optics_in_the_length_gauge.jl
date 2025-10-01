module MLGPresets
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

    function MLG_deriv(q, dir::Symbol)
        ind = (dir == :x ? 1 : 2)
        st = _d_sublattice_coupling(q)
        return [0 st[ind]; conj(st[ind]) 0]
    end

    MLG_nabla(q) = [MLG_deriv(q, :x), MLG_deriv(q, :y)]

    _sublattice_coupling(q) = t1 * sum([exp(1im*δs[j]'*  q)  for j in 1:3])
    """
    [∂H/∂kx, ∂H/∂ky] for MLG hamiltonian along the :x and .:y directions.
    """
    _d_sublattice_coupling(q) = 1im * t1 .* [sum([δs[j][1] * exp(1im*δs[j]'*  q) for j in 1:3]), 
        sum([δs[j][2] * exp(1im*δs[j]'*  q) for j in 1:3])]
    export MLG_hamiltonian, MLG_nabla, K1, K2
end

