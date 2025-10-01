@with_kw struct Computation_presets
    xbounds::SVector{2, Float64} # Integration limits along k-space x-axis
    ybounds::SVector{2, Float64} # Integration limits along k-space y-axis
    ωlist::Array{Float64}        # List of frequencies for numerical evaluation
    broadening::Float64          # Phenomenological δ-function broadening
    evals::Number                # Max number of evaluations in adaptive integration
end

@with_kw struct DOS_presets
    h::Function # Hamiltonian H(k)
    computation::Computation_presets
end

@with_kw struct JDOS_presets
    h::Function # Hamiltonian H(k)
    computation::Computation_presets
end

@with_kw struct σij_presets
    dirJ::Symbol # i'th direction of σij
    dirE::Symbol # j'th direction of σij
    h::Function # Hamiltonian H(k)
    nabla_h::Function # ∇_k H(k)
    computation::Computation_presets
end