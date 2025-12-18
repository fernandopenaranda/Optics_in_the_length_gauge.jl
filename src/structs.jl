# UNITS E -> eV, time -> s, length -> Å, temperature -> K

@with_kw struct Optical_computation_presets
    xbounds::SVector{2, Float64} # Integration limits along k-space x-axis
    ybounds::SVector{2, Float64} # Integration limits along k-space y-axis
    ωlist::Array{Float64}        # List of frequencies for numerical evaluation
    broadening::Float64          # Phenomenological δ-function broadening
    evals::Number                # Max number of evaluations in adaptive integration
end

@with_kw struct Transport_computation_presets
    xbounds::SVector{2, Float64} # Integration limits along k-space x-axis
    ybounds::SVector{2, Float64} # Integration limits along k-space y-axis
    evals::Number                # Max number of evaluations in adaptive integration
end

@with_kw struct DOS_presets
    a0::Float64
    h::Function # Hamiltonian H(k)
    computation::Union{Optical_computation_presets, Transport_computation_presets}
end

@with_kw struct JDOS_presets
    a0::Float64
    h::Function # Hamiltonian H(k)
    computation::Union{Optical_computation_presets, Transport_computation_presets}
end

@with_kw struct σij_presets
    a0::Float64
    dirJ::Symbol # i'th direction of σij
    dirE::Symbol # j'th direction of σij
    h::Function # Hamiltonian H(k)
    nabla_h::Function # ∇_k H(k)
    computation::Union{Optical_computation_presets, Transport_computation_presets}
end

# TRANSPORT 

@with_kw struct Drude_presets
    a0::Float64
    dirJ::Symbol # i'th direction of σij
    dirE::Symbol # j'th direction of σij
    h::Function # Hamiltonian H(k)
    dhi::Function # k-derivative of H(k) in dirE
    T::Float64
    τ::Float64
    computation::Union{Optical_computation_presets, Transport_computation_presets}
end

@with_kw struct AH_presets
    a0::Float64
    dirJ::Symbol # i'th direction of σij in i = {x, y}
    dirE::Symbol # j'th direction of σij in j = {y, x}
    h::Function # Hamiltonian H(k)
    dh::Function # x, y, k-derivative of H(k): [dhx, dhy]
    T::Float64
    computation::Union{Optical_computation_presets, Transport_computation_presets}
end

@with_kw struct Planar_σijk_presets_orbital
    a0::Float64
    dirJ::Symbol # i'th direction of σijk
    dirE::Symbol # j'th direction of σijk
    dirB::Symbol # k'th direction of σijk
    h::Function # Hamiltonian H(k)
    nabla_h::Function # ∇_k H(k)
    nabla_nabla_h::Function #∂^2/∂i^2 H(k) with i = x or y (planar)
    rz::Function # even thought it is generally independent on q, pass it in the form rz(q, ψ)
    τ::Float64 # scattering time
    T::Number # Temperature in K
    computation::Union{Optical_computation_presets, Transport_computation_presets}
    berry_contribution = true # berry_contribution to the Planar_σijk
    omm_contribution = true # omm contribution to the Planar_σijk
    fermi_surface = false # computes the jdos if true
    with_shift = true # corrects the Planar_σijk with the constant factor that comes from the B field dependency in the bands
end