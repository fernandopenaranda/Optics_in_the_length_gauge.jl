# Optics_in_the_length_gauge

<!-- [![Build Status](https://github.com/fernandopenaranda/Optics_in_the_length_gauge.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fernandopenaranda/Optics_in_the_length_gauge.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

Package to compute several spectral quantities and optical responses of general k-space Hamiltonians in the length gauge, including:

+ Density of states

+ Joint density of states
  + $\text{jdos}(\omega) = \sum_{nm} f_{nm} \delta(\omega-\omega_{mn})$

+ Linear optical responses $\vec{J}_a = σ_{ab}(ω) \vec{J}_b$ (absortive)
   + $\sigma_{ab}^{\text{abs}}(\omega) = \frac{\pi e^2}{\hbar} \sum_{n,m} \int \frac{d\vec k}{(2\pi)^d}(-\omega) f_{nm} r_{nm}^a r_{mn}^b \delta(\omega-\omega_{mn})$
   + 2D implementation -> Generalize to unbounded 3D (trivial)

    where $a$ and $b$ denote spatial directions unbounded dimensions 

+ Non Linear optical responses [updated in future versions]

And transport quantities:

+ Magneto conductivity
+ ...
  
Input: the Hamiltonian $H(\vec{k})$, its gradient $\nabla_{\vec{k}} H(\vec{k})$, and the Brillouin zone boundaries (the latter two will be computed automatically in future versions).

Comments: 

+ Input units: $\omega$ in eV, BZ boundaries: adimensional (multiplied by $a_0$ the lattice constant)
+ $\delta(\omega-\omega_{mn})$ is approximated by a Lorentzian broadened by $\Gamma$: $\frac{\Gamma}{(\omega-\omega_{mn})^2 + \Gamma^2}$
+ IMPORTANT, the adaptive integration in 2D assumes a rectangular BZ, take this in mind when passing the bounds `xbounds` and `ybounds`. Improvement needed.