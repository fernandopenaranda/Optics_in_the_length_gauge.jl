# Optics_in_the_length_gauge

<!-- [![Build Status](https://github.com/fernandopenaranda/Optics_in_the_length_gauge.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fernandopenaranda/Optics_in_the_length_gauge.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

Given a k space Hamiltonian, $H(\vec k)$, it computes several optical responses in the length gauge:

+ Joint density of states
  + $\text{jdos}(\omega) = \omega f_{nm} \delta(\omega-\omega_{mn})$
  + 2D implementation

+ Linear optical response $\vec{J}_a = σ_{ab}(ω) \vec{J}_b$ (absortive)
   + $\sigma_{ab}^{\text{abs}}(\omega) = \frac{\pi e^2}{\hbar} \sum_{n,m} \omega f_{nm} r_{nm}^a r_{mn}^b \delta(\omega-\omega_{mn})$
   + 2D implementation

    where $a$ and $b$ denote spatial directions unbounded dimensions 
+ ...

Input: $H(\vec k)$, $\vec \partial H(\vec k)/\partial \vec k$, and the boundaries of the BZ (the latter two in improved versions will be computed automatically).

Comments: $\delta(\omega-\omega_{mn})$ is approximated by a Lorentzian broadened by $\Gamma$: $\frac{\Gamma}{(\omega-\omega_{mn})^2 + \Gamma^2}$