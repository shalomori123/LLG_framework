# 5. Introduction

*(One paragraph per subsection, per the guidelines.)*

**Background.** When light interacts with matter, its electromagnetic field
carries both an electric and a magnetic component. Conventional optics treats
the interaction almost exclusively through the electric field, and the optical
magnetic field is neglected as insignificant. Recent work by our supervisors
shows that this neglect is not always justified: the optical magnetic field
alone, entered into the Landau–Lifshitz–Gilbert equation, produces a
helicity-dependent torque on the magnetization [1], and it contributes to the
Faraday and inverse Faraday effects [2]. Studying this interaction
computationally requires solving Maxwell's equations and the magnetization
dynamics as one coupled system, with numerical noise far below effects that
are many orders of magnitude smaller than the driving field.

**Objectives.** The goal of this project is to build a fast, validated
simulation framework for the coupled light–magnetization problem: (i) solve
Maxwell's equations and the Landau–Lifshitz–Gilbert (LLG) equation as a
two-way coupled system in one dimension; (ii) accelerate the computation on a
GPU so that production runs complete in minutes rather than hours; (iii)
resolve magneto-optical signatures at the 10⁻⁹ level of the field, choosing
the arithmetic precision by measuring its accuracy–performance trade-off; and
(iv) validate the results against the lab's previous CPU method and its
reference results.

**Scope.** The framework treats one-dimensional layered geometries — a thin
magnetic film (Co, Fe, or Ni) with a conductive cap (Au or Pt) in vacuum —
excited by a femtosecond–picosecond Gaussian pulse at λ = 784 nm. Two
simulation presets are provided: a subwavelength high-resolution grid for
in-material magnetization dynamics, and a propagation-scale grid for
reflection/transmission analysis. The implementation targets a single NVIDIA
GPU. Multi-dimensional geometries, thermal effects, and quantum corrections
are outside the scope of this work.

**Organization.** Chapter 6 reviews the relevant literature. Chapter 7
summarizes the theoretical background. Chapter 8 describes the numerical
method, the GPU implementation and how it developed. Chapter 9 presents the
testing methodology, validation, and performance results. Chapter 10
discusses the results and their limitations, and Chapter 11 concludes.

---

# 6. Literature and State-of-the-Art Review

*(1–2 paragraphs per the guidelines.)*

The finite-difference time-domain (FDTD) method is the standard direct solver
for Maxwell's equations in the time domain; our formulation follows
Sullivan's textbook treatment [3]. The Landau–Lifshitz–Gilbert equation
[4, 5] is the established description of magnetization dynamics in
ferromagnets. Their combination into a mutually coupled FDTD–LLG scheme —
where the optical magnetic field drives the magnetization and the
magnetization acts back on the field — was demonstrated by Choi et al. [6],
whose one-dimensional formulation is the methodological basis of this work.
Our supervisors showed that for ultrashort, circularly polarized pulses the
optical magnetic field in the LLG equation produces a helicity-dependent
torque whose strength scales as αγH/f_opt [1], and that the magnetic
component of light also contributes to the Faraday effect, accounting for
about 17% of the measured Verdet constant of terbium-gallium-garnet at
800 nm [2]. These results motivate resolving the coupling quantitatively in
layered metallic films. Within the lab, the coupled scheme was previously
implemented as a CPU (NumPy) simulation, which produced reference
reflection/transmission and light-induced magnetization values for the
material systems considered here; that method is our baseline and
validation target.

On the computational side, FDTD is a canonical GPU workload: the per-cell
update is local (nearest-neighbour stencil), making one-thread-per-cell
mappings natural, with performance bounded by memory bandwidth on large grids
and by launch overhead on small ones. However, the available GPU-FDTD
literature and tooling address the electromagnetic problem alone. To our
knowledge, no openly available implementation couples FDTD to the LLG
equation on the GPU with the accuracy needed for sub-10⁻⁹ magneto-optical
signatures, which is the gap this project addresses.

---

# 7. Theoretical Background

*(~1 page per the guidelines.)*

**The FDTD method (after Sullivan [3], Ch. 1).** For a plane wave
propagating along z, Maxwell's curl equations reduce to two scalar equations
per polarization. For (E_x, H_y):

∂E_x/∂t = −(1/ε₀ε_r) ∂H_y/∂z − (σ/ε₀ε_r) E_x,  ∂H_y/∂t = −(1/μ₀) ∂E_x/∂z

and the same pair, with the opposite sign, couples (E_y, H_x). Both the time
and the space derivatives are replaced by central differences. For these to
be centred, E and H are *interleaved in space and time*: E is evaluated at
whole cells k and half time steps, H at half cells k ± ½ and whole time
steps. Each field is then advanced from its own previous value and the two
neighbouring values of the other field, so the simulation alternates: all H
values, then all E values, then H again.

Sullivan normalizes the electric field, Ẽ = √(ε₀/μ₀) E, so that both
equations carry the same coefficient c₀Δt/Δz. Stability requires that a wave
not cross more than one cell per time step (the Courant condition). Like
Sullivan, we choose Δt = Δz/(2c₀), which makes the vacuum coefficient exactly
½. In our indexing H[k] sits half a cell to the left of E[k]:

H_x[k] += ½ (E_y[k] − E_y[k−1]),  H_y[k] −= ½ (E_x[k] − E_x[k−1])
E_y[k] = ca[k] E_y[k] + cb[k] (H_x[k+1] − H_x[k]),  E_x[k] = ca[k] E_x[k] − cb[k] (H_y[k+1] − H_y[k])

In a lossy dielectric the conduction term is averaged over the two time
levels, which gives the coefficients ca = (1−η)/(1+η) and
cb = ½/(ε_r(1+η)), with η = σΔt/(2ε₀ε_r) (Sullivan's lossy-medium section).
In vacuum ca = 1 and cb = ½; in the metal layers they make the film reflect
and absorb.

The grid must end without reflecting. With Δt = Δz/(2c₀) a wave needs
exactly two time steps to cross one cell, so the edge value is the
neighbour's value from two steps earlier, E^n(0) = E^{n−2}(1), and likewise
at the far end (Sullivan's one-dimensional absorbing boundary condition).

**The LLG equation.** The magnetization M of a ferromagnetic cell precesses
around the local magnetic field H and relaxes toward it with the
dimensionless Gilbert damping α:

dM/dt = −γμ₀ (M × H) + (α/|M|) (M × dM/dt)

where γ is the gyromagnetic ratio. We integrate the equivalent explicit
(Landau–Lifshitz) form, which has no dM/dt on the right-hand side:

dM/dt = −γ_LL μ₀ (M × H) − (λ μ₀/|M|) M × (M × H),  γ_LL = γ/(1+α²), λ = α γ_LL

Two properties matter here. First, the equation preserves |M| exactly (both
terms are perpendicular to M), which provides a sharp numerical invariant.
Second, under fast optical driving it is stiff, which motivates fourth-order
Runge–Kutta integration rather than the leapfrog scheme used for the fields.
M is stored at the positions of H, and each RK4 step is driven by the current
field H^n.

**Two-way coupling.** The optical magnetic field enters the LLG equation as
the driving field H; in turn, the magnetization enters Faraday's law through
B = μ₀(H + M), so a change of M during a step feeds back into H. It is this
back-action that produces the magneto-optical observables studied here: the
difference between a simulation with the back-action enabled ("closed") and
disabled ("open") isolates the effect of the magnetization on the light.

**Magnitude of the effect.** For the film thicknesses and pulse energies
considered, the open−closed field difference is of order 10⁻⁹ of the driving
amplitude. This sets how much numerical noise a run may accumulate over 10⁷
steps, and therefore which arithmetic precision is adequate for which
observable. Section 8.8 describes how precision is selected, and Section 9.5
compares the two choices.

**Polarization eigenmodes.** In a magnetized medium the circular polarization
states σ± are the propagation eigenmodes, with slightly different refractive
indices (circular birefringence). Linearly polarized light — an equal
superposition of σ+ and σ− — therefore experiences Faraday rotation of its
polarization plane, at a rate proportional to the index difference. Purely
circular excitation, as used in our presets, experiences the same physics as
a phase shift of the driven eigenmode rather than a rotation; both
observables are extracted in Section 9.5.
