# 5. Introduction

*(One paragraph per subsection, per the guidelines.)*

**Background.** When light interacts with matter, its electromagnetic field
carries both an electric and a magnetic component. Conventional optics treats
the interaction almost exclusively through the electric field, and the optical
magnetic field is neglected as insignificant. Recent research — including work
by our supervisors — shows that this neglect is not always justified: the
optical magnetic field couples to the magnetization of magnetic materials and
produces measurable magneto-optical responses, which can become significant in
specific regimes such as circularly polarized excitation. Studying this
interaction computationally requires solving Maxwell's equations and the
magnetization dynamics as one coupled system, at a numerical precision fine
enough to resolve effects that are many orders of magnitude smaller than the
driving field.

**Objectives.** The goal of this project is to build a fast, validated
simulation framework for the coupled light–magnetization problem: (i) solve
Maxwell's equations and the Landau–Lifshitz–Gilbert (LLG) equation as a
two-way coupled system in one dimension; (ii) accelerate the computation on a
GPU so that production runs complete in minutes rather than hours; (iii) do so
at double precision, so that magneto-optical signatures at the 10⁻⁹ level of
the field are resolved; and (iv) validate the results against the previous
project's published values.

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
method and the GPU implementation. Chapter 9 presents the testing methodology,
validation, and performance results. Chapter 10 discusses the results and
their limitations, and Chapter 11 concludes.

---

# 6. Literature and State-of-the-Art Review

*(1–2 paragraphs per the guidelines; insert exact advisor references.)*

The finite-difference time-domain (FDTD) method on the staggered Yee grid
[Yee, 1966] is the standard direct solver for Maxwell's equations, and the
Landau–Lifshitz–Gilbert equation [Landau & Lifshitz, 1935; Gilbert, 2004]
is the established description of magnetization dynamics in ferromagnets.
Their combination into a mutually coupled FDTD–LLG scheme — where the optical
magnetic field drives the magnetization and the magnetization acts back on
the field — was demonstrated by Choi et al. [2017], whose one-dimensional
formulation and published observables form the methodological basis of this
work. Recent work on the role of the optical magnetic field in light–matter
interaction, including publications by our supervisors [**refs: Capua group /
Assouline — insert exact citations**], motivates resolving these couplings
quantitatively. The direct predecessor of this project is the previous
project by O. Stein and R. Kolodny, which implemented the coupled scheme
on CPU and published validated reflection/transmission and light-induced
magnetization values for the material systems considered here.

On the computational side, FDTD is a canonical GPU workload: the per-cell
update is local (nearest-neighbor stencil), making one-thread-per-cell
mappings natural, with performance bounded by memory bandwidth on large grids
and by launch overhead on small ones. However, the available GPU-FDTD
literature and tooling address the electromagnetic problem alone. To our
knowledge, no openly available implementation couples FDTD to the LLG
equation on GPU at double precision — the regime required for sub-10⁻⁹
magneto-optical signatures — which is the gap this project addresses.

---

# 7. Theoretical Background

*(~1 page per the guidelines.)*

**Maxwell's equations on the Yee grid.** In a source-free 1-D geometry with
propagation along z, Maxwell's curl equations reduce to coupled first-order
equations for the transverse components (Eₓ, E_y) and (Hₓ, H_y):

∂H/∂t = −(1/μ₀) ∂E/∂z − ∂M/∂t,  ∂E/∂t = −(1/ε₀ε_r)(∂H/∂z + σE)

with material permittivity ε_r and conductivity σ. Discretizing on the Yee
grid places E and H on staggered half-cells and half-time-steps, so each
update uses a two-point spatial stencil. Stability requires the Courant
condition; we use Δt = Δz/2c, for which the vacuum update coefficient is 0.5
(Sullivan Eq. 1.9):

H_x[k] += 0.5 (E_y[k] − E_y[k−1]),  H_y[k] −= 0.5 (E_x[k] − E_x[k−1])
E_y[k] = ca[k] E_y[k] + cb[k] (H_x[k+1] − H_x[k]),  E_x[k] = ca[k] E_x[k] − cb[k] (H_y[k+1] − H_y[k])

Here H[k] sits half a cell to the left of E[k]. The coefficients
ca = (1−η)/(1+η) and cb = 0.5/(ε_r(1+η)), with η = σΔt/(2ε₀ε_r), add the
permittivity and conduction losses of the material (Sullivan Eq. 1.23–1.24);
in vacuum ca = 1 and cb = 0.5. The grid ends are absorbing: since a wave
crosses one cell in exactly two steps, the edge value is the neighbour's value
from two steps earlier, E^n(0) = E^{n−2}(1) (Mur boundary, Sullivan Eq. 1.12).

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
On the Yee time grid M lives with H, and each RK4 step is driven by the
current field H^n.

**Two-way coupling.** The optical magnetic field enters the LLG equation as
the driving field H; in turn, the magnetization enters Faraday's law through
B = μ₀(H + M), so a change of M during a step feeds back into H. It is this
back-action that produces the magneto-optical observables studied here: the
difference between a simulation with the back-action enabled ("closed") and
disabled ("open") isolates the effect of the magnetization on the light.

**Magnitude of the effect and precision.** For the film thicknesses and
pulse energies considered, the open−close field difference is of order 10⁻⁹
of the driving amplitude. Resolving a 10⁻⁹ relative signal after 10⁷
integration steps requires double precision: float32 arithmetic accumulates
rounding at the 10⁻⁷ level over such runs and buries the signal (quantified
in Section 9.5).

**Polarization eigenmodes.** In a magnetized medium the circular polarization
states σ± are the propagation eigenmodes, with slightly different refractive
indices (circular birefringence). Linearly polarized light — an equal
superposition of σ+ and σ− — therefore experiences Faraday rotation of its
polarization plane, at a rate proportional to the index difference. Purely
circular excitation, as used in our presets, experiences the same physics as
a phase shift of the driven eigenmode rather than a rotation; both
observables are extracted in Section 9.5.
