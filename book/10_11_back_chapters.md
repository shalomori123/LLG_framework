# 10. Discussion

*(2–3 paragraphs per the guidelines.)*

The coupled Maxwell–LLG problem suits the GPU well, but where the speedup comes
from depends on the grid size. On the small subwavelength grid one simulation
uses a single block, and the step time is mostly launch overhead. The gain
comes from replacing many per-step array operations with two kernels, and
batching independent cases along `blockIdx.y` uses the rest of the GPU. On the
large propagation grid the run is limited by memory bandwidth, and the speedup
approaches the GPU/CPU bandwidth ratio. The design followed from these two
regimes:

- the source and boundary moved inside the kernels to cut launches;
- an owner-thread rule made that safe on any number of blocks, so one kernel
  serves both grids;
- float64 became the default, so the 10⁻⁹ magneto-optical signal is resolved
  on the GPU while |M| is conserved to 2.35×10⁻¹³ over 2×10⁷ steps.

A second lesson concerns verification. Developing two GPU implementations in
parallel and diffing them on identical inputs revealed seven discrepancies
(Section 8.9) that neither implementation exposed on its own. Some of them,
such as a sign error in the Yee update or a coupling switch that was always on,
change the physics without making the program fail.

Some discrepancies remain and should be stated plainly. The
reflection/transmission validation reproduces the published values to their
full reported precision, but part of the light-induced magnetization table
agrees only to 71–93%: the previous team's exact extraction script was not
preserved, so our reconstruction of their post-processing (averaging windows
and normalization) is a plausible source of the residual gap, alongside minor
configuration drift accumulated over the two projects. Additionally, the
open−close energy-difference labels (ΔR vs ΔT) appear swapped relative to the
published convention — the magnitudes match at ~3×10⁻⁹ but the attribution
awaits a convention check. Finally, the full-length CPU validation of the
large-grid preset was truncated (30,000 steps) for practicality; the complete
~16-hour CPU run was never executed, so the large-grid CPU–GPU equivalence is
verified on the truncated window only.

Several improvements follow naturally:

- *Small-grid launch overhead.* CUDA Graphs (recording many steps and
  replaying them with one launch) or a persistent kernel that loops over time
  steps inside one block would remove most of the remaining launch cost.
  The device-only curve of Figure 9.7, measured by exactly such a graph
  replay, is the upper bound of this gain.
- *Validation tables.* Running the remaining material combinations, now a
  single batched run, would complete them.
- *Faraday rotation.* A linearly polarized preset (a one-line change in
  `sim_config.py`) would measure the rotation angle directly.
- *Higher dimensions.* The per-cell structure extends to 2-D and 3-D
  geometries, where the GPU's advantage grows with the cell count.

---

# 11. Conclusions

*(One paragraph per the guidelines.)*

We built and validated a GPU-accelerated framework for the mutually coupled
Maxwell–LLG equations, packaged as four files: one CUDA kernel source for
both precisions, any grid size and batched cases; a host loop with bounded
memory use; a configuration file; and a command-line runner. The implementation
reproduces the previous project's published reflection/transmission values to
five decimal places, conserves the magnetization magnitude to the
double-precision machine floor over twenty million integration steps, and —
using custom float64 CUDA kernels — resolves magneto-optical signatures at
10⁻⁹ of the driving field directly on the GPU, at speedups of 15–180×
over the CPU reference depending on problem size. Beyond the specific physics, the project's central
lesson is the interplay between hardware architecture and numerical physics:
the correct precision, memory layout, and synchronization strategy were as
essential to observing the physical effect as the equations themselves.

---

# 4. List of Abbreviations and Symbols

*(Place after the table of contents.)*

**Abbreviations**

| Abbrev. | Meaning |
|---|---|
| ABC | Absorbing Boundary Condition (first-order Mur) |
| NVRTC | NVIDIA runtime compiler (kernels compiled for the installed GPU) |
| SM | Streaming Multiprocessor (GPU core cluster that runs blocks) |
| CPU / GPU | Central / Graphics Processing Unit |
| CUDA | NVIDIA's GPU programming platform |
| FDTD | Finite-Difference Time-Domain |
| fp32 / fp64 | 32-/64-bit IEEE floating point |
| FWHM | Full Width at Half Maximum |
| LLG | Landau–Lifshitz–Gilbert (equation) |
| MO | Magneto-Optical |
| R / T | Reflection / Transmission energy ratio |
| RK4 | 4th-order Runge–Kutta |
| SoA | Structure of Arrays (memory layout) |

**Symbols**

| Symbol | Meaning | Units |
|---|---|---|
| E, H | Electric / magnetic field | V/m, A/m |
| M, Ms | Magnetization, saturation magnetization | A/m |
| α | Gilbert damping factor | — |
| γ | Gyromagnetic ratio | rad·s⁻¹·T⁻¹ |
| λ | Wavelength (784 nm) | m |
| σ | Electrical conductivity | S/m |
| ε_r | Relative permittivity | — |
| Δz, Δt | Spatial / temporal step | m, s |
| σ± | Circular polarization eigenmodes | — |
| θ_F | Faraday rotation angle | rad |

---

# 13. References (skeleton — IEEE style)

1. K. S. Yee, "Numerical solution of initial boundary value problems
   involving Maxwell's equations in isotropic media," *IEEE Trans. Antennas
   Propag.*, vol. 14, no. 3, pp. 302–307, 1966.
2. L. Landau and E. Lifshitz, "On the theory of the dispersion of magnetic
   permeability in ferromagnetic bodies," *Phys. Z. Sowjetunion*, vol. 8,
   pp. 153–169, 1935.
3. T. L. Gilbert, "A phenomenological theory of damping in ferromagnetic
   materials," *IEEE Trans. Magn.*, vol. 40, no. 6, pp. 3443–3449, 2004.
4. D. M. Sullivan, *Electromagnetic Simulation Using the FDTD Method*,
   2nd ed. Hoboken, NJ: Wiley–IEEE Press, 2013. **[verify edition used]**
5. **[Choi et al., 2017 — insert exact citation of the coupled FDTD–LLG
   paper used]**
6. **[Advisors' publications — request exact list from Prof. Capua /
   Dr. Assouline]**
7. O. Stein and R. Kolodny, *[previous project book title]*, final project
   book, The Hebrew University of Jerusalem, [year].
8. NVIDIA Corp., *CUDA C++ Programming Guide*, [version/year accessed].
9. R. Okuta et al., "CuPy: A NumPy-compatible library for NVIDIA GPU
   calculations," *Proc. 31st NeurIPS Workshop*, 2017.
10. NVIDIA Corp., *Nsight Systems* and *Nsight Compute* user guides,
    [version/year accessed].

---

# 14. Appendices (planned contents)

- **A. Simulation configurations** — the materials table and the two presets
  from `sim_config.py` (all physical and numerical inputs, plus the formulas for
  derived quantities).
- **B. Output bundle format** — the `.simdata` manifest structure and a
  loading example.
- **C. Reproduction guide** — for every figure and table in Chapter 9: the
  script and command that regenerates it (draft below).
- **D. Repository guide** — condensed README: layout, install, quickstart.
- **E. Supplementary results** — run overviews and videos (list below).

**C. Reproduction guide (draft).** [U] = unified framework repository
(`LLG_framework`), [A] = analysis repository
(`llg_framework_cpu_version_updated`). Bundles are produced first by
`python gpu_tools/sim_run.py --preset <full_width|rt> --state both` [U].

| Item | Script | Command / notes |
|---|---|---|
| 9.2 equivalence | `validation/compare_with_partner.py` [U] | `--precision fp64`, then `fp32` |
| Figure 9.1 | `analysis/advisor_view_compare.py` [A] | `<GPU Co_Pt_close.simdata> <orig_data/Co_Pt_close.npy>` |
| Table 9.1 | `comparison/compute_table3_RT.py` [A] | `rt` bundles, Co/Pt and Ni/Pt |
| Figure 9.2 | `analysis/transmit_reflect_analysis.py` [A] | `rt` Co/Pt bundle |
| Table 9.2 | `comparison/compute_table2_mz_my.py` [A] | `full_width` open + closed bundles |
| Figure 9.3 | `analysis/check_M_norm.py` [A] | numbers; plot = Figure 9.1 panel C as lines |
| Figure 9.4 | `analysis/make_poster_mo_signal.py` [A] | `full_width` open + closed, float64 |
| Figure 9.5 | `analysis/faraday_rotation.py` [A] | open + closed bundles |
| Figure 9.6 | `analysis/run_ms_sweep.py`, `analysis/compare_open_close_sweep.py` [A] | ~20 min per Ms scale |
| Table 9.3 | `benchmarks/profile_modes.py` [U] | full-length times: `gpu_tools/sim_run.py` runs, timed |
| Table 9.4 | `benchmarks/profile_modes.py` [U] | `--presets full_width rt --precisions fp64 --batches 1 --nsys --ncu` → `report.md` |
| Figures 9.7–9.8 | `benchmarks/scaling_sweep.py` [U] | default sweep 78 → 10⁶ cells; `--plot-only <dir>` to restyle |

**E. Supplementary results (planned).**

- One-page overview of every final run (12 materials × open/closed):
  `analysis/advisor_view.py` [A] on each bundle.
- Figure 9.1 for the open system: `analysis/advisor_view_compare.py` [A] on the
  `_open` bundles.
- Videos (linked or QR code, not printable): wave propagation
  `analysis/make_video.py`, magnetization profile `analysis/make_M_video.py`,
  open−closed field difference `analysis/make_open_close_diff_video.py` [A].
- Full mode profile tables (all presets × precisions × batch sizes):
  `report.md` of `benchmarks/profile_modes.py` [U].
