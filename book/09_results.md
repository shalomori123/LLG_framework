# 9. Testing and Results

*(All graphs are line plots without titles, per the guidelines. Figures to be
regenerated from final bundles of the unified framework. Items marked
**[re-run]** were measured with the per-preset backends of track B and must be
re-measured with the unified framework before submission; the physics is
identical, see 9.2. Each figure placeholder names the script that produces
it; Appendix C lists the exact commands. Several analysis scripts draw panel
titles, a suptitle or bars for quick review; the book versions drop the titles
and redraw bars as lines.)*

## 9.1 Testing procedures

Testing ran at three levels, from the arithmetic up to the physics.

1. **Implementation equivalence.** All six output arrays of a run (field time
   series at probe locations, full-grid snapshots, magnetization) were diffed
   element by element:
   - GPU against the CPU method, over full-length runs;
   - each optimization stage against the stage before it;
   - the unified framework against the per-preset backends of track B
     (`validation/compare_with_partner.py`).
2. **Physical invariants.** |M| must be conserved by the LLG equation, and
   R + T ≈ 1 must hold for energy. The building blocks (pulse, vacuum
   propagation, dielectric and conductive interfaces, free precession) were
   verified separately in notebooks.
3. **Lab reference results.** The reflection/transmission and magnetization
   values obtained with the lab's previous method, for Co/Pt and Ni/Pt.

## 9.2 Numerical correctness

**GPU vs CPU.** Over the full 2×10⁷-step `full_width` run, the float32 GPU
backend matches the float64 CPU method at the float32 floor: maximum
relative deviation ~5×10⁻⁶ in E, ~2×10⁻⁶ in H, 10⁻⁷–10⁻⁶ in M, with no NaN/Inf
and no drift. The two optimization stages (device-side output buffers;
boundaries inside the kernels) gave byte-identical output to their
predecessors. They moved work without changing the arithmetic.

**Unified framework.** The unified kernel is compared against the per-preset
backends on the same inputs (Co/Pt, open and closed, both presets). The
expected agreement is ~10⁻⁵ at float32 and ~10⁻¹² at float64. **[re-run]**

**The whole run at a glance.** Figure 9.1 overlays our GPU run (solid) on the
lab's previous CPU method (dashed) for Co/Pt in the closed system, with one
color per series. Four observables are shown together: the wave crossing the
grid, the precession it drives, the energy flow at the probes and the final
magnetization profile. If the two runs agree, each solid curve covers its
dashed twin. **[re-run]**

> **[Figure 9.1]** GPU (solid) vs the lab's previous CPU method (dashed),
> Co/Pt closed, `full_width`. (a) H_x(z) at five times [A/m vs cell index];
> (b) M_y(t), M_z(t) at the first magnetic cell [A/m vs ps]; (c) Poynting
> S_z(t) at the probes [W/m² vs ps]; (d) final M_y(z), M_z(z) in the magnetic
> layer [A/m vs nm]. The open system is in Appendix E.
> *Source:* `analysis/advisor_view_compare.py` (panels A, B, D, E). Book
> version: no panel titles or suptitle, legend label "lab previous method";
> panel C (per-cell |M| drift) becomes the line plot of Figure 9.3.

## 9.3 Validation against the lab's reference results

**Table 9.1 — Reflection (R) and transmission (T) energy ratios.**

| Material | R (lab reference) | R (this work) | T (lab reference) | T (this work) | R+T |
|---|---|---|---|---|---|
| Co/Pt | 0.91266 | 0.91267 | 0.08717 | 0.08717 | 0.9998 |
| Ni/Pt | 0.90383 | 0.90384 | 0.09598 | 0.09599 | 0.9998 |

Both systems agree to the five decimals of the reference values, in float32
and in float64 alike. The two systems differ in damping (α = 0.025 vs 0.06)
and conductivity, so the agreement does not rest on a single parameter set.

R and T are energy ratios: the Poynting flux S_z = E_x·H_y − E_y·H_x is
integrated in time at a probe before the material (incident and reflected
pulses, separated in time by the vacuum gap) and at a probe after it
(transmitted pulse), as in Figure 9.2.

> **[Figure 9.2]** S_z(t) at the entry and exit probes, Co/Pt, `rt` preset
> [W/m² vs ps], with the integration windows (t₁, t₂) of the incident,
> reflected and transmitted pulses marked.
> *Source:* `analysis/transmit_reflect_analysis.py`; numbers of Table 9.1 from
> `comparison/compute_table3_RT.py`.

**Table 9.2 — Magnetization observables, Co/Pt (`full_width`).**

| Observable | Lab reference | This work | Agreement |
|---|---|---|---|
| m_z/M (closed) ×10⁻⁵ | −0.0249 | −0.0232 (CPU float64) | 93% |
| Δ(m_y/M) open−closed ×10⁻⁹ | −2.237 | −2.057 (GPU float32) | 92% |

The remaining magnetization observables agree to 70–93%. Chapter 10 discusses
the gap: the post-processing script of the previous method was not
preserved, so part of it comes from reconstructing that post-processing.

## 9.4 Conservation of |M|

After 2×10⁷ steps at float64, the largest per-cell deviation
|M(z, t_final)|/Ms − 1 is 2.35×10⁻¹³, the double-precision floor (Figure 9.3).
RK4 with the per-cell |M| in the damping term keeps the LLG invariant over the
whole run.

> **[Figure 9.3]** Per-cell relative drift |M(z, t_final)|/Ms − 1 vs magnetic
> cell index, GPU and the lab's previous CPU method (two lines).
> *Source:* `analysis/check_M_norm.py` (numbers, worst transient per run) and
> panel C of `analysis/advisor_view_compare.py`, redrawn as lines.

## 9.5 The magneto-optical signal and the precision trade-off

The physical observable is the difference between the open system (M does
not act back on the field) and the closed system (two-way coupling). The
magnetization imprints a field difference ΔH_x(z) of order 10⁻⁹–10⁻¹⁰ of the
driving amplitude, growing across the magnetic layer (Figure 9.4).

The preset drives pure σ⁻ circular light, an eigenmode of the magnetized
medium. The effect therefore appears as a phase shift of that eigenmode rather
than a rotation, accumulating to ≈5×10⁻⁹ rad across the layer and constant
downstream (Figure 9.5).

> **[Figure 9.4]** ΔH_x(z) = H_x^open − H_x^closed at a fixed time [A/m vs
> cell index]; material region shaded.
> *Source:* `analysis/make_poster_mo_signal.py` (book version: no title).
>
> **[Figure 9.5]** Eigenmode phase shift Δarg F₋(z) [rad vs nm].
> *Source:* `analysis/faraday_rotation.py`.

**float32 or float64.** Because the signal is so small, the same runs were
made at both precisions (Section 8.8) and compared observable by observable
(Table 9.3). Rounding in float32 accumulates to ~10⁻⁷ of the field over
2×10⁷ steps. Every quantity at or above ~10⁻⁹ survives this; the smallest
open−closed differences do not.

**Table 9.3 — GPU float32 vs float64, Co/Pt.**

| Quantity | float32 | float64 |
|---|---|---|
| Max deviation from the CPU method (E / H) | ~5×10⁻⁶ / ~2×10⁻⁶ | ~10⁻¹² expected **[re-run]** |
| R / T (`rt`) | 0.91267 / 0.08717 | 0.91267 / 0.08717 |
| Δ(m_y/M) open−closed ×10⁻⁹ (reference −2.237) | −2.057 | **[re-run]** |
| Δ(m_z/M) open−closed | ~3600× the reference (rounding floor) | resolved **[re-run: value]** |
| ΔR / ΔT open−closed (`rt`) | at the rounding floor **[re-run]** | 3.0×10⁻⁹ / −1.9×10⁻⁹ |
| Run time, one case, `full_width` (RTX 4090) | ~8 min | ~8 min |
| Run time, one case, `rt` (RTX 4090) | **[re-run]** | 5 min |

The trade-off is therefore narrow. float32 is sufficient for the
reflection/transmission ratios, the magnetization dynamics and the larger
open−closed differences, and it is the faster choice wherever arithmetic or
memory traffic limits the step. float64 is needed only for the smallest
differences. On the small grid it costs nothing measurable, because the step
is launch-bound (both precisions take ~8 min); on the large grid it doubles
the bytes per cell (Section 9.6). Since the precision is a compile-time switch,
each study can choose.

**Scaling with the magnetization.** Raising Ms artificially makes the
open−closed difference visible by eye. At the physical Ms the two M(z)
profiles coincide, and they separate as Ms grows (Figure 9.6). The growth of
|Δm_y| and |Δm_z| with the scale factor shows that the difference is the
coupling term, not numerical noise. **[re-run]**

> **[Figure 9.6]** (a) Final M_y(z) open vs closed for Ms × 1, 10, 100, 1000,
> 10⁴ [A/m vs nm]; (b) |Δm_y| and |Δm_z| vs Ms scale factor (log–log).
> *Source:* `analysis/run_ms_sweep.py` (runs and panel b),
> `analysis/compare_open_close_sweep.py` (panel a).

## 9.6 Performance

**Table 9.4 — Run times, the lab's previous CPU method vs GPU (RTX 4090,
float64).** **[re-run]**

| Benchmark | Grid cells | CPU | GPU | Speedup |
|---|---|---|---|---|
| One material case (`full_width`) | 78 | 2 h 01 m | 8 min | ~15× |
| All 12 cases (`full_width`) | 78 | 24 h 32 m | 1 h 36 m | ~15× |
| R/T, one case (`return_transmission`) | 240,889 | 15 h 02 m | 5 min | ~180× |

**Two regimes.** The speedup depends on the grid size because the bottleneck
changes:

- *78 cells — launch-bound.* One simulation fills a single block of a GPU
  with 128 multiprocessors. A step costs ~28 µs, almost all of it launch
  overhead (~95% by CUDA-event profiling). The gain over the CPU comes from
  replacing dozens of NumPy operations per step with two kernels. This is
  also why batching the 12 cases into one launch (`blockIdx.y`) is expected to
  bring the 12-case batch close to the time of a single case. **[re-run]**
- *2.4×10⁵ cells — bandwidth-bound.* Each cell update performs a handful of
  operations per value loaded, so the run is limited by memory bandwidth
  (~1 TB/s on the GPU vs ~50 GB/s on the CPU), which matches the ~180× gain.

**One curve, not two claims.** The 15× and 180× figures come from two very
different grids. To show they are two points on one curve, the grid was
widened continuously: the `full_width` material stack and cell size are kept,
and only the vacuum gap grows, from 78 cells to 10⁶. The real `rt` preset
(Δz = 1 nm) is measured as an extra point. If it lands on the curve, the
cost depends on the number of cells and not on Δz.

Figure 9.7 shows time per step against grid size. We expect a flat GPU curve
up to N ≈ 10⁵ cells, where blocks first outnumber what the 128 SMs can hold.
Below that, extra cells occupy idle SMs at no cost, and the step time is the
launch floor. Above it, time grows linearly with N at the memory-bandwidth
limit. The CPU curve is linear almost from the start. The speedup (Figure 9.8)
is therefore the ratio of a flat curve to a rising one: it climbs from the
launch-bound value at 78 cells to the bandwidth ratio around the `rt` size,
then levels off. Batching (12 cases per launch) moves the small-grid end of
the GPU curve down by up to 12×, because the 12 cases share the launch floor.
The float32 and float64 curves should coincide on the flat part and separate
on the rising part, which is the performance side of Table 9.3.
**[re-run: replace "expect" with the measured crossover N and plateau]**

> **[Figure 9.7]** Time per step and case [µs] vs grid cells N (log–log): CPU
> method, GPU float64, GPU float32, GPU float64 device-only time (CUDA
> graph, no host cost) and GPU float64 with 12 cases per launch. Dotted
> verticals at N = 78 (`full_width`) and N = 240,889 (`rt`); star = real `rt`
> preset.
> *Source:* `benchmarks/scaling_sweep.py` → `fig_scaling_time.png`.
>
> **[Figure 9.8]** Speedup over the CPU method vs grid cells N (log–log),
> same series, with the float64 speedup annotated at the two presets.
> *Source:* `benchmarks/scaling_sweep.py` → `fig_scaling_speedup.png`
> (`fig_scaling_throughput.png`, cell updates per second, is the alternative
> if one performance figure must go).

The gap between the GPU float64 curve and its device-only curve in
Figure 9.7 is the host cost per step (Python loop, two launches, frame
capture). It dominates at small N and vanishes at large N, the quantitative
form of the two regimes above.

**Profiling.** Table 9.5 confirms the regime of each preset per step. Times
come from `benchmarks/profile_modes.py`, which runs the production loop
(frames and snapshots to disk), a bare loop of the two launches, and a
CUDA-graph replay that measures pure GPU time. Nsight Systems confirms where
the host time goes, and Nsight Compute places each kernel on the roofline.
**[re-run]**

**Table 9.5 — Profile of one time step, float64, one case per launch.**

| Metric | `full_width` (78 cells) | `rt` (240,889 cells) |
|---|---|---|
| Production time per step [µs] | | |
| GPU time per step, CUDA graph [µs] (`magnetic` / `electric`) | | |
| GPU busy fraction of step time [%] | | |
| Host overhead per step [µs] (launch / frame capture and I/O) | | |
| Thread blocks / SMs (SM fill) | | |
| Registers per thread, occupancy (`magnetic_kernel`) | | |
| Effective DRAM bandwidth [GB/s] (% of peak) | | |
| Arithmetic intensity [FLOP/byte] (estimated / Nsight Compute) | | |
| VRAM in use / host RAM peak [MB] | | |
| Verdict | launch-bound (expected) | bandwidth-bound (expected) |

*All rows except the Nsight Compute intensity are columns of `report.md`
written by `profile_modes.py`. The same report covers float32 and 12 cases per
launch, for the precision and batching discussions.*

**Control experiment.** The CPU code run at float32 is only ~1.01× faster than
at float64. The speedup therefore comes from the GPU architecture
(parallelism, fewer and fused operations, bandwidth), not from lower
precision.

## 9.7 Evaluation against the objectives

| Objective | Outcome |
|---|---|
| Accelerate the coupled Maxwell–LLG simulation | 15–180× depending on grid size (Table 9.4) |
| Stay equivalent to the CPU method | Agreement at the floor of the precision used; invariants at machine precision |
| Reproduce the lab's reference results | R/T to 5 decimals; magneto-optical Δ(m_y/M) to 92% |
| Resolve the 10⁻⁹ magneto-optical regime on the GPU | Resolved; precision chosen per observable from a measured trade-off (Section 9.5) |
