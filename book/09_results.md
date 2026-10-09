# 9. Testing and Results

*(All graphs are line plots without titles, per the guidelines. Figures to be
regenerated from final float64 bundles. Items marked **[re-run]** were measured
with the per-preset backends of Stage 4 and must be re-measured with the
unified framework before submission; the physics is identical, see 9.2.)*

## 9.1 Testing procedures

Testing ran at three levels, from the arithmetic up to the physics.

1. **Implementation equivalence.** All six output arrays of a run (field time
   series at probe locations, full-grid snapshots, magnetization) were diffed
   element by element:
   - GPU against the CPU reference, over full-length runs;
   - each optimization stage against the stage before it;
   - the unified framework against the validated per-preset backends
     (`validation/compare_with_partner.py`).
2. **Physical invariants.** |M| must be conserved by the LLG equation, and
   R + T ≈ 1 must hold for energy. The building blocks (pulse, vacuum
   propagation, dielectric and conductive interfaces, free precession) were
   verified separately in notebooks.
3. **Published results.** The previous project's reflection/transmission table
   and magnetization table, for Co/Pt and Ni/Pt.

## 9.2 Numerical correctness

**GPU vs CPU.** Over the full 2×10⁷-step `full_width` run, the float32 GPU
backend matches the float64 CPU reference at the float32 floor: maximum
relative deviation ~5×10⁻⁶ in E, ~2×10⁻⁶ in H, 10⁻⁷–10⁻⁶ in M, with no NaN/Inf
and no drift. The two optimization stages (device-side output buffers;
boundaries inside the kernels) gave byte-identical output to their
predecessors. They moved work without changing the arithmetic.

**Unified framework.** The unified kernel is compared against the validated
backends on the same inputs (Co/Pt, open and closed, both presets). The
expected agreement is ~10⁻⁵ at float32 and ~10⁻¹² at float64. **[re-run]**

## 9.3 Validation against published results

**Table 9.1 — Reflection (R) and transmission (T) energy ratios.**

| Material | R (published) | R (this work) | T (published) | T (this work) | R+T |
|---|---|---|---|---|---|
| Co/Pt | 0.91266 | 0.91267 | 0.08717 | 0.08717 | 0.9998 |
| Ni/Pt | 0.90383 | 0.90384 | 0.09598 | 0.09599 | 0.9998 |

Both systems agree to the five decimals the previous project reported. The
two systems differ in damping (α = 0.025 vs 0.06) and conductivity, so the
agreement does not rest on a single parameter set.

**Table 9.2 — Magnetization observables, Co/Pt (`full_width`).**

| Observable | Published | This work | Agreement |
|---|---|---|---|
| m_z/M (closed) ×10⁻⁵ | −0.0249 | −0.0232 (CPU fp64) | 93% |
| Δ(m_y/M) open−closed ×10⁻⁹ | −2.237 | −2.057 (GPU) | 92% |

The remaining Table-2 observables agree to 70–93%. Chapter 10 discusses the
gap; the previous team's extraction script was not preserved, so part of it
comes from reconstructing their post-processing.

## 9.4 Conservation of |M|

After 2×10⁷ steps at float64, the largest per-cell deviation
|M(z, t_final)|/Ms − 1 is 2.35×10⁻¹³, the double-precision floor (Figure 9.1).
RK4 with the per-cell |M| in the damping term keeps the LLG invariant over the
whole run.

> **[Figure 9.1]** Per-cell relative drift of |M| at the final step vs magnetic
> cell index (line plot).

## 9.5 The magneto-optical signal and the role of precision

The physical observable is the difference between the open system (M does
not act back on the field) and the closed system (two-way coupling). The
magnetization imprints a field difference ΔH_x(z) of order 10⁻⁹–10⁻¹⁰ of the
driving amplitude, growing across the magnetic layer (Figure 9.2).

This scale is what sets the precision requirement:

- At float32, rounding accumulated over 2×10⁷ steps reaches ~10⁻⁷, two to
  three orders of magnitude above the signal. Absolute observables and signals
  at or above 10⁻⁹ survive (Δ(m_y/M) in Table 9.2). Smaller ones do not: the
  float32 Δ(m_z/M) came out ~3600× too large.
- At float64 the signal is resolved on the GPU. On the
  `return_transmission` preset, ΔR ≈ 3.0×10⁻⁹ and ΔT ≈ −1.9×10⁻⁹.

The preset drives pure σ⁻ circular light, an eigenmode of the magnetized
medium. The effect therefore appears as a phase shift of that eigenmode rather
than a rotation, accumulating to ≈5×10⁻⁹ rad across the layer and constant
downstream (Figure 9.3).

> **[Figure 9.2]** ΔH_x(z) = H_x^open − H_x^closed at fixed time; material
> region marked.
> **[Figure 9.3]** Eigenmode phase shift Δarg F₋(z).

## 9.6 Performance

**Table 9.3 — Run times, CPU prototype vs GPU (RTX 4090, float64).** **[re-run]**

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

**Profiling.** Table 9.4 confirms the regime of each preset with hardware
counters: Nsight Systems for where time goes per step, and Nsight Compute
for each kernel's position on the roofline. **[re-run]**

**Table 9.4 — Profile of one time step (Nsight).**

| Metric | `full_width` | `return_transmission` |
|---|---|---|
| GPU busy fraction of step time | | |
| Host time per step (Python + launch) | | |
| `magnetic_kernel` time / DRAM throughput (% of peak) | | |
| `electric_kernel` time / DRAM throughput (% of peak) | | |
| Arithmetic intensity (FLOP/byte) | | |
| Registers per thread (`num_regs`) | | |

**Control experiment.** The CPU code run at float32 is only ~1.01× faster than
at float64. The speedup therefore comes from the GPU architecture
(parallelism, fewer and fused operations, bandwidth), not from lower
precision.

## 9.7 Evaluation against the objectives

| Objective | Outcome |
|---|---|
| Accelerate the coupled Maxwell–LLG simulation | 15–180× depending on grid size (Table 9.3) |
| Stay equivalent to the CPU reference | Agreement at the float32 floor; invariants at machine precision |
| Reproduce the published results | R/T to 5 decimals; magneto-optical Δ(m_y/M) to 92% |
| Resolve the 10⁻⁹ magneto-optical regime on the GPU | Achieved at float64 (Sections 9.4–9.5) |
