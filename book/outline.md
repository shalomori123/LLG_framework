# Project Book — Outline
### "Parallel GPU Computing of Light–Magnetism Interaction"
Asif Levi & Shalom-Ori Danino · Supervisors: Prof. Amir Capua, Dr. Benjamin Assouline

Structure maps 1:1 to `docs/Book_guidelines_Nirel_v2.pdf` (points in parentheses).
Grading: 90 pts for sections + 10 pts readability. **Results = 20 pts, double
any other section — write it first and best.**

Global figure rules (from the guidelines — apply everywhere):
- **Line plots only** for 2-D graphs (no bar/pie). Scalar comparisons → tables.
- **No titles on graphs** — numbered captions below, referenced from the text
  ("Figure 3 shows…").
- Proper x/y axis labels with units; legends where needed.

---

## 1. Title Page
Title, both names + IDs + emails, advisors + affiliation (HUJI, Applied
Physics / Electrical Engineering), submission date.

## 2. Abstract (10 pts, 150–300 words)
**Status: DONE** — use the drafted abstract ("Parallel GPU Computing of
Light–Magnetism Interaction", ~250 words), already incorporating advisor
feedback (no "weak effect" framing; generic "published results" instead of
naming Choi).

## 3. Table of Contents
Auto-generated at the end.

## 4. List of Abbreviations and Symbols
FDTD, LLG, GPU, CUDA, CPU, RK4, ABC (absorbing boundary condition), fp32/fp64,
R/T, MO (magneto-optical), SoA; symbols: E, H, M, Ms, α (Gilbert damping),
γ (gyromagnetic ratio), λ, σ (conductivity), ε_r, Δz, Δt, θ_F.

## 5. Introduction (10 pts — one paragraph each)
- **Background:** light–matter interaction; the optical magnetic field usually
  neglected; advisors' research + new findings show it matters (large effects
  for circular polarization). *Source: poster Motivation + abstract ¶1.*
- **Objectives:** build an accelerated framework for coupled Maxwell–LLG
  simulation, resolving the ~10⁻⁹ MO signal (fp64) at practical runtimes.
- **Scope:** 1-D geometry; Co/Fe/Ni × Au/Pt thin films; two presets
  (subwavelength full_width, propagation-scale return_transmission); single-GPU
  CUDA; validation against the previous project's published tables.
- **Organization:** one sentence per chapter.

## 6. Literature and state-of-the-art review (10 pts — 1–2 paragraphs)
- FDTD (Yee) as standard EM solver; LLG for magnetization dynamics; the coupled
  Maxwell–LLG scheme (Choi et al. 2017 — cite here even though the abstract
  stays generic).
- Advisors' publications on optical-magnetic-field effects (get exact refs
  from Amir/Benjamin).
- Previous project (Or & Reem): CPU implementation + published R/T and
  magnetization tables — the baseline and validation target.
- GPU-accelerated FDTD literature: memory-bandwidth-bound stencils, one thread
  per cell. **Gap:** no available GPU implementation of the *coupled* fp64
  Maxwell–LLG system → this project.

## 7. Theoretical Background (10 pts — ~1 page)
- Maxwell curl equations → 1-D Yee staggering of E and H.
- LLG equation (precession + Gilbert damping term), conserves |M|.
- Two-way coupling: optical H drives M (LLG); M feeds back into Faraday's law.
- The MO observable: open-vs-close difference at ~10⁻⁹ of the field → motivates
  fp64. Circular polarization = eigenmodes of the magnetized medium (rotation
  vs eigenmode-phase-shift distinction — feeds Results §9 Faraday part).
- Numerical stability: Courant condition, dt = dz/2c; RK4 for the stiff LLG.
*Source: poster physics section + script §3; keep to 1 page per guidelines.*

## 8. Methodology and Implementation (10 pts — 4–6 pages)
**Status: DRAFTED** (`08_methodology.md`). Every design choice is tied to the
problem that forced it.
- **Development path:** 5 stages from git history: CPU reference → isolated
  LLG kernel → two-kernel step → validation/performance/fp64 campaign →
  cross-validation and unification.
- **Algorithm:** `magnetic_kernel` (pulse, ABC, LLG on H^n, Faraday, coupling)
  → `electric_kernel` (ABC buffer, Ampère with ca/cb). In-place updates, and
  why that is safe.
- **Parallel mapping:** one thread per cell; grid `(blocks_x, C)` with
  `blockIdx.y` = case (batching).
- **Memory:** monolithic SoA, warp shuffle → shared → global, chunked
  memmap output.
- **Boundaries in-kernel:** launch-count argument (35 → 8 min) and the
  owner-thread rule (one kernel for every grid size).
- **Precision macro**, **4-file software structure**.
- **Challenges & solutions:** two tables, design challenges and the 8
  discrepancies found by cross-validation.
- **Block diagrams:** update the two poster diagrams: add the `blockIdx.y`
  case dimension, and replace the double-buffer swap with in-place updates.

## 9. Testing and Results (20 pts — 4–6 pages, MOST IMPORTANT)
Order: correctness → physics → performance.
- **9.1 Testing procedures:** CPU-vs-GPU bit-level diffs over full runs;
  conservation checks; reproduction of previous team's published tables;
  verification notebooks (pulse, dielectric, conductivity, LLG).
- **9.2 Numerical correctness:** CPU↔GPU agreement (~10⁻⁶ at fp32 floor,
  RESULTS.md §3); optimization stages byte-identical.
- **9.3 Validation vs published results:**
  - **Table 1 (book):** R/T vs previous team — Co_Pt 0.91267/0.91266,
    0.08717/0.08717; Ni_Pt likewise (5 decimals). TABLE, not a plot.
  - **Table 2 (book):** magnetization ratios mz/M, Δ(my/M) vs published values
    (93% / 92%). Include the honest gap discussion pointer to §10.
- **9.4 |M| conservation:** line plot of per-cell relative drift (~10⁻¹³),
  final frame. *(Remake as LINE plot from the existing per-cell data — the
  two-bar poster figure is not guideline-compliant.)*
- **9.5 The magneto-optical signal (the physics payoff):** open−close field
  difference ΔH(z) at ~10⁻⁹ (line plot); resolved on GPU **only at fp64** —
  fp32 comparison showing the floor (numbers from RESULTS.md §4). Faraday
  eigenmode phase shift Δarg F₋(z) ≈ 5×10⁻⁹ rad (line plot from
  `analysis/faraday_rotation.py`) + the σ± eigenmode explanation.
- **9.6 Performance:** runtime table CPU vs GPU per preset; launch-bound vs
  bandwidth-bound regimes; batching; Nsight profile table (9.4); CPU-fp32
  control experiment (~1.01×).
- **9.7 Performance evaluation vs objectives:** goals met per Introduction.

## 10. Discussion (10 pts — 2–3 paragraphs)
- Interpretation: GPU fp64 makes the previously CPU-only precision regime
  practical; speedup regime depends on problem size (launch-bound vs
  bandwidth-bound).
- Discrepancies (required): Table-2 residual gap (methodology drift, original
  extraction script not preserved); R↔T Δ-label swap caveat; rt CPU full-run
  never executed (16 hr) — truncated 30k-step validation only.
- Improvements/future work: kernel fusion / CUDA Graphs; remaining 10 material
  combos; 2-D/3-D extension; HDF5 saving; GUI.

## 11. Conclusions (10 pts — one paragraph)
Coupled Maxwell–LLG framework on GPU at fp64; validated to 5 decimals against
published results; resolves 10⁻⁹ MO signatures; 15–180× speedup (pending final
fp64 numbers); reflection on what was learned (GPU architecture ↔ physics
precision interplay).

## 13. References (IEEE style)
Choi et al. 2017; advisors' papers (ask for exact list); Yee 1966; LLG
originals (Landau–Lifshitz 1935, Gilbert 2004); Or & Reem project book; CUDA
programming guide; CuPy.

## 14. Appendices
- `sim_config.py` materials table + presets (full_width + rt) — shows reproducibility.
- `.simdata` bundle format.
- Repository guide (README summary) + how to reproduce each figure/table
  (script + command per item).

---

# Open items the book is waiting on
| # | Item | Needed for | Effort |
|---|---|---|---|
| 1 | ~~fp64 runtime measurements~~ **DONE** — taken from the poster benchmark chart (15×/15×/180×, RESULTS.md §9.2) | §9.6 | — |
| 2 | R↔T Δ-label convention check vs book | §9.5, §10 | ~1 hr |
| 3 | Linear-polarization run (`phase_diff_over_pi: 0`) for true θ_F(z) | §9.5 (optional, nice) | 1 Colab run |
| 4 | Advisors' exact paper references | §6, §13 | email |
| 5 | Remake conservation figure as line plot | §9.4 | ½ hr |
| 6 | Trusted fp64 bundles for final figures (lab PC), produced by the unified framework | §9.4–9.5 | reuse #7 session |
| 7 | Lab session with the unified framework: `validation/compare_with_partner.py` (fp32 and fp64), `benchmarks/profile_modes.py` (all modes), re-time Table 9.3 | §9.2, §9.6 (all **[re-run]** marks) | ½ day |
| 8 | `benchmarks/profile_modes.py --presets full_width rt --precisions fp64 --batches 1 --nsys --ncu`, fill Table 9.4 | §9.6 | 1–2 hr |
| 10 | `benchmarks/scaling_sweep.py` → Figures 9.7–9.8; replace "expect" in §9.6 with the measured crossover and plateau | §9.6 | ~1 hr GPU (CPU points dominate) |
| 11 | Book versions of the analysis figures (no titles, bars → lines): `advisor_view_compare.py` panels → Figures 9.1 + 9.3, `make_poster_mo_signal.py` → Figure 9.4 | §9.2–9.5 | 1–2 hr |
| 9 | Update poster diagrams (batch dimension, in-place updates) | §8 Figures 8.1–8.2 | 1 hr |
