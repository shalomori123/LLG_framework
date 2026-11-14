# 8. Methodology and Implementation

*(Target 4–6 pages. Figures 8.1–8.2 are the two poster diagrams, updated to the
final design. Section 8.1 follows the version history of the two development
tracks; Section 8.4 collects the architecture decisions and the evidence
behind each; Section 8.10 lists the problems met on the way.)*

## 8.1 Development history

The work ran on two parallel tracks, each in its own git repository, with one
commit per change. Track A (Shalom-Ori) built the CUDA kernels from the
bottom up, starting with the LLG step alone. Track B (Asif) ported the whole
CPU loop to CuPy first, built the analysis and validation tooling, and
produced the per-preset production backends. The tracks met in October, when
each was used to check the other (Section 8.10), and the result was merged
into one framework.

| Period | Track A — kernels (Shalom-Ori) | Track B — port, analysis, validation (Asif) |
|---|---|---|
| Dec–Jan | Start from the lab's previous CPU simulation, older versions rebased into the history. Refactor names, the LLG derivative and the pulse generator. **First bottleneck:** every step copied all field arrays (`np.copy`); swapping references instead removed most of the loop cost and gave a fair CPU baseline. | — |
| Jan | Stand-alone RK4 LLG kernel, checked against the CPU. Constants pre-computed and pre-negated on the host, so the derivative reduces to fused multiply-adds and the kernel takes fewer arguments and registers. The real part of the complex H is read directly from the complex buffer with a stride, with no host-side copy. Repeated sub-expressions hoisted. | — |
| Feb–Mar | First `maxwell_kernel.cu`. | Analysis of the CPU outputs; parameter checks with the TA (Ms scale). CuPy port of the CPU loop; all materials run together; a batched CPU version (the first form of case batching). |
| Apr | Stencil split into left/right helpers. **Two-kernel baseline** with warp-shuffle stencils, written in float32 to test whether single precision suffices. Layout: interleaved N×3 → SoA per component → **monolithic SoA** (one buffer per field, offsets computed in the kernel), which also fixed an LLG indexing error. Sullivan's absorbing boundary. **Deadlock fix:** out-of-range threads skipped barriers and shuffles; replaced by a `valid` mask. LLG folded into `magnetic_kernel` as a device function, so H^n is read once and one launch disappears. M–H coupling. M updated in place and coupled from registers; each field value loaded once into registers. `extern "C"` for CuPy name lookup; forced inlining removed, leaving it to the compiler. | First Maxwell kernels. Boundary moved to the host, then back into the kernels after a step profiler showed launch overhead dominating. Per-step |M| in the kernel. A float32 CPU version as a control experiment. Backends for the propagation-scale (`rt`) preset. Scripts that reproduce the lab's reference R/T and magnetization tables. |M|-norm check, field plots, videos, the five-panel run overview. |
| May–Jul | — | float64 GPU backends for both presets, to compare with float32. Ms sweep (supervisor suggestion). YAML configuration, repository clean-up, automatic mode selection. |
| Oct | Cross-validation against track B: seven discrepancies, each fixed in its own commit. Restructure into four files. Precision macro (one source, float32 or float64). Chunked output to disk. Case batching on `blockIdx.y`. In-kernel source with a temporary host-boundary switch to A/B-benchmark both placements, then the host variant removed. | Track B's backends are the reference for the comparison script. |
| Nov | Benchmark suite: per-mode profile and space-scaling sweep (Section 9.6). | — |

Three habits from this history shaped the final design:

- **Measure before moving work.** The boundary moved host → kernel → host →
  kernel across the two tracks. It stayed in the kernel only after a profile
  (track B) and an A/B benchmark (track A) showed what each placement cost.
- **Keep two implementations.** Neither track's own tests caught the errors
  that the cross-comparison exposed (Section 8.10).
- **Small, reversible commits.** Every optimization was a separate commit
  whose output could be diffed against its predecessor (Section 9.2).

## 8.2 Tools

| Tool | Role |
|---|---|
| Python 3, NumPy | Host orchestration, CPU method, pulse pre-computation |
| CuPy (`RawModule`) | Device arrays; runtime compilation (NVRTC) and launch of the CUDA kernels |
| CUDA C++ | `magnetic_kernel`, `electric_kernel` |
| NVIDIA Nsight Systems / Nsight Compute, CUDA graphs, cProfile | Timeline profiling (CPU vs GPU, launch gaps), per-kernel hardware metrics, pure device time, host function profile |
| RTX 4090 (lab PC), Colab T4 | Production and portability runs |
| git | Two repositories, one commit per change: each fix traceable and revertible |

## 8.3 One time step

A time step is two kernel launches (Figure 8.1). Each grid cell is handled by
one GPU thread.

**`magnetic_kernel`**
1. *Source:* the thread of cell 3 adds the pre-computed pulse value (Ex, Ey).
2. *Absorbing boundary:* the threads of cells 0 and N−1 take their values from a
   2-slot buffer, E^n(0) = E^{n−2}(1) (Section 7). The two-step delay follows
   from Δt = Δz/2c: a wave needs two steps to cross one cell.
3. *LLG (material cells):* RK4 step of M, driven by the current field H^n, with
   the per-cell |M| in the damping term.
4. *Faraday:* H update from the left E neighbour.
5. *Coupling (closed system only):* H ← H + (M_old − M_new), which keeps
   B = μ₀(H + M) continuous across the update.

**`electric_kernel`**
6. *Boundary buffer:* the threads of cells 1 and N−2 store their E^n for use two
   steps later.
7. *Ampère:* E update from the right H neighbour with per-cell coefficients
   ca, cb (Section 7). This is what makes the metal layers reflect and absorb.

The launch boundary between the two kernels is the only global
synchronization. The fields are updated **in place**: `magnetic_kernel` writes
only H and M and reads E, while `electric_kernel` writes only E and reads H. No
thread therefore reads a value another thread is writing in the same launch.
This halves the field memory compared with double buffering and removes the
pointer swap.

> **[Figure 8.1]** Per-step algorithm: the same cell program instanced over all
> cells (threads) and cases (blocks in y); the sequential time loop on the
> host.

## 8.4 Design decisions

| Decision | Options considered | Choice | Why / evidence |
|---|---|---|---|
| Where the parallelism goes | One thread per case (whole run serial); CPU processes per case; one thread per cell | One thread per cell on `blockIdx.x`, **and** one case per `blockIdx.y` | Within a half step all cells are independent. On the 78-cell grid one case fills one block of 128 SMs, so the independent cases are the second axis (track B's batched CPU version first, `blockIdx.y` in the final kernel) |
| How the step is split into kernels | One kernel per operation (pulse, ABC, Faraday, LLG, coupling, Ampère); one persistent kernel with grid-wide sync; two kernels | Two kernels: magnetic (pulse, ABC, LLG, Faraday, coupling) and electric (ABC buffer, Ampère) | The E update needs the new H of a neighbour that may live in another block, so a grid-wide barrier is unavoidable between them; the launch boundary provides it. Everything that touches the same cells and data between barriers is fused (LLG as a device function in `magnetic_kernel`) |
| Buffers | Double buffering with swap; in place | In place | Each kernel writes only fields the other reads (Section 8.3); half the memory |
| Memory layout | Interleaved N×3 (AoS); SoA per component; monolithic SoA | Monolithic SoA per field and case, complex as `real2` pairs | Coalesced access; one pointer per field with offsets in the kernel; extends directly to the case dimension |
| Neighbour values | Global memory for all; shared-memory tile; warp shuffle | Shuffle → shared at warp edges → global at block edges | Cheapest level that holds the value (Section 8.6) |
| Source and boundary placement | Host CuPy operations; in kernel with `__syncthreads()`; in kernel with owner threads | In kernel, owner-thread rule | Host placement costs ~7 extra launches per step, which dominate the small grid (35 → 8 min); owner threads make it correct for any block count (Section 8.7) |
| Host loop | Python loop per step; CUDA graph; persistent kernel | Python loop of two launches; pulse pre-computed in float64 for 65,536 steps at a time; constants pre-computed on the host | Simple and portable; the cost is the launch floor measured in Section 9.6, and CUDA graphs are the next step (Chapter 10) |
| Saving data | Python lists in RAM; device buffers copied every step; device buffers flushed in chunks | Device frame buffers flushed every K frames into memory-mapped `.npy` files of a self-describing bundle, manifest written last | RAM and VRAM independent of run length (a full run is ~4 GB per case); an interrupted run is never mistaken for a finished one |
| Arithmetic precision | float32; float64 | One kernel source, compile-time choice | An accuracy–performance trade-off that depends on the observable; both were measured (Section 8.8, 9.5) |
| Configuration and compilation | One runner script per preset and precision; YAML files; Python presets | Python presets with derived quantities in one file; kernels compiled at run time for the installed GPU | One place to change a parameter; no hard-coded GPU architecture |

## 8.5 Parallel mapping

The kernels are launched on a 2-D grid of blocks, `(blocks_x, C)`, with 256
threads per block:

- `blockIdx.x` splits the cells of one simulation across blocks.
- `blockIdx.y = c` selects the simulation case (material × open/closed). At the
  top of each kernel the case's arrays are found by offsetting the pointers
  (`E += c·2N`, `H += c·3N`, `M += c·3N_mat`, …). Per-case scalars (damping,
  open/closed) are small arrays indexed by `c`.

C = 1 is an ordinary single run; C = 12 advances all material cases in the same
launches. On the 78-cell grid one simulation occupies a single block and leaves
almost the whole GPU idle, so batching fills it at no extra launch cost.

## 8.6 Memory

**Layout.** Each field is one contiguous SoA array per case: E as [Ex | Ey], H
as [Hx | Hy | Hz], complex values as `real2` pairs, and M as [Mx | My | Mz] over
the material cells. Consecutive threads read consecutive addresses, so every
global access is coalesced. Each thread loads its field values into registers
once, updates them there, and writes them back once.

**Neighbour exchange (Figure 8.2).** Each update needs one neighbouring value.
It is fetched from the cheapest level that holds it:

- *within a warp:* register-to-register shuffle (`__shfl_up_sync` /
  `__shfl_down_sync`), with no memory traffic;
- *between warps of a block:* a small `__shared__` array of warp-edge values;
- *between blocks:* global memory, which is consistent because the other field
  is not written in this launch.

All threads, including those past the end of the grid, take part in the
shuffles and barriers; out-of-range threads are masked only for memory access.
An early `return` there was the cause of the deadlock in April.

**Output.** Fixed-size buffers on the device collect K output frames. When a
buffer fills, it is flushed into `.npy` files opened as memory maps inside the
result bundle. RAM and VRAM use is therefore independent of run length. A
full-length `full_width` run produces about 4 GB of output per case at float64,
which earlier versions held entirely in memory.

> **[Figure 8.2]** Stencil neighbour exchange: warp shuffle → shared memory at
> warp edges → global memory at block edges; ABC buffers per case.

## 8.7 Boundary and source handling inside the kernels

The source and the boundary each touch only one or two cells, but on the host
each such operation is a separate kernel launch. On the small grid, launches
are the main cost: moving six boundary operations from the host into the
kernels reduced a 20-million-step run from 35 to 8 minutes.

The first in-kernel version used `__syncthreads()`, which only synchronizes
threads of one block, so it was correct only for single-block grids. The large
grid kept a separate host-side copy. The final design avoids block
synchronization with an **owner-thread rule**: only the thread that owns a cell
reads or writes its boundary buffer.

- Cell 0's thread swaps in the boundary value as its own register value. Cell 1
  receives it by warp shuffle, which is possible because cells 0 and 1 are
  always in the same warp.
- Cell N−1 is needed by no other cell in that launch.
- Cells 1 and N−2 fill the buffers in the next kernel, after the launch
  boundary.

The same kernel is therefore correct for any number of blocks. One source
file serves both the 78-cell and the 2.4×10⁵-cell grids.

The pulse itself is computed on the CPU, in float64, for 65,536 steps at a
time. The kernel only adds one pre-computed value in one thread. This keeps the
transcendental functions off the GPU and does not increase the register count,
which is set by the heaviest path (the RK4 step).

## 8.8 Precision

Precision was treated as a parameter of the study, to be chosen by
measurement. The CPU method ran in float64. When the full time step was first
assembled on the GPU, it was written in float32 on purpose: consumer GPUs
execute float32 far faster than float64 (1/64 of the float32 rate on the
RTX 4090), and the question was whether single precision suffices for this
problem. A float32 CPU version served as a control, and float64 GPU versions
of both presets were then built to compare.

The answer depends on the observable (Section 9.5). float32 reproduces the
reflection/transmission ratios, the magnetization dynamics and the larger
open−closed differences; the smallest open−closed differences need float64.
The cost of float64 depends on the regime (Section 9.6): on the launch-bound
small grid the arithmetic is not the limit, while on the bandwidth-bound large
grid each cell moves twice the bytes. **[re-run: measured float64/float32
time ratio for both presets, from `profile_modes.py`]**

The final kernel is therefore written once against two type names, and the
precision is a compile-time choice:

```cpp
#ifdef USE_FP32
typedef float real;   typedef float2 real2;
#else
typedef double real;  typedef double2 real2;
#endif
```

The default is float64, because the open/closed studies are the main use of the
framework; `--precision fp32` compiles with `-DUSE_FP32` and fast math for
quick exploration.

## 8.9 Software structure

| File | Content |
|---|---|
| `maxwell_kernel.cu` | Device code: LLG RK4, stencils, the two kernels |
| `simulation_gpu.py` | Host: device state, coefficients, time loop, chunked output, bundle I/O |
| `sim_config.py` | Physical constants, materials table, the two presets and the formulas that derive Δt, grid indices and run length |
| `sim_run.py` | Command line: preset, materials, open/closed, precision, batching, truncated runs, benchmark |

Each run produces a self-describing bundle (`<material>_<state>.simdata/`:
`manifest.json` plus one `.npy` per array). `manifest.json` is written last, so
an interrupted run is never mistaken for a finished one. The `benchmarks/`
folder (Section 9.6) and the comparison script (Section 9.2) sit beside these
four files.

## 8.10 Challenges and solutions

**Performance and correctness of the GPU design**

| Challenge | Solution |
|---|---|
| CPU loop copied all field arrays every step | Reference swapping on the CPU; in-place updates on the GPU (Section 8.3) |
| Kernel hung when N was not a multiple of the block size: out-of-range threads returned early and skipped `__syncthreads()` / shuffles | All threads take part in barriers and shuffles; out-of-range threads are masked (`valid` flag) only for memory access |
| Neighbour reads of the stencil | Warp shuffle → shared memory → global memory hierarchy (Section 8.6) |
| Repeated global reads of H and M inside a step | LLG fused into `magnetic_kernel`; each value loaded once into registers, coupling applied from registers |
| Small grid is launch-bound: profiling attributed ~95% of step time to launch overhead | Source and boundary moved into the kernels (35 → 8 min); batching via `blockIdx.y` fills the idle GPU |
| `__syncthreads()` cannot order work across blocks, so in-kernel boundaries only worked on one block | Owner-thread rule; one kernel for all grid sizes (Section 8.7) |
| Output time series grew with run length and exhausted RAM | Fixed device buffers flushed to memory-mapped files (Section 8.6) |
| Six runners and four near-identical kernels (precision × grid size) | One kernel source, four files, Python presets |

**Discrepancies found by cross-validating the two development tracks**

Running both implementations on identical inputs turned up the following
differences. Each was resolved against Sullivan's FDTD text or the CPU
method, and fixed in its own commit:

| Discrepancy | Effect | Resolution |
|---|---|---|
| Open/closed switch was a compile-time `#define`, always on | Every "open" run was actually closed | Runtime per-case flag |
| Out-of-range threads read and wrote past the field arrays | Corrupted Ey near the array start | Guarded loads and stores |
| One sign in each E/H update pair was flipped | Not a wave equation: pulses grow instead of propagating | Signs from Sullivan Ch. 1; stencil direction matched to the CPU (H uses left E, E uses right H) |
| Constant 0.5 coefficient in Ampère's law | Conductors behaved like vacuum | Per-cell ca/cb (Section 7) |
| LLG driven by H^{n+1} instead of H^n | Half-step shift relative to the CPU method | LLG evaluated before the Faraday update |
| Fixed |M| in the damping term | Slow drift of |M| over long runs | |M| computed per cell from M at the start of each step |
| 1-slot boundary buffer | Partial reflection at the boundary for Δt = Δz/2c | 2-slot buffer, E^n(0) = E^{n−2}(1) |
| Initial magnetization 10⁶× too large | Divergence at the pulse peak through the coupling term | Physical Ms = 3×10⁵ A/m |
