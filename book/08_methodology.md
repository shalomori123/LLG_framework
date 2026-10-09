# 8. Methodology and Implementation

*(Target 4–6 pages. Figures 8.1–8.2 are the two poster diagrams, updated to the
final design. Every design choice below is paired with the problem that forced
it; Section 8.9 collects them.)*

## 8.1 Development path

The project went through five stages. Each stage was driven by a specific
limitation of the previous one.

1. **Inherited CPU reference (Dec–Jan).** We started from the previous team's
   NumPy implementation (float64) and refactored it for readability. The first
   bottleneck appeared here: every time step copied the full field arrays
   (`np.copy`). Swapping buffer references instead removed most of the loop
   cost and gave us a fair CPU baseline.
2. **Isolated LLG kernel (Jan).** The RK4 LLG step was ported to CUDA alone
   and checked against the CPU. Here we learned to precompute constants on the
   host (so the device code reduces to fused multiply-adds) and to remove
   repeated sub-expressions.
3. **Full two-kernel time step (Feb–Apr).** The Maxwell updates joined the LLG
   step in two kernels. This stage set the architecture: one thread per cell,
   structure-of-arrays (SoA) memory, warp-shuffle stencils, the absorbing
   boundary and the M→H coupling.
4. **Validation and performance campaign (Apr–Sep).** A parallel track ran the
   GPU code against the CPU reference, reproduced the previous team's published
   tables, moved the boundary handling into the kernels (35 → 8 min per run),
   profiled the step, and promoted the kernels to float64 to resolve the 10⁻⁹
   magneto-optical signal.
5. **Unification (Oct).** Cross-validating the two development tracks exposed
   seven discrepancies (Section 8.9). After fixing them, the code was merged
   into a single kernel source that covers both precisions, both grid sizes and
   batched runs, plus three Python files.

## 8.2 Tools

| Tool | Role |
|---|---|
| Python 3, NumPy | Host orchestration, CPU reference, pulse precomputation |
| CuPy (`RawModule`) | Device arrays; runtime compilation (NVRTC) and launch of the CUDA kernels |
| CUDA C++ | `magnetic_kernel`, `electric_kernel` |
| NVIDIA Nsight Systems / Nsight Compute | Timeline profiling (CPU vs GPU, launch gaps) and per-kernel hardware metrics (roofline, stall reasons) |
| RTX 4090 (lab PC), Colab T4 | Production and portability runs |
| git | Version control: one commit per change, which made each fix traceable and revertible |

## 8.3 One time step

A time step is two kernel launches (Figure 8.1). Each grid cell is handled by
one GPU thread.

**`magnetic_kernel`**
1. *Source:* the thread of cell 3 adds the precomputed pulse value (Ex, Ey).
2. *Absorbing boundary:* the threads of cells 0 and N−1 take their values from a
   2-slot buffer, E^n(0) = E^{n−2}(1) (Mur, Sullivan Eq. 1.12). The two-step
   delay follows from Δt = Δz/2c: a wave needs two steps to cross one cell.
3. *LLG (material cells):* RK4 step of M, driven by the current field H^n, with
   the per-cell |M| in the damping term.
4. *Faraday:* H update from the left E neighbour.
5. *Coupling (closed system only):* H ← H + (M_old − M_new), which keeps
   B = μ₀(H + M) continuous across the update.

**`electric_kernel`**
6. *Boundary buffer:* the threads of cells 1 and N−2 store their E^n for use two
   steps later.
7. *Ampère:* E update from the right H neighbour with per-cell coefficients
   ca = (1−η)/(1+η), cb = 0.5/(ε_r(1+η)), η = σΔt/(2ε₀ε_r) (Sullivan Eq.
   1.23–1.24). This is what makes the metal layers reflect and absorb.

The launch boundary between the two kernels is the only global
synchronization. The fields are updated **in place**: `magnetic_kernel` writes
only H and M and reads E, while `electric_kernel` writes only E and reads H. No
thread therefore reads a value another thread is writing in the same launch.
This halves the field memory compared with double buffering and removes the
pointer swap.

> **[Figure 8.1]** Per-step algorithm: the same cell program instanced over all
> cells (threads) and cases (blocks in y); the sequential time loop on the
> host.

## 8.4 Parallel mapping

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

## 8.5 Memory

**Layout.** Each field is one contiguous SoA array per case: E as [Ex | Ey], H
as [Hx | Hy | Hz], complex values as `real2` pairs, and M as [Mx | My | Mz] over
the material cells. Consecutive threads read consecutive addresses, so every
global access is coalesced.

**Neighbour exchange (Figure 8.2).** Each update needs one neighbouring value.
It is fetched from the cheapest level that holds it:

- *within a warp:* register-to-register shuffle (`__shfl_up_sync` /
  `__shfl_down_sync`), with no memory traffic;
- *between warps of a block:* a small `__shared__` array of warp-edge values;
- *between blocks:* global memory, which is consistent because the other field
  is not written in this launch.

**Output.** Fixed-size buffers on the device collect K output frames. When a
buffer fills, it is flushed into `.npy` files opened as memory maps inside the
result bundle. RAM and VRAM use is therefore independent of run length. A
full-length `full_width` run produces about 4 GB of output per case at float64,
which earlier versions held entirely in memory.

> **[Figure 8.2]** Stencil neighbour exchange: warp shuffle → shared memory at
> warp edges → global memory at block edges; ABC buffers per case.

## 8.6 Boundary and source handling inside the kernels

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
time. The kernel only adds one precomputed value in one thread. This keeps the
transcendental functions off the GPU and does not increase the register count,
which is set by the heaviest path (the RK4 step).

## 8.7 Precision

The kernel is written once against two type names:

```cpp
#ifdef USE_FP32
typedef float real;   typedef float2 real2;
#else
typedef double real;  typedef double2 real2;
#endif
```

float64 is the default because the magneto-optical signal is ~10⁻⁹ of the
field (Section 9.5). `--precision fp32` compiles with `-DUSE_FP32` and fast
math, for quick exploration. No GPU architecture is hard-coded: the kernels are
compiled at run time for the installed GPU.

## 8.8 Software structure

| File | Content |
|---|---|
| `maxwell_kernel.cu` | Device code: LLG RK4, stencils, the two kernels |
| `simulation_gpu.py` | Host: allocation, coefficients, time loop, chunked output, bundle I/O |
| `sim_config.py` | Physical constants, materials table, the two presets and the formulas that derive Δt, grid indices and run length |
| `sim_run.py` | Command line: preset, materials, open/closed, precision, batching, truncated runs, benchmark |

Each run produces a self-describing bundle (`<material>_<state>.simdata/`:
`manifest.json` plus one `.npy` per array). `manifest.json` is written last, so
an interrupted run is never mistaken for a finished one.

## 8.9 Challenges and solutions

**Performance and correctness of the GPU design**

| Challenge | Solution |
|---|---|
| CPU loop copied all field arrays every step | Buffer swapping on the CPU; in-place updates on the GPU (Section 8.3) |
| Kernel hung when N was not a multiple of the block size: out-of-range threads returned early and skipped `__syncthreads()` / shuffles | All threads take part in barriers and shuffles; out-of-range threads are masked (`valid` flag) only for memory access |
| Neighbour reads of the stencil | Warp shuffle → shared memory → global memory hierarchy (Section 8.5) |
| Small grid is launch-bound: profiling attributed ~95% of step time to launch overhead | Source and boundary moved into the kernels (35 → 8 min); batching via `blockIdx.y` fills the idle GPU |
| `__syncthreads()` cannot order work across blocks, so in-kernel boundaries only worked on one block | Owner-thread rule; one kernel for all grid sizes (Section 8.6) |
| float32 rounding over 2×10⁷ steps (~10⁻⁷) buries the 10⁻⁹ signal | float64 default, float32 kept as a compile-time option (Section 8.7) |
| Output time series grew with run length and exhausted RAM | Fixed device buffers flushed to memory-mapped files (Section 8.5) |
| Six runners and four near-identical kernels (precision × grid size) | One kernel source, four files, Python presets |

**Discrepancies found by cross-validating the two development tracks**

Running both implementations on identical inputs turned up the following
differences. Each was resolved against Sullivan's FDTD text or the CPU
reference, and fixed in its own commit:

| Discrepancy | Effect | Resolution |
|---|---|---|
| Open/closed switch was a compile-time `#define`, always on | Every "open" run was actually closed | Runtime per-case flag |
| Out-of-range threads read and wrote past the field arrays | Corrupted Ey near the array start | Guarded loads and stores |
| One sign in each Yee pair was flipped | Not a wave equation: pulses grow instead of propagating | Signs from Sullivan Eq. 1.9; stencil direction matched to the CPU (H uses left E, E uses right H) |
| Constant 0.5 coefficient in Ampère's law | Conductors behaved like vacuum | Per-cell ca/cb (Sullivan Eq. 1.23–1.24) |
| LLG driven by H^{n+1} instead of H^n | Half-step shift vs the CPU reference and published tables | LLG evaluated before the Faraday update |
| Fixed |M| in the damping term | Slow drift of |M| over long runs | |M| computed per cell from M at the start of each step |
| 1-slot boundary buffer | Partial reflection at the boundary for Δt = Δz/2c | 2-slot buffer, E^n(0) = E^{n−2}(1) |
| Initial magnetization 10⁶× too large | Divergence at the pulse peak through the coupling term | Physical Ms = 3×10⁵ A/m |
