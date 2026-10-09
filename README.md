# LLG_framework

1D FDTD Maxwell solver (Yee grid, Sullivan convention) coupled to the
Landau-Lifshitz-Gilbert equation (RK4), for a magnetic layer + conductive cap
excited by a circularly polarized laser pulse.

## GPU simulation (`gpu_tools/`)

| File | Role |
| --- | --- |
| `maxwell_kernel.cu` | Device code: `magnetic_kernel` (pulse, ABC, Faraday, LLG, coupling) and `electric_kernel` (Ampere, ABC buffer) |
| `simulation_gpu.py` | Host: memory setup, time loop, chunked saving into `.simdata` bundles |
| `sim_config.py` | Physical constants, materials table, presets (`full_width`, `rt`) |
| `sim_run.py` | Command line entry point |

Requires CuPy (`pip install cupy-cuda12x`), NumPy and tqdm.

```
cd gpu_tools
python sim_run.py                                     # full_width, all 6 materials, open + close
python sim_run.py --batch                             # same, as many cases per launch as fit in VRAM
python sim_run.py --preset rt --materials Co_Pt       # return & transmission preset
python sim_run.py --steps 30000 --output-dir quick    # truncated run for quick checks
python sim_run.py --precision fp32                    # faster, but can't resolve open vs close
python sim_run.py --benchmark --materials Co_Pt       # steps/s and kernel registers
```

Results go to `simulation_result/<preset>/<material>_<state>.simdata/`
(`manifest.json` + one `.npy` per array). Load with
`simulation_gpu.load_simulation_data(path)`; arrays are memory-mapped.

### Design notes

- **Precision**: `real`/`real2` typedefs, fp64 by default, `-DUSE_FP32` with `--precision fp32`.
- **Batching**: launch grid `(blocks_x, C)`. `blockIdx.x` splits the cells,
  `blockIdx.y` picks the case; C = 1 is a sequential run.
- **Boundaries**: Mur ABC with a 2-step delay (Sullivan eq. 1.12, Courant 0.5)
  and the soft source run inside the kernels; only the thread owning a cell
  touches it, so any number of blocks is safe. The pulse is computed on the CPU
  in fp64, one chunk of steps at a time.
- **Memory**: per-location frames are buffered on the GPU and flushed every K
  frames into memmap files, so RAM and VRAM don't grow with the run length.

## Validation

`validation/compare_with_partner.py` runs the same truncated Co_Pt cases
through this code and the partner repo's GPU backends and prints the max
relative difference per output array.

## Other folders

`python_core/` (CPU reference and verification notebooks) and
`matlab_attempts/` (early MATLAB prototypes).
