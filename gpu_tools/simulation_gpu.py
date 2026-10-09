import json
import os
import shutil
from functools import lru_cache
from pathlib import Path

import numpy as np
import cupy as cp
from tqdm import tqdm

from sim_config import EPSILON_0, MU_0, SPEED_OF_LIGHT, GAMMA, grid_size as compute_grid_size

# Host dtypes matching the kernel's real / real2 typedefs
DTYPES = {
    "fp64": (np.float64, np.complex128),
    "fp32": (np.float32, np.complex64),
}

# Device memory for the per-location frame buffer, flushed to the output when full
FRAME_BUFFER_BYTES = 256 * 2**20
THREADS_PER_BLOCK = 256
# Steps of source pulse computed on the CPU (fp64, vectorized) and uploaded at once
PULSE_CHUNK_STEPS = 65536


def pulse_gen(max_E: float, phase: float, FWHM: float, dt: float, 
              wave_len=800e-9, angle=0, safety_start=4):
    """Generates the CPU function for the complex E-field components."""
    omega = 2 * np.pi * SPEED_OF_LIGHT / wave_len
    angle -= np.pi/4
    
    R = np.array([
        [np.cos(angle), -np.sin(angle)], 
        [np.sin(angle), np.cos(angle)]
    ])
    
    A_vec = np.array([1, np.exp(-phase * 1j)]).reshape(2) 
    pulse_amp = max_E * R @ A_vec 
    sigma = (0.5 / (np.log(2)**0.5)) * FWHM
    peak_loc = safety_start * sigma

    def out_pulse(n) -> np.array:
        """n: step index or array of K indices -> (2,) or (K, 2) complex [Ex, Ey]."""
        t = np.asarray(n) * dt
        envelope = np.exp(- (t - peak_loc)**2 / (2 * sigma**2)) * np.exp(1j * omega * t)
        return envelope[..., None] * pulse_amp
        
    return out_pulse


@lru_cache(maxsize=None)
def load_kernels(precision: str = "fp64"):
    """Compile maxwell_kernel.cu once per process and precision.
    No -arch: NVRTC targets the installed GPU."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, "maxwell_kernel.cu"), "r") as f:
        cuda_code = f.read()
    options = ('-DUSE_FP32', '-use_fast_math') if precision == "fp32" else ()
    module = cp.RawModule(code=cuda_code, options=options)
    return module.get_function("magnetic_kernel"), module.get_function("electric_kernel")


def case_device_bytes(grid_size: int, material_size: int, precision: str = "fp64") -> int:
    """VRAM of one case's state: E, H, coefficients, M, ABC buffers."""
    real_t, complex_t = DTYPES[precision]
    r, c = np.dtype(real_t).itemsize, np.dtype(complex_t).itemsize
    return 5 * grid_size * c + 2 * grid_size * r + 3 * material_size * r + 8 * c


def max_batch(grid_size: int, material_size: int, precision: str = "fp64") -> int:
    """How many cases fit in the free VRAM at once (frame buffer is shared)."""
    free_bytes, _ = cp.cuda.Device().mem_info
    usable = 0.9 * free_bytes - FRAME_BUFFER_BYTES
    return max(1, int(usable // case_device_bytes(grid_size, material_size, precision)))


def simulation_cases(cases: list[dict], z_indices: list[int], M0: np.array,
                     max_E: float, phase_diff: float, pulse_width: float,
                     save_locations, save_time_indices,
                     dt: float, total_time_steps: int, wavelength=800e-9, polarization_angle=0, safety_factor=4,
                     time_interval=1, precision="fp64", output_paths=None):
    """
    Run C cases that share the grid and pulse in the same kernel launches
    (blockIdx.y = case). Each case is a dict with eps_r, conductivity, damping,
    closed_system. C = 1 is a plain sequential run.

    With output_paths=None the results are returned in RAM; otherwise case c
    streams into <output_paths[c]>.simdata (memmap .npy files), so neither RAM
    nor VRAM grows with the number of time steps.

    Pulse injection and ABC run inside magnetic_kernel/electric_kernel: as host
    CuPy ops they cost ~7 extra launches per step, which dominate small grids.
    """
    real_t, complex_t = DTYPES[precision]
    C = len(cases)
    save_locations = np.asarray(save_locations, dtype=np.int64)
    save_time_indices = np.asarray(save_time_indices, dtype=np.int64)

    # --- 1. Grid Definition ---
    material_begin = z_indices[0]
    material_end   = z_indices[-1]
    material_size  = material_end - material_begin
    grid_size = compute_grid_size(z_indices, dt, wavelength)

    # --- 2. LLG Physics Parameters (per case) ---
    # The damping term is divided by the per-cell |M|, computed in the kernel each step
    damping = np.array([case["damping"] for case in cases], dtype=np.float64)
    gamma_ll = GAMMA / (1 + damping**2)
    lambda_ll = gamma_ll * damping
    neg_gamma_LL_gpu = cp.asarray((-gamma_ll * MU_0).astype(real_t))
    neg_lambda_mu_gpu = cp.asarray((-lambda_ll * MU_0).astype(real_t))
    closed_system_gpu = cp.asarray([int(case["closed_system"]) for case in cases], dtype=cp.int32)

    # --- 3. GPU Memory Allocation (Monolithic SoA, one row per case) ---
    # Complex arrays for E and H (real2 in CUDA)
    E_gpu = cp.zeros((C, 2 * grid_size), dtype=complex_t)
    H_gpu = cp.zeros((C, 3 * grid_size), dtype=complex_t)

    # ABC Buffers - 2-step delay at Courant 0.5: [2 slots][Ex, Ey]
    abc_left_gpu = cp.zeros((C, 4), dtype=complex_t)
    abc_right_gpu = cp.zeros((C, 4), dtype=complex_t)

    # --- 4. Initialize M on GPU ---
    M_cpu_init = np.zeros((material_size, 3), dtype=real_t)
    for i in range(len(z_indices) - 1):
        start_idx = z_indices[i] - material_begin
        end_idx = z_indices[i+1] - material_begin
        for k in range(3):
            M_cpu_init[start_idx:end_idx, k] = M0[i, k] / (end_idx - start_idx)

    # Convert AoS to Monolithic SoA and upload to GPU
    M_soa_cpu = np.concatenate([M_cpu_init[:, 0], M_cpu_init[:, 1], M_cpu_init[:, 2]])
    M_gpu = cp.asarray(np.tile(M_soa_cpu, (C, 1)))

    # --- 5. Injection Pulse ---
    # Computed on the CPU in fp64, one chunk of steps at a time (refreshed in the loop)
    pulse_func = pulse_gen(max_E, phase_diff, pulse_width, dt, wavelength, polarization_angle, safety_factor)
    pulse_chunk_gpu = None

    # --- 6. Load CUDA Kernels ---
    magnetic_kernel, electric_kernel = load_kernels(precision)

    # Block and Grid dimensions: x splits the cells, y selects the case
    blocks_per_grid = ((grid_size + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK, C)

    # Faraday coefficient (vacuum, Courant 0.5)
    coeff = real_t(0.5)

    # Per-cell Ampere coefficients for lossy dielectric [Sullivan eq. 3.8]
    coeff_E_cpu = np.ones((C, grid_size), dtype=real_t)
    coeff_H_cpu = np.full((C, grid_size), 0.5, dtype=real_t)
    for c, case in enumerate(cases):
        for i in range(len(case["eps_r"])):
            eaf = dt * case["conductivity"][i] / (2 * EPSILON_0 * case["eps_r"][i])
            coeff_E_cpu[c, z_indices[i]:z_indices[i+1]] = (1 - eaf) / (1 + eaf)
            coeff_H_cpu[c, z_indices[i]:z_indices[i+1]] = 0.5 / (case["eps_r"][i] * (1 + eaf))
    coeff_E_gpu = cp.asarray(coeff_E_cpu)
    coeff_H_gpu = cp.asarray(coeff_H_cpu)

    # --- 7. Outputs per case (RAM arrays, or memmap files inside the bundle) ---
    total_saved_frames = int(total_time_steps / time_interval + 1)
    n_loc = len(save_locations)
    n_snap = len(save_time_indices)
    shapes = {
        "E in locations": ((total_saved_frames, n_loc, 2), complex_t),
        "H in locations": ((total_saved_frames, n_loc, 2), complex_t),
        "ms in locations": ((total_saved_frames, material_size, 3), real_t),
        "E times picture": ((n_snap, grid_size, 2), complex_t),
        "H times picture": ((n_snap, grid_size, 2), complex_t),
        "ms time picture": ((n_snap, material_size, 3), real_t),
    }
    writers = [BundleWriter(path) for path in output_paths] if output_paths is not None else [None] * C
    outs = [{key: (w.array(key, shape, dtype) if w else np.zeros(shape, dtype))
             for key, (shape, dtype) in shapes.items()} for w in writers]

    # Device frame buffer: K frames x C cases of [x-comp locations, y-comp locations] and SoA M
    E_loc_index = cp.asarray(np.concatenate([save_locations, grid_size + save_locations]))
    frame_bytes = C * (4 * n_loc * np.dtype(complex_t).itemsize + 3 * material_size * np.dtype(real_t).itemsize)
    chunk_frames = int(max(1, min(total_saved_frames, FRAME_BUFFER_BYTES // frame_bytes)))
    E_loc_buf = cp.zeros((chunk_frames, C, 2 * n_loc), dtype=complex_t)
    H_loc_buf = cp.zeros((chunk_frames, C, 2 * n_loc), dtype=complex_t)
    M_buf = cp.zeros((chunk_frames, C, 3 * material_size), dtype=real_t)
    buffered, flushed = 0, 0

    def flush():
        nonlocal buffered, flushed
        if buffered == 0:
            return
        sl = slice(flushed, flushed + buffered)
        E_host, H_host, M_host = E_loc_buf[:buffered].get(), H_loc_buf[:buffered].get(), M_buf[:buffered].get()
        for c, out in enumerate(outs):
            out["E in locations"][sl] = E_host[:, c].reshape(buffered, 2, n_loc).transpose(0, 2, 1)
            out["H in locations"][sl] = H_host[:, c].reshape(buffered, 2, n_loc).transpose(0, 2, 1)
            out["ms in locations"][sl] = M_host[:, c].reshape(buffered, 3, material_size).transpose(0, 2, 1)
        flushed += buffered
        buffered = 0

    save_time_set = set(save_time_indices.tolist())
    save_index_counter = 0

    # --- 8. THE MAIN FDTD LOOP ---
    for time_idx in tqdm(range(0, total_time_steps), desc=f"GPU Simulation Process ({C} case{'s' * (C > 1)}): "):
        
        # 1. Source pulse for this step: [Ex, Ey] row of the current chunk
        step_in_chunk = time_idx % PULSE_CHUNK_STEPS
        if step_in_chunk == 0:
            steps = np.arange(time_idx, min(time_idx + PULSE_CHUNK_STEPS, total_time_steps))
            pulse_chunk_gpu = cp.asarray(pulse_func(steps).astype(complex_t))
        pulse_gpu = pulse_chunk_gpu[step_in_chunk]

        # 2. Magnetic Kernel (Injects pulse, applies ABC, updates H, steps LLG, applies Coupling)
        magnetic_kernel(
            blocks_per_grid, (THREADS_PER_BLOCK,),
            (E_gpu, H_gpu, M_gpu, 
             np.int32(grid_size), coeff, 
             np.int32(material_begin), np.int32(material_end), 
             real_t(dt), neg_gamma_LL_gpu, neg_lambda_mu_gpu,
             closed_system_gpu, abc_left_gpu, abc_right_gpu, pulse_gpu)
        )

        # 3. Electric Kernel (Updates E, fills ABC buffers)
        electric_kernel(
            blocks_per_grid, (THREADS_PER_BLOCK,),
            (E_gpu, H_gpu, np.int32(grid_size), coeff_E_gpu, coeff_H_gpu, abc_left_gpu, abc_right_gpu)
        )

        # --- 9. Data Extraction ---
        # Per-location frames stay on the GPU until the buffer is full
        if not time_idx % time_interval:
            cp.take(E_gpu, E_loc_index, axis=1, out=E_loc_buf[buffered])
            cp.take(H_gpu, E_loc_index, axis=1, out=H_loc_buf[buffered])
            M_buf[buffered] = M_gpu
            buffered += 1
            if buffered == chunk_frames:
                flush()

        # Full-grid snapshots are rare: copy straight out
        if time_idx in save_time_set:
            E_host, H_host, M_host = E_gpu.get(), H_gpu[:, :2 * grid_size].get(), M_gpu.get()
            for c, out in enumerate(outs):
                out["E times picture"][save_index_counter] = E_host[c].reshape(2, grid_size).T
                out["H times picture"][save_index_counter] = H_host[c].reshape(2, grid_size).T
                out["ms time picture"][save_index_counter] = M_host[c].reshape(3, material_size).T
            save_index_counter += 1

    flush()

    results = []
    for case, out, writer in zip(cases, outs, writers):
        data = dict(out)
        data.update({
            "time picture": save_time_indices,
            "locations": save_locations,
            "time intervals": time_interval,
            "matitral location": z_indices, "e_r": case["eps_r"], "conductivity": case["conductivity"],
            "gilbert damping factor": case["damping"], "initial magnetization": M0,
            "max magnetic field": max_E, "polarization phase": phase_diff, "pulse width (FWHM)": pulse_width, "mean wave lenght": wavelength,
            "systen status": case["closed_system"], "dt": dt, "safety_start": safety_factor
        })
        results.append(load_simulation_data(writer.close(data)) if writer else data)
    return results


def simulation(z_indices: list[int], eps_r: list[float], conductivity: list[float], 
                   damping: float, M0: np.array,
                   max_E: float, phase_diff: float, pulse_width: float,
                   save_locations, save_time_indices, closed_system: bool,
                   dt: float, total_time_steps: int, wavelength=800e-9, polarization_angle=0, safety_factor=4, time_interval=1,
                   precision="fp64", output_path=None):
    """One Maxwell-LLG run (a batch of one, see simulation_cases)."""
    case = {"eps_r": eps_r, "conductivity": conductivity, "damping": damping, "closed_system": closed_system}
    return simulation_cases(
        [case], z_indices, M0, max_E, phase_diff, pulse_width, save_locations, save_time_indices,
        dt, total_time_steps, wavelength, polarization_angle, safety_factor, time_interval,
        precision, None if output_path is None else [output_path],
    )[0]


# =============================================================================
# BUNDLE I/O - <name>.simdata/ directory: manifest.json + one .npy per array
# (same format as the partner repo's load_simulation_data)
# =============================================================================
def _to_jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(item) for item in value]
    return value


def bundle_complete(bundle_path) -> bool:
    return (Path(bundle_path).with_suffix(".simdata") / "manifest.json").exists()


class BundleWriter:
    """Arrays are created as .npy memmaps inside <name>.simdata/. manifest.json is
    written last, so a crashed run never looks finished (see bundle_complete)."""

    def __init__(self, bundle_path):
        self.bundle_path = Path(bundle_path).with_suffix(".simdata")
        if self.bundle_path.exists():
            shutil.rmtree(self.bundle_path)
        self.bundle_path.mkdir(parents=True)
        self.manifest = {"arrays": {}, "metadata": {}}
        self.memmaps = []

    def _next_file(self, key):
        file_name = f"array_{len(self.manifest['arrays']):02d}.npy"
        self.manifest["arrays"][key] = file_name
        return self.bundle_path / file_name

    def array(self, key, shape, dtype):
        if int(np.prod(shape)) == 0:  # np.memmap can't map an empty file
            empty = np.zeros(shape, dtype)
            np.save(self._next_file(key), empty)
            return empty
        mm = np.lib.format.open_memmap(self._next_file(key), mode="w+", dtype=dtype, shape=shape)
        self.memmaps.append(mm)
        return mm

    def close(self, data: dict) -> Path:
        """Store the remaining entries of data (small arrays, metadata) and finalize."""
        for mm in self.memmaps:
            mm.flush()
        self.memmaps.clear()
        for key, value in data.items():
            if key in self.manifest["arrays"]:
                continue
            if isinstance(value, np.ndarray):
                np.save(self._next_file(key), value, allow_pickle=False)
            else:
                self.manifest["metadata"][key] = _to_jsonable(value)

        with (self.bundle_path / "manifest.json").open("w") as fh:
            json.dump(self.manifest, fh)
        return self.bundle_path


def save_simulation_bundle(bundle_path, data: dict) -> Path:
    return BundleWriter(bundle_path).close(data)


def load_simulation_data(bundle_path, mmap: bool = True) -> dict:
    """Arrays are memory-mapped by default, so large bundles don't fill RAM."""
    bundle_path = Path(bundle_path)
    with (bundle_path / "manifest.json").open() as fh:
        manifest = json.load(fh)
    data = dict(manifest["metadata"])
    for key, file_name in manifest["arrays"].items():
        data[key] = np.load(bundle_path / file_name, mmap_mode="r" if mmap else None)
    return data
