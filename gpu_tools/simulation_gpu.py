import json
import os
import shutil
from functools import lru_cache
from pathlib import Path

import numpy as np
import cupy as cp
from tqdm import tqdm

from sim_config import EPSILON_0, MU_0, SPEED_OF_LIGHT, GAMMA


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

    def out_pulse(n: int) -> np.array:
        t = n * dt
        return pulse_amp * np.exp(- (t - peak_loc)**2 / (2 * sigma**2)) * np.exp(1j * omega * t)
        
    return out_pulse


@lru_cache(maxsize=None)
def load_kernels():
    """Compile maxwell_kernel.cu once per process."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, "maxwell_kernel.cu"), "r") as f:
        cuda_code = f.read()
    module = cp.RawModule(code=cuda_code, options=('-use_fast_math', '-arch=sm_89'))
    return module.get_function("magnetic_kernel"), module.get_function("electric_kernel")


def simulation(z_indices: list[int], eps_r: list[float], conductivity: list[float], 
                   damping: float, M0: np.array,
                   max_E: float, phase_diff: float, pulse_width: float,
                   save_locations, save_time_indices, closed_system: bool,
                   dt: float, total_time_steps: int, wavelength=800e-9, polarization_angle=0, safety_factor=4, time_interval=1):
    
    # --- 1. Grid Definition ---
    spatial_step = 2 * SPEED_OF_LIGHT * dt
    material_begin = z_indices[0]
    material_end   = z_indices[-1]
    material_size  = material_end - material_begin
    grid_size = material_end + min(int(wavelength / spatial_step), int(0.5 * material_size))

    # --- 2. LLG Physics Parameters ---
    gamma_ll = GAMMA / (1 + damping**2)
    lambda_ll = gamma_ll * damping

    # The damping term is divided by the per-cell |M|, computed in the kernel each step
    neg_gamma_LL_kernel = -gamma_ll * MU_0
    neg_lambda_mu = -lambda_ll * MU_0

    # --- 3. GPU Memory Allocation (Monolithic SoA) ---
    # Complex arrays for E and H (float2 in CUDA)
    E_gpu = cp.zeros(2 * grid_size, dtype=cp.complex64)
    H_gpu = cp.zeros(3 * grid_size, dtype=cp.complex64)

    # ABC Buffers - 2-step delay at Courant 0.5: [2 slots][Ex, Ey]
    abc_left_gpu = cp.zeros(4, dtype=cp.complex64)
    abc_right_gpu = cp.zeros(4, dtype=cp.complex64)

    # --- 4. Initialize M on GPU ---
    M_cpu_init = np.zeros((material_size, 3), dtype=np.float32)
    for i in range(len(eps_r)):
        start_idx = z_indices[i] - material_begin
        end_idx = z_indices[i+1] - material_begin
        for k in range(3):
            M_cpu_init[start_idx:end_idx, k] = M0[i, k] / (end_idx - start_idx)

    # Convert AoS to Monolithic SoA and upload to GPU
    M_soa_cpu = np.concatenate([M_cpu_init[:, 0], M_cpu_init[:, 1], M_cpu_init[:, 2]])
    M_gpu = cp.asarray(M_soa_cpu)

    # --- 5. Pre-calculate the Injection Pulse ---
    # Doing this on the CPU once prevents a massive GPU/CPU sync bottleneck during the loop
    pulse_func = pulse_gen(max_E, phase_diff, pulse_width, dt, wavelength, polarization_angle, safety_factor)
    pulse_array = np.array([pulse_func(t) for t in range(total_time_steps)], dtype=np.complex64)
    
    pulse_x_gpu = cp.asarray(pulse_array[:, 0])
    pulse_y_gpu = cp.asarray(pulse_array[:, 1])

    # --- 6. Load CUDA Kernels ---
    magnetic_kernel, electric_kernel = load_kernels()

    # Block and Grid dimensions
    threads_per_block = 256
    blocks_per_grid = (grid_size + threads_per_block - 1) // threads_per_block

    # --- 7. Output Buffers (CPU) ---
    total_saved_frames = int(total_time_steps / time_interval + 1)
    E_space_return = np.zeros((total_saved_frames, len(save_locations), 2), dtype=np.complex128)
    H_space_return = np.zeros((total_saved_frames, len(save_locations), 2), dtype=np.complex128)
    M_space_return = np.zeros((total_saved_frames, material_size, 3), dtype=np.float64)

    E_time_return = np.zeros((len(save_time_indices), grid_size, 2), dtype=np.complex128)
    H_time_return = np.zeros((len(save_time_indices), grid_size, 2), dtype=np.complex128)
    M_time_return = np.zeros((len(save_time_indices), material_size, 3), dtype=np.float64)
    
    save_index_counter = 0

    # Faraday coefficient (vacuum, Courant 0.5)
    coeff = cp.float32(0.5)

    # Per-cell Ampere coefficients for lossy dielectric [Sullivan eq. 3.8]
    coeff_E_cpu = np.ones(grid_size, dtype=np.float32)
    coeff_H_cpu = np.full(grid_size, 0.5, dtype=np.float32)
    for i in range(len(eps_r)):
        eaf = dt * conductivity[i] / (2 * EPSILON_0 * eps_r[i])
        coeff_E_cpu[z_indices[i]:z_indices[i+1]] = (1 - eaf) / (1 + eaf)
        coeff_H_cpu[z_indices[i]:z_indices[i+1]] = 0.5 / (eps_r[i] * (1 + eaf))
    coeff_E_gpu = cp.asarray(coeff_E_cpu)
    coeff_H_gpu = cp.asarray(coeff_H_cpu)

    # --- 8. THE MAIN FDTD LOOP ---
    for time_idx in tqdm(range(0, total_time_steps), desc="GPU Simulation Process: "):
        
        # 1. Inject the input signal directly on the GPU
        E_gpu[3] += pulse_x_gpu[time_idx]
        E_gpu[grid_size + 3] += pulse_y_gpu[time_idx]

        # 2. Magnetic Kernel (Applies ABC, updates H, steps LLG, applies Coupling)
        magnetic_kernel(
            (blocks_per_grid,), (threads_per_block,),
            (E_gpu, H_gpu, M_gpu, 
             np.int32(grid_size), coeff, 
             np.int32(material_begin), np.int32(material_end), 
             cp.float32(dt), cp.float32(neg_gamma_LL_kernel), cp.float32(neg_lambda_mu),
             np.int32(closed_system), abc_left_gpu, abc_right_gpu)
        )

        # 3. Electric Kernel (Updates E, fills ABC buffers)
        electric_kernel(
            (blocks_per_grid,), (threads_per_block,),
            (E_gpu, H_gpu, np.int32(grid_size), coeff_E_gpu, coeff_H_gpu, abc_left_gpu, abc_right_gpu)
        )

        # --- 9. Data Extraction ---
        # (We pull data back to the CPU ONLY when needed to avoid bottlenecking the GPU)
        if time_idx in save_time_indices or (time_idx % time_interval == 0):
            # Sync GPU before reading
            cp.cuda.Stream.null.synchronize()
            
            # Reconstruct AoS from SoA
            E_cpu = np.stack([E_gpu[:grid_size].get(), E_gpu[grid_size:2*grid_size].get()], axis=1)
            H_cpu = np.stack([H_gpu[:grid_size].get(), H_gpu[grid_size:2*grid_size].get()], axis=1)
            M_cpu = np.stack([M_gpu[:material_size].get(), 
                              M_gpu[material_size:2*material_size].get(), 
                              M_gpu[2*material_size:].get()], axis=1)

            if time_idx in save_time_indices:
                H_time_return[save_index_counter, :, :] = H_cpu
                E_time_return[save_index_counter, :, :] = E_cpu
                M_time_return[save_index_counter, :, :] = M_cpu
                save_index_counter += 1

            if not time_idx % time_interval:
                frame_idx = int(time_idx / time_interval)
                H_space_return[frame_idx, :, :] = H_cpu[save_locations, :]
                E_space_return[frame_idx, :, :] = E_cpu[save_locations, :]
                M_space_return[frame_idx, :, :] = M_cpu
            
    return {
        "H times picture": H_time_return,
        "E times picture": E_time_return,
        "ms time picture": M_time_return,
        "time picture": save_time_indices,
        "H in locations": H_space_return,
        "E in locations": E_space_return,
        "ms in locations": M_space_return,
        "locations": save_locations,
        "time intervals": time_interval,
        "matitral location": z_indices, "e_r": eps_r, "conductivity": conductivity, "gilbert damping factor": damping, "initial magnetization": M0,
        "max magnetic field": max_E, "polarization phase": phase_diff, "pulse width (FWHM)": pulse_width, "mean wave lenght": wavelength,
        "systen status": closed_system, "dt": dt, "safety_start": safety_factor
    }


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


def save_simulation_bundle(bundle_path, data: dict) -> Path:
    bundle_path = Path(bundle_path).with_suffix(".simdata")
    temp_path = bundle_path.with_name(bundle_path.name + ".tmp")
    if temp_path.exists():
        shutil.rmtree(temp_path)
    temp_path.mkdir(parents=True)

    manifest = {"arrays": {}, "metadata": {}}
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            file_name = f"array_{len(manifest['arrays']):02d}.npy"
            np.save(temp_path / file_name, value, allow_pickle=False)
            manifest["arrays"][key] = file_name
        else:
            manifest["metadata"][key] = _to_jsonable(value)

    with (temp_path / "manifest.json").open("w") as fh:
        json.dump(manifest, fh)
    if bundle_path.exists():
        shutil.rmtree(bundle_path)
    temp_path.rename(bundle_path)
    return bundle_path


def load_simulation_data(bundle_path, mmap: bool = True) -> dict:
    """Arrays are memory-mapped by default, so large bundles don't fill RAM."""
    bundle_path = Path(bundle_path)
    with (bundle_path / "manifest.json").open() as fh:
        manifest = json.load(fh)
    data = dict(manifest["metadata"])
    for key, file_name in manifest["arrays"].items():
        data[key] = np.load(bundle_path / file_name, mmap_mode="r" if mmap else None)
    return data
