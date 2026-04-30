import numpy as np
import cupy as cp
from tqdm import tqdm
import os

# Physical Constants 
EPSILON_0 = 8.854e-12
MU_0 = 4 * np.pi * 1e-7
SPEED_OF_LIGHT = 1 / (EPSILON_0 * MU_0) ** 0.5
IMPEDANCE = (MU_0 / EPSILON_0) ** 0.5

ELEMENTARY_CHARGE = 1.60217646e-19
LANDAU_FACTOR = 2
ELECTRON_MASS = 9.1093821545e-31
GAMMA_FACTOR = 1
GAMMA = GAMMA_FACTOR * LANDAU_FACTOR * ELEMENTARY_CHARGE / (2 * ELECTRON_MASS)


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
    
    # Calculate average magnitude per cell for the damping scalar
    # (Matches your logic of dividing M0 by the segment size)
    cell_M_mag = np.linalg.norm(M0[0, :]) / material_size 
    
    neg_gamma_LL_kernel = -gamma_ll * MU_0
    neg_coeff_damp = -(lambda_ll * MU_0) / cell_M_mag

    # --- 3. GPU Memory Allocation (Monolithic SoA) ---
    # Complex arrays for E and H (float2 in CUDA)
    E_gpu = cp.zeros(2 * grid_size, dtype=cp.complex64)
    H_gpu = cp.zeros(3 * grid_size, dtype=cp.complex64)
    
    # Real array for M (float in CUDA)
    M_gpu = cp.zeros(3 * material_size, dtype=cp.float32)

    # ABC Buffers
    abc_left_gpu = cp.zeros(2, dtype=cp.complex64)
    abc_right_gpu = cp.zeros(2, dtype=cp.complex64)

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
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, "maxwell_kernel.cu"), "r") as f:
        cuda_code = f.read()

    # We pass '-I.' so nvcc finds your 'llg_device.cuh' file in the same folder
    module = cp.RawModule(code=cuda_code, options=('-I.', '-use_fast_math', '-arch=sm_89'))
    magnetic_kernel = module.get_function("magnetic_kernel")
    electric_kernel = module.get_function("electric_kernel")

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

    # Constant Coefficient (Note: To use array coefficients, update the kernel signature)
    coeff = cp.float32(0.5)

    # --- 8. THE MAIN FDTD LOOP ---
    for time_idx in tqdm(range(0, total_time_steps), desc="GPU Simulation Process: "):
        
        # 1. Inject the input signal directly on the GPU
        E_gpu[3] += pulse_x_gpu[time_idx]
        E_gpu[grid_size + 3] += pulse_y_gpu[time_idx]

        # 2. Magnetic Kernel (Updates H, steps LLG, applies Coupling)
        magnetic_kernel(
            (blocks_per_grid,), (threads_per_block,),
            (E_gpu, H_gpu, M_gpu, 
             grid_size, coeff, 
             material_begin, material_end, 
             cp.float32(dt), cp.float32(neg_gamma_LL_kernel), cp.float32(neg_coeff_damp))
        )

        # 3. Electric Kernel (Updates E, manages ABC boundaries)
        electric_kernel(
            (blocks_per_grid,), (threads_per_block,),
            (E_gpu, H_gpu, grid_size, coeff, abc_left_gpu, abc_right_gpu)
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

if __name__ == '__main__':
    print("Initializing GPU Simulation...")
    # Add your test variables here and run simulation_gpu(...)