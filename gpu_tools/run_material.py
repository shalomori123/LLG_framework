import numpy as np
import os
import simulation_gpu as sim

# =============================================================================
# CO-PT MATERIAL PROPERTIES & CONSTANTS
# =============================================================================
epsr_list = [1, 1]              # Co_er, Pt_er
sigma_list = [1.7e7, 7e6]       # Co_sigma, Pt_sigma
alpha_val = 0.025               # Gilbert damping for Co-Pt

mean_wave_lenght = 784e-9       # [m]
B_opt = 2.4e5                   # [A/m]
fie = np.pi / 2
sefe_start = 3

print("=== Starting Co-Pt Focused Simulations ===")

# =============================================================================
# SIMULATION 1: FULL WIDTH PULSE (Co-Pt Only)
# =============================================================================
print("\n--- Running Simulation 1: Full Width Pulse ---")

FWHM_1 = 1.1e-12 
path_full_width = "simulation_result/Full_Width_Pulse/"
os.makedirs(path_full_width, exist_ok=True)

# Grid and Time Steps
dz_1 = 2e-9 / 8
dt_1 = dz_1 / (2 * sim.SPEED_OF_LIGHT)
time_cycal_1 = mean_wave_lenght / sim.SPEED_OF_LIGHT

magntic_steps_1 = 10e-9 // dz_1
conductive_steps_1 = 2e-9 // dz_1

# Pulse properties
sigma_pulse_1 = ((0.5 / np.log(2)) ** 0.5) * FWHM_1
peak_enter_time_1 = sefe_start * sigma_pulse_1

# Initial Magnetization (Specific to Part 1)
ms0_1 = 3e5 * np.array([1, 0, 0, 0, 0, -0.0001]).reshape((2, 3))

# Spatial Bounds
z1_1 = 8
z_end_1 = int(z1_1 + magntic_steps_1 + conductive_steps_1)
z_indexes_1 = [z1_1, int(z1_1 + magntic_steps_1), z_end_1]
save_loc_1 = np.array([4] + z_indexes_1)

void_steps_1 = 2 * z1_1
exit_steps_1 = int(time_cycal_1 / dt_1)
sim_time_1 = int(3 * (peak_enter_time_1 // dt_1 + void_steps_1 + 4 * (magntic_steps_1 + conductive_steps_1)) + exit_steps_1)
save_time_1 = np.linspace(sim_time_1 * 0.5, sim_time_1 * 0.8, 20).astype("int")

for is_closed in [False, True]:
    state = "close" if is_closed else "open"
    print(f"Running Co_Pt (Sim 1) - State: {state.upper()}")
    
    file_name = f"{path_full_width}Co_Pt_{state}.npy"
    
    data = sim.simulation(
        z_indexes_1, epsr_list, sigma_list, alpha_val, ms0_1,
        B_opt, fie, FWHM_1, save_loc_1, save_time_1, is_closed, dt_1, 
        sim_time_1, wavelength=mean_wave_lenght, safety_factor=sefe_start, time_interval=7
    )
    np.save(file_name, data)


# =============================================================================
# SIMULATION 2: RETURN & TRANSMISSION INCLUDED (Co-Pt Only)
# =============================================================================
print("\n--- Running Simulation 2: Return & Transmission ---")

FWHM_2 = 0.22e-12 
path_return_trans = "simulation_result/Return_Tranmision_Incloud/"
os.makedirs(path_return_trans, exist_ok=True)

# Grid and Time Steps
dz_2 = 2e-9 / 2
dt_2 = dz_2 / (2 * sim.SPEED_OF_LIGHT)
time_cycal_2 = mean_wave_lenght / sim.SPEED_OF_LIGHT

magntic_steps_2 = 10e-9 // dz_2
conductive_steps_2 = 2e-9 // dz_2

# Pulse properties
sigma_pulse_2 = ((0.5 / np.log(2)) ** 0.5) * FWHM_2
peak_enter_time_2 = sefe_start * sigma_pulse_2

# Initial Magnetization (Specific to Part 2)
ms0_2 = 3e5 * np.array([1, 0, 0, 0, 0, -0.00001]).reshape((2, 3))

# Spatial Bounds
z1_2 = int(2.15 * sigma_pulse_2 / dt_2)
z_end_2 = int(z1_2 + magntic_steps_2 + conductive_steps_2)
z_indexes_2 = [z1_2, int(z1_2 + magntic_steps_2), z_end_2]

# Exact save_loc logic from your original notebook
save_loc_2 = np.array([4] + z_indexes_2 + [0, z_end_2 + 1])
save_loc_2[3:5] = save_loc_2[2:4]
save_loc_2[2] = int(z1_2 + magntic_steps_2 - 1)

void_steps_2 = 2 * z1_2
exit_steps_2 = int(time_cycal_2 / dt_2)
sim_time_2 = int(2.5 * (peak_enter_time_2 // dt_2 + void_steps_2 + 4 * (magntic_steps_2 + conductive_steps_2)) + exit_steps_2)

for is_closed in [False, True]:
    state = "close" if is_closed else "open"
    print(f"Running Co_Pt (Sim 2) - State: {state.upper()}")
    
    file_name = f"{path_return_trans}Co_Pt_{state}.npy"

    data = sim.simulation(
        z_indexes_2, epsr_list, sigma_list, alpha_val, ms0_2,
        B_opt, fie, FWHM_2, save_loc_2, [], is_closed, dt_2, 
        sim_time_2, wavelength=mean_wave_lenght, safety_factor=sefe_start, time_interval=6
    )
    np.save(file_name, data)

print("\n=== Co-Pt Simulations Completed Successfully ===")