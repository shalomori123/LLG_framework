import numpy as np
import matplotlib.pyplot as plt
import simulation_gpu as sim
import scipy.special as scp
from ..python_core import analyze_tool as an
import os

# =============================================================================
# MATERIAL PROPERTIES
# =============================================================================
# Cobalt
Co_er = 1
Co_sigma = 1.7e7 
Co_er_complex = -16 + 23.3j 

# Iron
Fe_er = 1
Fe_sigma = 1e7

# Nickel
Ni_er = 1
Ni_sigma = 1.4e7 

alpae_renge = (0.025, 0.0025) # general range
gilbert_damping_time = 600e-12 # [sec] 
precession_frequency = 8.5e9 # [Hz]

# Additional materials property
Pt_er = 1
Pt_sigma = 7e6 # [1/ohm*m]

Au_er = 1
Au_sigma = 1e6

condac = {"Co": Co_sigma, "Fe": Fe_sigma, "Ni": Ni_sigma}
alpha = {"Co_Au": 0.02, "Co_Pt": 0.025, "Fe_Au": 0.02, "Fe_Pt": 0.025, "Ni_Au": 0.05, "Ni_Pt": 0.06}

key_name = [
    "H times picture", "E times picture", "ms time picture", "time picture",
    "H in locations", "E in locations", "ms in locations", "locations",
    "time intervals", "matitral location", "e_r", "conductivity", 
    "gilbert damping factor", "initial magnetization", "max magnetic field",
    "polarization phase", "pulse width (FWHM)", "systen status", "dt", "safety_start"                
]

# =============================================================================
# SIMULATION 1: FULL WIDTH PULSE
# =============================================================================
print("Starting Simulation 1: Full Width Pulse...")

mean_wave_lenght = 784e-9 # [m]
FWHM = 1.1e-12 # [sec]
B_opt = 2.4e5 # [A/m] = 0.3[T]

w = 2 * np.pi / (mean_wave_lenght / sim.SPEED_OF_LIGHT)
path_full_width = "C:\\maxwell-LLG\\ws\\simulation_result\\Full_Width_Pulse\\"
os.makedirs(path_full_width, exist_ok=True)

# Basic stuff
time_cycal = mean_wave_lenght / sim.SPEED_OF_LIGHT
dz = 2e-9 / 8
dt = dz / (2 * sim.SPEED_OF_LIGHT)

# Material length according to article
magntic_steps = 10e-9 // dz
conductive_steps = 2e-9 // dz

# Pulse properties
sefe_start = 3
sigma_pulse = ((0.5 / np.log(2)) ** 0.5) * FWHM
peak_enter_time = sefe_start * sigma_pulse
fie = np.pi / 2

# Material properties arrays
magntic_name = ["Co", "Fe", "Ni"]
conductive_name = ["Au", "Pt"]

magntic_eps = [Co_er, Fe_er, Ni_er]
conductive_eps = [Au_er, Pt_er]

magntic_sigma = [Co_sigma, Fe_sigma, Ni_sigma]
conductive_sigma = [Au_sigma, Pt_sigma]

ms0 = 3e5 * np.array([1, 0, 0, 0, 0, -0.0001]).reshape((2, 3))

# Simulation building blocks
z1 = 8
z_end = int(z1 + magntic_steps + conductive_steps)
z_indexes = [z1, int(z1 + magntic_steps), z_end]

save_loc = np.array([4] + z_indexes)
save_time = []

void_steps = 2 * z1
exit_steps = int(time_cycal / dt)

# Simulation time
simulation_time = int(3 * (peak_enter_time // dt + void_steps + 4 * (magntic_steps + conductive_steps)) + exit_steps)
save_time = np.linspace(simulation_time * 0.5, simulation_time * 0.8, 20).astype("int")

for mag in range(3):
    for cond in range(2):
        for close_sistem in [False, True]:
            file_name = path_full_width + magntic_name[mag] + "_" + conductive_name[cond] 
            file_name += "_close.npy" if close_sistem else "_open.npy"

            epsr_list = [magntic_eps[mag], conductive_eps[cond]]
            sigma_list = [magntic_sigma[mag], conductive_sigma[cond]]

            data = sim.simulation(
                z_indexes, epsr_list, sigma_list, 
                alpha[magntic_name[mag] + "_" + conductive_name[cond]], ms0,
                B_opt, fie, FWHM, save_loc, save_time, close_sistem, dt, 
                simulation_time, wavelength=mean_wave_lenght, safety_factor=sefe_start, time_interval=7
            )
            np.save(file_name, data)


# =============================================================================
# SIMULATION 2: RETURN & TRANSMISSION INCLUDED
# =============================================================================
print("Starting Simulation 2: Return Transmission Included...")

FWHM = 0.22e-12 # [sec]
path_return_trans = "C:\\maxwell-LLG\\ws\\simulation_result\\Return_Tranmision_Incloud\\"
os.makedirs(path_return_trans, exist_ok=True)

# Basic stuff
lamda = 8e-7 
dz = 2e-9 / 2
dt = dz / (2 * sim.SPEED_OF_LIGHT)

magntic_steps = 10e-9 // dz
conductive_steps = 2e-9 // dz

# Pulse properties
sigma_pulse = ((0.5 / np.log(2)) ** 0.5) * FWHM
peak_enter_time = sefe_start * sigma_pulse

# Note: Magnetic order changed here
magntic_name = ["Co", "Ni", "Fe"] 
magntic_eps = [Co_er, Ni_er, Fe_er]
magntic_sigma = [Co_sigma, Ni_sigma, Fe_sigma]

ms0 = 3e5 * np.array([1, 0, 0, 0, 0, -0.00001]).reshape((2, 3))

z1 = int(2.15 * sigma_pulse / dt)
z_end = int(z1 + magntic_steps + conductive_steps)
z_indexes = [z1, int(z1 + magntic_steps), z_end]

save_loc = np.array([4] + z_indexes + [0, z_end + 1])
save_loc[3:5] = save_loc[2:4]
save_loc[2] = int(z1 + magntic_steps - 1)
save_time = []

void_steps = 2 * z1
exit_steps = int(time_cycal / dt)

simulation_time = int(2.5 * (peak_enter_time // dt + void_steps + 4 * (magntic_steps + conductive_steps)) + exit_steps)

for mag in range(3):
    for cond in range(2):
        for close_sistem in [False, True]:
            file_name = path_return_trans + magntic_name[mag] + "_" + conductive_name[cond]
            file_name += "_close.npy" if close_sistem else "_open.npy"
            
            epsr_list = [magntic_eps[mag], conductive_eps[cond]]
            sigma_list = [magntic_sigma[mag], conductive_sigma[cond]]

            data = sim.simulation(
                z_indexes, epsr_list, sigma_list, 
                alpha[magntic_name[mag] + "_" + conductive_name[cond]], ms0,
                B_opt, fie, FWHM, save_loc, save_time, close_sistem, dt, 
                simulation_time, wavelength=mean_wave_lenght, safety_factor=sefe_start, time_interval=6
            )
            np.save(file_name, data)


# =============================================================================
# DATA ANALYSIS & SUMMARIZE
# =============================================================================
print("Analyzing Data...")
ms_to_ms_for_all = np.zeros((12, z_indexes[1] - z_indexes[0], 3))
fie_for_all = np.zeros((12, z_indexes[1] - z_indexes[0]))

i = 0
for mag in range(3):
    for cond in range(2):
        for close_sistem in [False, True]:
            file_name = path_return_trans + magntic_name[mag] + "_" + conductive_name[cond]
            file_name += "_close.npy" if close_sistem else "_open.npy"
 
            data = np.load(file_name, allow_pickle=True).item()
            ms = np.real(data["ms in locations"])
            z2_in_ms = data["matitral location"][1] - data["matitral location"][0]

            ms_to_ms_for_all[i, :, :] = ms[-1, :z2_in_ms, :]
            i += 1

data_summ = {"mz": ms_to_ms_for_all[:, :, 2], "my": ms_to_ms_for_all[:, :, 1], "mx": ms_to_ms_for_all[:, :, 0]}
np.save(path_return_trans + "summrize_true_FWHM.npy", data_summ)

# =============================================================================
# PLOTTING
# =============================================================================
mz_to_ms_for_all = np.zeros((12, z_indexes[1] - z_indexes[0]))
H = data["H in locations"]
ms = data["ms in locations"]

sim_cycle_steps = int(time_cycal / (data["dt"] * data["time intervals"]))
mz_to_m0_in_every_time_step_every_dz = np.real(ms[-2*sim_cycle_steps:, :z2_in_ms, 2]) / np.real(ms[0, :z2_in_ms, 0])
mz_to_m0_in_every_dz = np.mean(mz_to_m0_in_every_time_step_every_dz, 0) 

ms_end_itration = np.mean(ms[-2*sim_cycle_steps:, :z2_in_ms, :], 0)
fie_of_mz = (180 / np.pi) * np.arctan(((ms_end_itration[:, 0]**2 + ms_end_itration[:, 1]**2)**0.5) / ms_end_itration[:, 2])

t = np.arange(np.size(ms, 0)) * (data["dt"] * data["time intervals"])

plt.figure()
plt.plot(t, np.real(ms[:, 0, 2]), label="z")
plt.plot(t, np.real(ms[:, 0, 1]), label="y")
plt.legend()
plt.show()

# # =============================================================================
# # TRANSMISSION & REFLECTION CALCULATION
# # =============================================================================
# print("Calculating Transmission & Reflection...")
# S_transmission = np.zeros((6, 2))
# S_return = np.zeros((6, 2))
# S_pass = np.zeros((6, 2))
# S_after_material = np.zeros((6, 2))
# S_bitween_materials = np.zeros((6, 2))

# i = 0
# for mag in range(3):
#     for cond in range(2):
#         for close_sistem in [False, True]:
#             file_name = path_return_trans + magntic_name[mag] + "_" + conductive_name[cond]
#             if close_sistem:
#                 file_name += "_close.npy"
#                 c = 0
#             else:
#                 file_name += "_open.npy"
#                 c = 1
            
#             data = np.load(file_name, allow_pickle=True).item()
#             data_dt = data["dt"] * data["time intervals"]
#             matiral_l = data["matitral location"]
#             save_indxes = list(data["locations"])

#             z_enter = 0
#             z_matirel = save_indxes.index(matiral_l[0])
#             z_out = save_indxes.index(matiral_l[-1] + 1)
#             z_end_magnet = save_indxes.index(matiral_l[1])
            
#             E_file = data["E in locations"]
#             void_steps = 2 * data["matitral location"][0] / data["time intervals"]
#             t1 = int(peak_enter_time / data_dt + void_steps)
#             t_hm = int((sefe_start * 0.5 * (np.log(2)**-0.5) * FWHM) / data_dt)
 
#             E = np.zeros((np.size(E_file, 0), np.size(E_file, 1), 3))
#             H = np.zeros((np.size(E_file, 0), np.size(E_file, 1), 3))

#             # Note: assuming sim.heta exists in your simulation_cpu based on the notebook
#             if hasattr(sim, 'heta'):
#                 E[:, :, 0:2] = sim.heta * np.real(E_file)
#             else:
#                 # Fallback if it was a typo for IMPEDANCE in the original code
#                 E[:, :, 0:2] = sim.IMPEDANCE * np.real(E_file)
                
#             H[:, :, 0:2] = np.real(data["H in locations"])

#             S = np.cross(E, H)
#             S_in_enter_point = S[:, z_enter, 2]
#             S_in_z1 = S[:, z_matirel, 2]
             
#             time = np.arange(0, np.size(E_file, 0)) * data_dt

#             S_transmission[i, c] = sum(S_in_enter_point[:t1]) * data_dt
#             S_return[i, c] = sum(S_in_enter_point[t1:]) * data_dt
#             S_pass[i, c] = sum(S_in_z1) * data_dt
#             S_after_material[i, c] = sum(S[:, z_out, 2]) * data_dt
#             S_bitween_materials[i, c] = sum(S[:, z_end_magnet, 2]) * data_dt
                        
#         i += 1

# transmition = abs(S_pass / S_transmission)
# Return = abs(S_return / S_transmission)
# chack = (transmition + Return).T

# open_vs_close_T = (transmition[:, 1] - transmition[:, 0])
# open_vs_close_R = (Return[:, 1] - Return[:, 0])

# print("**********Return**********\n", Return.T)
# print("**********Transmition**********\n", transmition.T)
# print("**********Open - Close T**********\n", open_vs_close_T.T)
# print("**********Open - Close R**********\n", open_vs_close_R.T)
# print("**********T+R**********\n", chack)

# d = {"Transmition": open_vs_close_T, "Return": open_vs_close_R, "a": [0.02, 0.025, 0.02, 0.025, 0.05, 0.06]}
# np.save(path_return_trans + "plot_RT_vs_alpha.npy", d)

# plt.figure()
# plt.plot(transmition, open_vs_close_T, label="Transmission Delta")
# plt.plot(Return, open_vs_close_R, label="Return Delta")
# plt.legend()
# plt.show()

# print("Script completed successfully.")