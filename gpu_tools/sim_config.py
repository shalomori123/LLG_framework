"""
All simulation configuration: physical constants, materials and presets.

Presets hold the primitive inputs; build_preset() derives dt, grid indexes,
save locations/times and the simulation length (same formulas as the original
run_material.py). Edit a value here and rerun sim_run.py.
"""
import numpy as np

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================
EPSILON_0 = 8.854e-12                               # [F/m]
MU_0 = 4 * np.pi * 1e-7                             # [H/m]
SPEED_OF_LIGHT = 1 / (EPSILON_0 * MU_0) ** 0.5      # [m/s]
IMPEDANCE = (MU_0 / EPSILON_0) ** 0.5               # [Ohm]

ELEMENTARY_CHARGE = 1.60217646e-19                  # [C]
LANDAU_FACTOR = 2
ELECTRON_MASS = 9.1093821545e-31                    # [kg]
GAMMA_FACTOR = 1
GAMMA = GAMMA_FACTOR * LANDAU_FACTOR * ELEMENTARY_CHARGE / (2 * ELECTRON_MASS)  # [rad/sec*T]

# =============================================================================
# MATERIALS - magnetic layer + conductive cap
# =============================================================================
MAGNETIC = {                # eps_r, sigma [S/m]
    "Co": (1, 1.7e7),
    "Fe": (1, 1e7),
    "Ni": (1, 1.4e7),
}
CONDUCTIVE = {
    "Au": (1, 1e6),
    "Pt": (1, 7e6),
}
DAMPING = {                 # Gilbert damping per pair
    "Co_Au": 0.02,  "Co_Pt": 0.025,
    "Fe_Au": 0.02,  "Fe_Pt": 0.025,
    "Ni_Au": 0.05,  "Ni_Pt": 0.06,
}
ALL_MATERIALS = list(DAMPING)
STATES = {"open": False, "close": True}


def material_case(name: str) -> dict:
    """'Co_Pt' -> per-layer eps_r and conductivity lists, plus damping."""
    magnetic, conductive = name.split("_")
    (eps_m, sigma_m), (eps_c, sigma_c) = MAGNETIC[magnetic], CONDUCTIVE[conductive]
    return {"eps_r": [eps_m, eps_c], "conductivity": [sigma_m, sigma_c], "damping": DAMPING[name]}


# =============================================================================
# PRESETS
# =============================================================================
COMMON = {
    "mean_wavelength": 784e-9,          # [m]
    "b_opt": 2.4e5,                     # [A/m] peak optical field
    "phase_diff": np.pi / 2,            # circular polarization
    "safety_start": 3,                  # pulse peak enters after this many sigma
    "magnetic_thickness": 10e-9,        # [m]
    "conductive_thickness": 2e-9,       # [m]
    "ms0_amplitude": 3e5,               # [A/m]
}

PRESETS = {
    # Subwavelength grid, resolves the in-material magnetization dynamics
    "full_width": dict(COMMON,
        fwhm=1.1e-12, dz=2e-9 / 8, z1=8,
        time_interval=7, sim_time_factor=3,
        ms0_vectors=[[1, 0, 0], [0, 0, -0.0001]],
        save_time_fracs=(0.5, 0.8, 20),
    ),
    # Return & transmission included: long vacuum gap around the material
    "rt": dict(COMMON,
        fwhm=0.22e-12, dz=2e-9 / 2, z1_sigma_factor=2.15,
        time_interval=6, sim_time_factor=2.5,
        ms0_vectors=[[1, 0, 0], [0, 0, -0.00001]],
        save_time_fracs=(0.0, 1.0, 300),
    ),
}


def grid_size(z_indexes, dt: float, wavelength: float) -> int:
    """Number of cells: material end + a vacuum tail after it."""
    spatial_step = 2 * SPEED_OF_LIGHT * dt
    material_size = z_indexes[-1] - z_indexes[0]
    return z_indexes[-1] + min(int(wavelength / spatial_step), int(0.5 * material_size))


def build_preset(name: str, total_time_steps: int = None) -> dict:
    """Derive all run parameters of a preset. total_time_steps truncates the run
    (save times are spread over the truncated length)."""
    p = PRESETS[name]
    dz = p["dz"]
    dt = dz / (2 * SPEED_OF_LIGHT)
    time_cycle = p["mean_wavelength"] / SPEED_OF_LIGHT

    magnetic_steps = p["magnetic_thickness"] // dz
    conductive_steps = p["conductive_thickness"] // dz

    sigma_pulse = ((0.5 / np.log(2)) ** 0.5) * p["fwhm"]
    peak_enter_time = p["safety_start"] * sigma_pulse

    z1 = p["z1"] if "z1" in p else int(p["z1_sigma_factor"] * sigma_pulse / dt)
    z_end = int(z1 + magnetic_steps + conductive_steps)
    z_indexes = [z1, int(z1 + magnetic_steps), z_end]

    if name == "rt":
        save_locations = np.array([4] + z_indexes + [0, z_end + 1])
        save_locations[3:5] = save_locations[2:4]
        save_locations[2] = int(z1 + magnetic_steps - 1)
    else:
        save_locations = np.array([4] + z_indexes)

    void_steps = 2 * z1
    exit_steps = int(time_cycle / dt)
    simulation_time = int(p["sim_time_factor"] * (peak_enter_time // dt + void_steps
                          + 4 * (magnetic_steps + conductive_steps)) + exit_steps)
    if total_time_steps is not None:
        simulation_time = int(total_time_steps)

    start, end, count = p["save_time_fracs"]
    save_time = np.linspace(simulation_time * start, simulation_time * end, count).astype(int)
    save_time = save_time[save_time < simulation_time]

    return {
        "preset": name,
        "dt": dt,
        "z_indexes": z_indexes,
        "grid_size": grid_size(z_indexes, dt, p["mean_wavelength"]),
        "save_locations": save_locations,
        "save_time": save_time,
        "simulation_time": simulation_time,
        "time_interval": p["time_interval"],
        "ms0": p["ms0_amplitude"] * np.array(p["ms0_vectors"], dtype=np.float64),
        "max_E": p["b_opt"],
        "phase_diff": p["phase_diff"],
        "fwhm": p["fwhm"],
        "wavelength": p["mean_wavelength"],
        "safety_start": p["safety_start"],
    }
