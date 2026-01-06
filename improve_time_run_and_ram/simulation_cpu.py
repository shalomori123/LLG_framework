import numpy as np
from tqdm import tqdm

# Physical Constants
EPSILON_0 = 8.854e-12  # [F/m]
MU_0 = 4 * np.pi * 1e-7  # [H/m]
SPEED_OF_LIGHT = 1 / (EPSILON_0 * MU_0) ** 0.5
IMPEDANCE = (MU_0 / EPSILON_0) ** 0.5

ELEMENTARY_CHARGE = 1.60217646e-19  # Elementary charge [Coulombs]
LANDAU_FACTOR = 2  # Landau factor
ELECTRON_MASS = 9.1093821545e-31  # Electron mass [kg]
GAMMA_FACTOR = 1
GAMMA = GAMMA_FACTOR * LANDAU_FACTOR * ELEMENTARY_CHARGE / (2 * ELECTRON_MASS)  # [rad/sec*T]


def pulse_gen(max_E, phase_diff, pulse_width, dt, wavelength=800e-9, polarization_angle=0, safety_factor=4):
    """
    input-
    :param max_E: max value of the fields [A/m]
    :param phase_diff: phase difference between E_x and E_y
    :param pulse_width: time from half power to half of optical autorotation [sec].
    :param dt: time step size.
    :param wavelength: wave len in void to the AC amplified
    :param polarization_angle: angle of the central axis of the polarization [rad]
    
    :return: function out_pulse(n) 
    """
    # move to angular frequency
    angular_freq = 2 * np.pi * SPEED_OF_LIGHT / wavelength
    polarization_angle -= np.pi / 4
    rotation_matrix = np.array([np.cos(polarization_angle), -np.sin(polarization_angle), 
                                np.sin(polarization_angle), np.cos(polarization_angle)])
    rotation_matrix = rotation_matrix.reshape((2, 2))
    
    pulse_vector = max_E * rotation_matrix @ np.array([1, np.exp(-phase_diff * 1j)]).reshape(2)
    
    sigma = (0.5 / (np.log(2) ** 0.5)) * pulse_width
    peak_location = safety_factor * sigma

    def out_pulse(time_index: int) -> np.array:
        """
        input-
        :param time_index: index of step in time
        :return: 2D complex np.array of the inject pulse (E_x E_y) in this time step 
        """
        current_time = time_index * dt
        result = pulse_vector * np.exp(-(current_time - peak_location) ** 2 / (2 * sigma ** 2)) * np.exp(1j * angular_freq * current_time)
        return result
    return out_pulse


def LLG_step(M: np.array, H: np.array, dt: float, alpha: float) -> np.array:
    """
    Calculates the next magnetization vector M using the LLG equation in Runge-Kutta 4th Order (RK4 scheme).
    
    :param M: Current magnetization (N x 3).
    :param H: Total effective field (N x 3).
    :param dt: Time step size [s].
    :param alpha: Gilbert damping parameter.
    
    :return: Updated M vector.
    """
    # Magnitude of M (used in the denominator of the damping term)
    # Calculates |M| for all N points (N x 1) and keeps dimensions for broadcasting.
    M_mag = np.linalg.norm(M, axis=1, keepdims=True) # GPU potential

    # Effective LLG terms
    gamma_ll = GAMMA / (1 + alpha**2)
    lambda_ll = gamma_ll * alpha

    def M_derivative(M_guess: np.array) -> np.array: # GPU potential
        """Calculates the instantaneous slope dM/dt for the LLG equation."""
        M_cross_H = np.cross(M_guess, H, axis=1)

        # Rotation term: -gamma_ll * mu * (M x H)
        rotate_term = -gamma_ll * MU_0 * M_cross_H
        # Damping term: -(lambda_ll * mu / M0) * (M x (M x H))
        damping_term = -(lambda_ll * MU_0 / M_mag) * np.cross(M_guess, M_cross_H, axis=1)
        
        return rotate_term + damping_term
    
    # --- RK4 Method Implementation ---
    
    # k1 (a_n): Slope at the beginning (M)
    a_n = M_derivative(M)
    
    # k2 (b_n): Slope at the midpoint (M + dt/2 * k1)
    M_k2 = M + (dt / 2) * a_n
    b_n = M_derivative(M_k2)
    
    # k3 (c_n): Slope at the midpoint (M + dt/2 * k2)
    M_k3 = M + (dt / 2) * b_n
    c_n = M_derivative(M_k3)
    
    # k4 (d_n): Slope at the end (M + dt * k3)
    M_k4 = M + dt * c_n
    d_n = M_derivative(M_k4)
    
    # Final RK4 update: M_next = M + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    M_next = M + (dt / 6) * (a_n + 2*b_n + 2*c_n + d_n) # GPU potential
    return M_next


def simulation(z_indices: list[int], eps_r: list[float], conductivity: list[float], 
               damping: float, M0: np.array,
               max_E: float, phase_diff: float, pulse_width: float,
               save_locations, save_time_indices, closed_system: bool,
               dt: float, total_time_steps: int, wavelength=800e-9, polarization_angle=0, safety_factor=4, time_interval=1):
    """
    inputs-
    :param z_indices: start index of material
    :param eps_r: relative permittivity
    :param conductivity: conductivity of material [Oham/m]
    :param damping: damping param of moments
    :param M0: initial BIG MAGNETIZATION (sum of moments) configuration in [x,y,z] -> vector size 3
    :param max_E: max value of the fields [A/m]
    :param phase_diff: phase difference between E_x and E_y
    :param pulse_width: the width of the gaussian, time when the pulse is greater than half max power [sec].
    :param save_locations: interval of index in space to save all value in the for all time.
    :param save_time_indices: interval of index in time to save all value at the entire space.
    :param closed_system: whether will use the reflection of M on H on E or not.
    :param dt: time step size [s]
    :param total_time_steps: number for time steps to do [int]
    :param wavelength: amplified wave length [m]
    :param polarization_angle: angle of the central axis of the polarization [rad]
    
    :return: two list in len 3
        first list - save parmaters in all space in the save_time stepes
        secnd list - saved paarmates in all time in the seve_location ind on space
        parm are: [H, E, msi ] ( msi given only in the relvent space if axis) 
                    all of them in the shep of [time_step]* [space_step]* [3 demntion]
    """

    # system definition 
    spatial_step = 2 * SPEED_OF_LIGHT * dt
    grid_size = z_indices[-1] + min(int(wavelength / spatial_step), int(0.5 * (z_indices[-1] - z_indices[0])))

    # material property
    coeff_E = np.ones(grid_size - 1)
    coeff_H = 0.5 * np.ones(grid_size - 1)
    for i in range(len(eps_r)):
        current_eps_r = np.ones(z_indices[i+1] - z_indices[i]) 
        current_eps_r += (eps_r[i] - 1) 
        
        # eaf = effective attenuation factor (likely)
        effective_factor = dt / 2 / EPSILON_0 * (conductivity[i]) / current_eps_r
        
        coeff_E[z_indices[i]:z_indices[i+1]] = ((1 - effective_factor) / (1 + effective_factor))
        coeff_H[z_indices[i]:z_indices[i+1]] = (0.5 / current_eps_r / (1 + effective_factor))
    
    # place for run the simulation
    pulse_func = pulse_gen(max_E, phase_diff, pulse_width, dt, wavelength, polarization_angle, safety_factor)

    # defining the work space for the simulation
    E_next = np.zeros((grid_size, 3)).astype(complex)
    E_current = np.zeros((grid_size, 3)).astype(complex)
    
    boundary_buffer_right = [np.zeros((3)), np.zeros((3))]
    boundary_buffer_left = [np.zeros((3)), np.zeros((3))]
    
    H_next = np.zeros((grid_size, 3)).astype(complex)
    H_current = np.zeros((grid_size, 3)).astype(complex)
    
    M_current = np.zeros((z_indices[-1] - z_indices[0], 3))
    
    # initial the beginning magnetization
    for i in range(len(eps_r)):
        mag_segment = M0[i, :]
        start_idx = z_indices[i] - z_indices[0]
        end_idx = z_indices[i+1] - z_indices[0]
        for k in range(3):
            M_current[start_idx:end_idx, k] = mag_segment[k] / (end_idx - start_idx)

    # defining returning array:
    # for space indexes
    total_saved_frames = int(total_time_steps / time_interval + 1)
    E_space_return = np.zeros((total_saved_frames, len(save_locations), 2)).astype(complex)
    H_space_return = np.zeros((total_saved_frames, len(save_locations), 2)).astype(complex)
    M_space_return = np.zeros((total_saved_frames, z_indices[-1] - z_indices[0], 3)).astype(complex)

    # for time indexes
    E_time_return = np.zeros((len(save_time_indices), grid_size, 2)).astype(complex)
    H_time_return = np.zeros((len(save_time_indices), grid_size, 2)).astype(complex)
    M_time_return = np.zeros((len(save_time_indices), z_indices[-1] - z_indices[0], 3)).astype(complex)
    
    save_index_counter = 0
    
    # propagate
    for time_idx in tqdm(range(0, total_time_steps), desc="simulation process: "):
        E_current[3, 0:2] += pulse_func(time_idx) # inject the input signal

        # set boundary condition for not reflection
        r_val = boundary_buffer_right.pop(0)
        l_val = boundary_buffer_left.pop(0)

        E_current[-1, :] = r_val
        E_current[0, :] = l_val
    
        # H equation
        H_next[1:, 0] = H_current[1:, 0] + 0.5 * (E_current[1:, 1] - E_current[:-1, 1])        
        H_next[1:, 1] = H_current[1:, 1] - 0.5 * (E_current[1:, 0] - E_current[:-1, 0])
        
        # LLG and change in the magnetization        
        M_next = LLG_step(M_current, np.real(H_current[z_indices[0]:z_indices[-1], :]), dt, damping)
        
        if closed_system: # integration line for close system B = H + M
            H_next[z_indices[0]:z_indices[-1], :] += M_current - M_next
            
        # E equation
        E_next[:-1, 1] = coeff_E * E_current[:-1, 1] + coeff_H * (H_next[1:, 0] - H_next[:-1, 0])
        E_next[:-1, 0] = coeff_E * E_current[:-1, 0] - coeff_H * (H_next[1:, 1] - H_next[:-1, 1])

        # get ready to next time step
        boundary_buffer_right.append(np.copy(E_current[-2, :]))
        boundary_buffer_left.append(np.copy(E_current[1, :]))
        
        E_current = np.copy(E_next)
        H_current = np.copy(H_next)
        M_current = np.copy(M_next)

        # save to return
        # save picture of parameters in time
        if time_idx in save_time_indices:
            H_time_return[save_index_counter, :, :] = H_current[:, :2]
            E_time_return[save_index_counter, :, :] = E_current[:, :2]
            M_time_return[save_index_counter, :, :] = M_current

            save_index_counter += 1

        # save the parameters in specific location in all times
        if not time_idx % time_interval:
            frame_idx = int(time_idx / time_interval)
            H_space_return[frame_idx, :, :] = H_current[save_locations, :2]
            E_space_return[frame_idx, :, :] = E_current[save_locations, :2]
            M_space_return[frame_idx, :, :] = M_current
            
    all_data = {
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
    return all_data

 

if __name__ == '__main__':
    print("i do this")
    wavelength = 8e-7 
    time_cycle = wavelength / SPEED_OF_LIGHT
    bwhp = 4 * time_cycle 
    max_E = 3e5
    phase_diff = np.pi / 2