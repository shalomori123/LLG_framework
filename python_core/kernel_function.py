import numpy as np
from numba import cuda

# Physical constants
q = 1.60217646e-19    # Elementary charge [Coulombs]
miu = 4 * np.pi * 1e-7    # Magnetic permeability [H/m]
g = 2    # Landau factor
me = 9.1093821545e-31    # Electron mass [kg]
gma_factor = 1
gma = gma_factor * g * q / (2 * me)    # Gyromagnetic ratio [rad/sec*T]

# Define the GPU kernel for Runge-Kutta LLG algorithm
@cuda.jit
def runge_kutta_LLG_kernel(M, H, dt, alpha, result):
    """
    CUDA kernel to perform one step of the Runge-Kutta LLG algorithm on GPU.
    
    Parameters:
    -----------
    M : numpy.ndarray
        Array of initial magnetization vectors (3D vectors).
    H : numpy.ndarray
        Array of applied magnetic field vectors (3D vectors).
    dt : float
        Time step size.
    alpha : float
        Damping coefficient.
    result : numpy.ndarray
        Array to store the resulting magnetization after the RK step.
    """
    idx = cuda.grid(1)
    
    if idx < M.shape[0]:
        # Extract components of the vectors
        Mx, My, Mz = M[idx]
        Hx, Hy, Hz = H[idx]
        
        M0 = np.sqrt(Mx**2 + My**2 + Mz**2)    # Compute magnitude of M
        
        gma_LL = gma / (1 + alpha**2)
        LL_lamda = gma * alpha / (1 + alpha**2)
        
        # Compute intermediate steps of RK4 method
        a1x, a1y, a1z = -gma_LL * miu * (My * Hz - Mz * Hy), -gma_LL * miu * (Mz * Hx - Mx * Hz), -gma_LL * miu * (Mx * Hy - My * Hx) - (LL_lamda * miu / M0) * (My * (My * Hz - Mz * Hy) - Mz * (Mz * Hx - Mx * Hz))
        b1x, b1y, b1z = -gma_LL * miu * ((My + dt/2 * a1y) * Hz - (Mz + dt/2 * a1z) * Hy), -gma_LL * miu * ((Mz + dt/2 * a1z) * Hx - (Mx + dt/2 * a1x) * Hz), -gma_LL * miu * ((Mx + dt/2 * a1x) * Hy - (My + dt/2 * a1y) * Hx) - (LL_lamda * miu / M0) * ((My + dt/2 * a1y) * ((My + dt/2 * a1y) * Hz - (Mz + dt/2 * a1z) * Hy) - (Mz + dt/2 * a1z) * ((Mz + dt/2 * a1z) * Hx - (Mx + dt/2 * a1x) * Hz))
        c1x, c1y, c1z = -gma_LL * miu * ((My + dt/2 * b1y) * Hz - (Mz + dt/2 * b1z) * Hy), -gma_LL * miu * ((Mz + dt/2 * b1z) * Hx - (Mx + dt/2 * b1x) * Hz), -gma_LL * miu * ((Mx + dt/2 * b1x) * Hy - (My + dt/2 * b1y) * Hx) - (LL_lamda * miu / M0) * ((My + dt/2 * b1y) * ((My + dt/2 * b1y) * Hz - (Mz + dt/2 * b1z) * Hy) - (Mz + dt/2 * b1z) * ((Mz + dt/2 * b1z) * Hx - (Mx + dt/2 * b1x) * Hz))
        d1x, d1y, d1z = -gma_LL * miu * ((My + dt * c1y) * Hz - (Mz + dt * c1z) * Hy), -gma_LL * miu * ((Mz + dt * c1z) * Hx - (Mx + dt * c1x) * Hz), -gma_LL * miu * ((Mx + dt * c1x) * Hy - (My + dt * c1y) * Hx) - (LL_lamda * miu / M0) * ((My + dt * c1y) * ((My + dt * c1y) * Hz - (Mz + dt * c1z) * Hy) - (Mz + dt * c1z) * ((Mz + dt * c1z) * Hx - (Mx + dt * c1x) * Hz))
        
        # Update magnetization using RK4 integration
        result[idx, 0] = Mx + dt/6 * (a1x + 2 * b1x + 2 * c1x + d1x)
        result[idx, 1] = My + dt/6 * (a1y + 2 * b1y + 2 * c1y + d1y)
        result[idx, 2] = Mz + dt/6 * (a1z + 2 * b1z + 2 * c1z + d1z)





@cuda.jit
def update_H_kernel(H, E, lenz):
    """
    CUDA kernel to update H using E and lenz.
    
    Parameters:
    -----------
    H : numpy.ndarray
        Array of H vectors to be updated (2D array).
    E : numpy.ndarray
        Array of E vectors used for update (2D array).
    lenz : int
        Length of H and E vectors.
    """
    idx = cuda.grid(1)
    
    if idx < lenz - 1:
        # Update H(1, 2:lenz, n+1)
        H[0, idx + 1] = H[0, idx + 1] + 0.5 * (E[1, idx + 1] - E[1, idx])
        
        # Update H(2, 2:lenz, n+1)
        H[1, idx + 1] = H[1, idx + 1] - 0.5 * (E[0, idx + 1] - E[0, idx])



@cuda.jit
def update_E_kernel(H, E, lenz):
    """
    CUDA kernel to update H using E and lenz.
    
    Parameters:
    -----------
    H : numpy.ndarray
        Array of H vectors to be updated (2D array).
    E : numpy.ndarray
        Array of E vectors used for update (2D array).
    lenz : int
        Length of H and E vectors.
    """
    idx = cuda.grid(1)
    
    if idx < lenz - 1:
        # Update H(1, 2:lenz, n+1)
        H[0, idx + 1] = H[0, idx + 1] + 0.5 * (E[1, idx + 1] - E[1, idx])
        
        # Update H(2, 2:lenz, n+1)
        H[1, idx + 1] = H[1, idx + 1] - 0.5 * (E[0, idx + 1] - E[0, idx])


# Host (CPU) code
def LLG_step(M : np.array , H : np.array, dt : float, alpha : float):
    # Number of 3D vectors (N)
    N = M.size(2)
    

    # Allocate memory on the GPU
    M_device = cuda.to_device(M)
    H_device = cuda.to_device(H)
    result_device = cuda.device_array_like(M)  # Allocate output array on the GPU
    
    # Configure the kernel
    threads_per_block = 32
    blocks_per_grid = (N + threads_per_block - 1) // threads_per_block
    
    # Launch the kernel
    runge_kutta_LLG_kernel[blocks_per_grid, threads_per_block](M_device, H_device, dt, alpha, result_device)
    
    # Copy the result back to the host
    return result_device.copy_to_host()


if __name__ == '__main__':
    LLG_step()