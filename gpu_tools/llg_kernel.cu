#include <cuda_runtime.h>

// Helper: Cross Product (Register-only)
__device__ __forceinline__ void cross_product(
    double ax, double ay, double az,
    double bx, double by, double bz,
    double &rx, double &ry, double &rz) 
{
    rx = ay * bz - az * by;
    ry = az * bx - ax * bz;
    rz = ax * by - ay * bx;
}

// Helper: LLG Derivative Logic
__device__ __forceinline__ void compute_llg_derivative(
    double Mx, double My, double Mz,      
    double Hx, double Hy, double Hz,      
    double m_gamma_LL, double coeff_damp,  // Constants
    double &kx, double &ky, double &kz)     // Output 
{
    // M x H
    double c1x, c1y, c1z;
    cross_product(Mx, My, Mz, Hx, Hy, Hz, c1x, c1y, c1z);

    // M x (M x H)
    double c2x, c2y, c2z;
    cross_product(Mx, My, Mz, c1x, c1y, c1z, c2x, c2y, c2z);

    kx = m_gamma_LL * c1x - coeff_damp * c2x;
    ky = m_gamma_LL * c1y - coeff_damp * c2y;
    kz = m_gamma_LL * c1z - coeff_damp * c2z;
}

// Main Kernel
__global__ void LLG_RK4_kernel(
    const double* __restrict__ M_curr, 
    const double* __restrict__ H_curr, 
    double* __restrict__ M_next,       
    int N,                            
    double dt,                         
    double alpha,                      
    double gamma,                      
    double M0,
    int stride 
) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // Constants
    double factor = 1.0 / (1.0 + alpha * alpha);
    double m_gamma_LL = -gamma * factor;  // -gamma_LL (negative)
    double coeff_damp = (gamma * alpha * factor) / M0; 

    // Load initial state to registers
    int ptr = idx * stride; 
    
    double Mx = M_curr[ptr];
    double My = M_curr[ptr + 1];
    double Mz = M_curr[ptr + 2];

    double Hx = H_curr[ptr];
    double Hy = H_curr[ptr + 1];
    double Hz = H_curr[ptr + 2];

    // RK4 Accumulators
    double kx, ky, kz;
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;

    // Step 1
    compute_llg_derivative(Mx, My, Mz, Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Step 2
    compute_llg_derivative(Mx + 0.5*dt*kx, My + 0.5*dt*ky, Mz + 0.5*dt*kz, 
                           Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz); 
    sum_x += 2.0*kx; sum_y += 2.0*ky; sum_z += 2.0*kz;

    // Step 3
    compute_llg_derivative(Mx + 0.5*dt*kx, My + 0.5*dt*ky, Mz + 0.5*dt*kz, 
                           Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz);
    sum_x += 2.0*kx; sum_y += 2.0*ky; sum_z += 2.0*kz;

    // Step 4
    compute_llg_derivative(Mx + dt*kx, My + dt*ky, Mz + dt*kz, 
                           Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Final Update
    double dt_div_6 = dt / 6.0;
    M_next[ptr]     = Mx + dt_div_6 * sum_x;
    M_next[ptr + 1] = My + dt_div_6 * sum_y;
    M_next[ptr + 2] = Mz + dt_div_6 * sum_z;
}