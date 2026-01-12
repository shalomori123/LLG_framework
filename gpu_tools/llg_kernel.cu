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
    double neg_gamma_LL, double neg_coeff_damp,  // Constants
    double &kx, double &ky, double &kz)     // Output 
{
    // M x H
    double c1x, c1y, c1z;
    cross_product(Mx, My, Mz, Hx, Hy, Hz, c1x, c1y, c1z);

    // M x (M x H)
    double c2x, c2y, c2z;
    cross_product(Mx, My, Mz, c1x, c1y, c1z, c2x, c2y, c2z);

    kx = neg_gamma_LL * c1x + neg_coeff_damp * c2x;
    ky = neg_gamma_LL * c1y + neg_coeff_damp * c2y;
    kz = neg_gamma_LL * c1z + neg_coeff_damp * c2z;
}

// Main Kernel
__global__ void LLG_RK4_kernel(
    const double* __restrict__ M_curr,  // Real
    const double* __restrict__ H_curr,  // Complex, only real part (double) needed
    double* __restrict__ M_next,       
    int material_size,                            
    double dt,                         
    double neg_gamma_LL,  // -gamma_LL (negative)
    double neg_coeff_damp // -(gamma * alpha * factor) / M0
) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= material_size) return;

    // Load initial state to registers
    int ptr_m = idx * 3; 
    double Mx = M_curr[ptr_m];
    double My = M_curr[ptr_m + 1];
    double Mz = M_curr[ptr_m + 2];

    int ptr_h = idx * 6; // complex
    double Hx = H_curr[ptr_h];
    double Hy = H_curr[ptr_h + 2]; // next real part
    double Hz = H_curr[ptr_h + 4];

    // RK4 Accumulators
    double kx, ky, kz;
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double half_dt = 0.5*dt;

    // Step 1
    compute_llg_derivative(Mx, My, Mz, Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Step 2
    compute_llg_derivative(Mx + half_dt*kx, My + half_dt*ky, Mz + half_dt*kz, 
                           Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz); 
    sum_x += 2.0*kx; sum_y += 2.0*ky; sum_z += 2.0*kz;

    // Step 3
    compute_llg_derivative(Mx + half_dt*kx, My + half_dt*ky, Mz + half_dt*kz, 
                           Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz);
    sum_x += 2.0*kx; sum_y += 2.0*ky; sum_z += 2.0*kz;

    // Step 4
    compute_llg_derivative(Mx + dt*kx, My + dt*ky, Mz + dt*kz, 
                           Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Final Update
    double dt_div_6 = dt / 6.0;
    M_next[ptr_m]     = Mx + dt_div_6 * sum_x;
    M_next[ptr_m + 1] = My + dt_div_6 * sum_y;
    M_next[ptr_m + 2] = Mz + dt_div_6 * sum_z;
}
