#include <cuda_runtime.h>

// Helper: Cross Product (Register-only)
__device__ __forceinline__ void cross_product(
    float ax, float ay, float az,
    float bx, float by, float bz,
    float &rx, float &ry, float &rz) 
{
    rx = ay * bz - az * by;
    ry = az * bx - ax * bz;
    rz = ax * by - ay * bx;
}

// Helper: LLG Derivative Logic
__device__ __forceinline__ void compute_llg_derivative(
    float Mx, float My, float Mz,
    float Hx, float Hy, float Hz,
    float m_gamma_LL, float coeff_damp,  // Constants
    float &kx, float &ky, float &kz)     // Output
{
    // M x H
    float c1x, c1y, c1z;
    cross_product(Mx, My, Mz, Hx, Hy, Hz, c1x, c1y, c1z);

    // M x (M x H)
    float c2x, c2y, c2z;
    cross_product(Mx, My, Mz, c1x, c1y, c1z, c2x, c2y, c2z);

    kx = m_gamma_LL * c1x - coeff_damp * c2x;
    ky = m_gamma_LL * c1y - coeff_damp * c2y;
    kz = m_gamma_LL * c1z - coeff_damp * c2z;
}

// Main Kernel
__global__ void LLG_RK4_kernel(
    const float* __restrict__ M_curr, 
    const float* __restrict__ H_curr, 
    float* __restrict__ M_next,
    int N,
    float dt,
    float alpha,
    float gamma,
    float M0,
    int stride
) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // Pre-calculate constants
    float factor = 1.0f / (1.0f + alpha * alpha);
    float m_gamma_LL = -gamma * factor;  // -gamma_LL (negative)
    float coeff_damp = (gamma * alpha * factor) / M0; 

    // Load initial state to registers
    int ptr = idx * stride; 
    
    float Mx = M_curr[ptr];
    float My = M_curr[ptr + 1];
    float Mz = M_curr[ptr + 2];

    float Hx = H_curr[ptr];
    float Hy = H_curr[ptr + 1];
    float Hz = H_curr[ptr + 2];

    // RK4 Accumulators
    float kx, ky, kz;
    float sum_x = 0.0f, sum_y = 0.0f, sum_z = 0.0f;

    // Step 1
    compute_llg_derivative(Mx, My, Mz, Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Step 2
    compute_llg_derivative(Mx + 0.5f*dt*kx, My + 0.5f*dt*ky, Mz + 0.5f*dt*kz, 
                           Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz); 
    sum_x += 2.0f*kx; sum_y += 2.0f*ky; sum_z += 2.0f*kz;

    // Step 3
    compute_llg_derivative(Mx + 0.5f*dt*kx, My + 0.5f*dt*ky, Mz + 0.5f*dt*kz, 
                           Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz);
    sum_x += 2.0f*kx; sum_y += 2.0f*ky; sum_z += 2.0f*kz;

    // Step 4
    compute_llg_derivative(Mx + dt*kx, My + dt*ky, Mz + dt*kz, 
                           Hx, Hy, Hz, m_gamma_LL, coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Final Update
    float dt_div_6 = dt / 6.0f;
    M_next[ptr]     = Mx + dt_div_6 * sum_x;
    M_next[ptr + 1] = My + dt_div_6 * sum_y;
    M_next[ptr + 2] = Mz + dt_div_6 * sum_z;
}
