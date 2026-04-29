#pragma once
#include <cuda_runtime.h>
#include <cuComplex.h>

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
    float neg_gamma_LL, float neg_coeff_damp,  // Constants
    float &kx, float &ky, float &kz)     // Output 
{
    // M x H
    float c1x, c1y, c1z;
    cross_product(Mx, My, Mz, Hx, Hy, Hz, c1x, c1y, c1z);

    // M x (M x H)
    float c2x, c2y, c2z;
    cross_product(Mx, My, Mz, c1x, c1y, c1z, c2x, c2y, c2z);

    kx = neg_gamma_LL * c1x + neg_coeff_damp * c2x;
    ky = neg_gamma_LL * c1y + neg_coeff_damp * c2y;
    kz = neg_gamma_LL * c1z + neg_coeff_damp * c2z;
}

// Main Device Function
__device__ __forceinline__ void LLG_RK4_calculation(
    const float2* __restrict__ H_curr, 
    const float* __restrict__ M_curr,
    float* __restrict__ M_next,       
    int idx,                                     // Passed from the calling kernel
    int material_start, int material_end, int N,
    float dt,                         
    float neg_gamma_LL,  // -gamma_LL (negative)
    float neg_coeff_damp // -(gamma * alpha * factor) / M0
) 
{
    // Boundary check for the specific thread
    if (idx >= material_end || idx < material_start) return;

    int material_size = material_end - material_start;

    // Load initial state
    float Mx = M_curr[idx - material_start];
    float My = M_curr[material_size + idx - material_start];
    float Mz = M_curr[2 * material_size + idx - material_start];

    // Complex H field, M only interacts with the real part
    float Hx = H_curr[idx].x;
    float Hy = H_curr[N + idx].x;
    float Hz = H_curr[2 * N + idx].x;

    // RK4 Accumulators
    float kx, ky, kz;
    float sum_x = 0.0f, sum_y = 0.0f, sum_z = 0.0f;
    float half_dt = 0.5f * dt;

    // Step 1
    compute_llg_derivative(Mx, My, Mz, Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Step 2
    compute_llg_derivative(Mx + half_dt*kx, My + half_dt*ky, Mz + half_dt*kz, 
                           Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz); 
    sum_x += 2.0f*kx; sum_y += 2.0f*ky; sum_z += 2.0f*kz;

    // Step 3
    compute_llg_derivative(Mx + half_dt*kx, My + half_dt*ky, Mz + half_dt*kz, 
                           Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz);
    sum_x += 2.0f*kx; sum_y += 2.0f*ky; sum_z += 2.0f*kz;

    // Step 4
    compute_llg_derivative(Mx + dt*kx, My + dt*ky, Mz + dt*kz, 
                           Hx, Hy, Hz, neg_gamma_LL, neg_coeff_damp, kx, ky, kz);
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Final Update
    float dt_div_6 = dt / 6.0f;
    M_next[idx - material_start]                     = Mx + dt_div_6 * sum_x;
    M_next[material_size + idx - material_start]     = My + dt_div_6 * sum_y;
    M_next[2 * material_size + idx - material_start] = Mz + dt_div_6 * sum_z;
}
