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
    float Hx_real, float Hy_real, float Hz_real,
    float neg_gamma_LL, float neg_coeff_damp,  // Constants
    float &kx, float &ky, float &kz)     // Output 
{
    // M x H
    float c1x, c1y, c1z;
    cross_product(Mx, My, Mz, Hx_real, Hy_real, Hz_real, c1x, c1y, c1z);

    // M x (M x H)
    float c2x, c2y, c2z;
    cross_product(Mx, My, Mz, c1x, c1y, c1z, c2x, c2y, c2z);

    kx = neg_gamma_LL * c1x + neg_coeff_damp * c2x;
    ky = neg_gamma_LL * c1y + neg_coeff_damp * c2y;
    kz = neg_gamma_LL * c1z + neg_coeff_damp * c2z;
}

// Main Device Function
__device__ __forceinline__ void LLG_RK4_calculation(
    float Mx, float My, float Mz,
    float Hx_real, float Hy_real, float Hz_real,
    float dt, float neg_gamma_LL, float neg_coeff_damp,
    float &M_next_x, float &M_next_y, float &M_next_z
) 
{
    // RK4 Accumulators
    float kx, ky, kz;
    float sum_x = 0.0f, sum_y = 0.0f, sum_z = 0.0f;
    float half_dt = 0.5f * dt;

    // Step 1
    compute_llg_derivative(
        Mx, My, Mz, 
        Hx_real, Hy_real, Hz_real, 
        neg_gamma_LL, neg_coeff_damp, 
        kx, ky, kz
    );
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Step 2
    compute_llg_derivative(
        Mx + half_dt*kx, My + half_dt*ky, Mz + half_dt*kz, 
        Hx_real, Hy_real, Hz_real, 
        neg_gamma_LL, neg_coeff_damp, 
        kx, ky, kz
    ); 
    sum_x += 2.0f*kx; sum_y += 2.0f*ky; sum_z += 2.0f*kz;

    // Step 3
    compute_llg_derivative(
        Mx + half_dt*kx, My + half_dt*ky, Mz + half_dt*kz, 
        Hx_real, Hy_real, Hz_real, 
        neg_gamma_LL, neg_coeff_damp, 
        kx, ky, kz
    );
    sum_x += 2.0f*kx; sum_y += 2.0f*ky; sum_z += 2.0f*kz;

    // Step 4
    compute_llg_derivative(
        Mx + dt*kx, My + dt*ky, Mz + dt*kz, 
        Hx_real, Hy_real, Hz_real, 
        neg_gamma_LL, neg_coeff_damp, 
        kx, ky, kz
    );
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Final Update
    float dt_div_6 = dt / 6.0f;
    M_next_x = Mx + dt_div_6 * sum_x;
    M_next_y = My + dt_div_6 * sum_y;
    M_next_z = Mz + dt_div_6 * sum_z;
}
