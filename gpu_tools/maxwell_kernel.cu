#include <cuda_runtime.h>
#include <cuComplex.h>

// Precision: fp64 by default, compile with -DUSE_FP32 for single precision
#ifdef USE_FP32
typedef float real;
typedef float2 real2;
#define make_real2 make_float2
#else
typedef double real;
typedef double2 real2;
#define make_real2 make_double2
#endif

#define WARPS_PER_BLOCK 32
#define ZERO make_real2((real)0, (real)0)
#define PULSE_CELL 3

// HOST_BOUNDARY (temporary, for benchmarking): pulse and ABC run as host CuPy ops instead

// ============================================================================
// LLG - Landau-Lifshitz-Gilbert step (RK4)
// ============================================================================

// Helper: Cross Product (Register-only)
__device__ void cross_product(
    real ax, real ay, real az,
    real bx, real by, real bz,
    real &rx, real &ry, real &rz 
) {
    rx = ay * bz - az * by;
    ry = az * bx - ax * bz;
    rz = ax * by - ay * bx;
}

// Helper: LLG Derivative Logic
__device__ void compute_llg_derivative(
    real Mx, real My, real Mz,      
    real Hx_real, real Hy_real, real Hz_real,
    real neg_gamma_LL, real neg_coeff_damp,  // Constants
    real &kx, real &ky, real &kz     // Output 
) {
    // M x H
    real c1x, c1y, c1z;
    cross_product(Mx, My, Mz, Hx_real, Hy_real, Hz_real, c1x, c1y, c1z);

    // M x (M x H)
    real c2x, c2y, c2z;
    cross_product(Mx, My, Mz, c1x, c1y, c1z, c2x, c2y, c2z);

    kx = neg_gamma_LL * c1x + neg_coeff_damp * c2x;
    ky = neg_gamma_LL * c1y + neg_coeff_damp * c2y;
    kz = neg_gamma_LL * c1z + neg_coeff_damp * c2z;
}

// Main Device LLG Function
__device__ void LLG_RK4_calculation(
    real Mx, real My, real Mz,
    real Hx_real, real Hy_real, real Hz_real,
    real dt, real neg_gamma_LL, real neg_coeff_damp,
    real &M_next_x, real &M_next_y, real &M_next_z
) 
{
    // RK4 Accumulators
    real kx, ky, kz;
    real sum_x = 0, sum_y = 0, sum_z = 0;
    real half_dt = (real)0.5 * dt;

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
    sum_x += (real)2*kx; sum_y += (real)2*ky; sum_z += (real)2*kz;

    // Step 3
    compute_llg_derivative(
        Mx + half_dt*kx, My + half_dt*ky, Mz + half_dt*kz, 
        Hx_real, Hy_real, Hz_real, 
        neg_gamma_LL, neg_coeff_damp, 
        kx, ky, kz
    );
    sum_x += (real)2*kx; sum_y += (real)2*ky; sum_z += (real)2*kz;

    // Step 4
    compute_llg_derivative(
        Mx + dt*kx, My + dt*ky, Mz + dt*kz, 
        Hx_real, Hy_real, Hz_real, 
        neg_gamma_LL, neg_coeff_damp, 
        kx, ky, kz
    );
    sum_x += kx; sum_y += ky; sum_z += kz;

    // Final Update
    real dt_div_6 = dt / (real)6;
    M_next_x = Mx + dt_div_6 * sum_x;
    M_next_y = My + dt_div_6 * sum_y;
    M_next_z = Mz + dt_div_6 * sum_z;
}

// ============================================================================
// Yee stencils - neighbour cell values via warp shuffle / shared memory
// ============================================================================

// Helper: Manage memory performance to get the left neighbor cell value
__device__ void get_left_stencil(
    real2 curr, real2& left, real2* warp_edges, 
    real2* global_data, int idx
) {
    const int lane_id = threadIdx.x % warpSize;
    const int warp_id = threadIdx.x / warpSize;

    // Inside a warp: shuffle
    left.x = __shfl_up_sync(0xffffffff, curr.x, 1);
    left.y = __shfl_up_sync(0xffffffff, curr.y, 1);

    // Export right boundary to Shared MEM
    if (lane_id == warpSize - 1) warp_edges[warp_id] = curr;
    __syncthreads(); // prevent race condition

    // Warp boundaries logic
    if (lane_id == 0) {
        if (warp_id > 0) {
            left = warp_edges[warp_id - 1]; // From previous warp
        } else if (idx > 0) {
            left = global_data[idx - 1];    // Global boundary check
        } else {
            left = ZERO;                // Padding
        }
    }

    // Sync to validate reading before future calls
    __syncthreads();
}

// Helper: Manage memory performance to get the right neighbor cell value
__device__ void get_right_stencil(
    real2 curr, real2& right, real2* warp_edges, 
    real2* global_data, int N, int idx
) {
    const int lane_id = threadIdx.x % warpSize;
    const int warp_id = threadIdx.x / warpSize;

    // Inside a warp: shuffle
    right.x = __shfl_down_sync(0xffffffff, curr.x, 1);
    right.y = __shfl_down_sync(0xffffffff, curr.y, 1);

    // Export left boundary to Shared MEM
    if (lane_id == 0) warp_edges[warp_id] = curr;
    __syncthreads(); // prevent race condition

    // Warp boundaries logic
    if (lane_id == warpSize - 1) {
        if (warp_id < (blockDim.x / warpSize) - 1) {
            right = warp_edges[warp_id + 1]; // From next warp
        } else if (idx < N - 1) {
            right = global_data[idx + 1];    // Global boundary check
        } else {
            right = ZERO;                // Padding
        }
    }

    // Sync to validate reading before future calls
    __syncthreads();
}

// ============================================================================
// Kernels - one time step = magnetic_kernel then electric_kernel
// Launch grid (blocks_x, C): blockIdx.x splits the cells of one simulation,
// blockIdx.y = c selects the case. Case c's arrays start at c * (case size),
// per-case scalars are arrays indexed by c. C = 1 is a single run.
// ============================================================================

// Main Magnetic Kernel - ABC & Faraday & LLG & Coupling
extern "C" __global__ void magnetic_kernel(
    real2* E, real2* H,                     // Per case: 2xN / 3xN Complex Matrices
    real* M,                                // Per case: 3x(material_size) Real Matrix
    int N, real coeff,                      // Maxwell coefficient
    int material_start, int material_end,   // Material spatial bounds
    real dt,
    const real* neg_gamma_LL, const real* neg_lambda_mu, // Per-case LLG constants
    const int* closed_system,               // Per case, 1: M feeds back into H
    real2* abc_left, real2* abc_right,      // Per case ABC buffers [2 slots][Ex, Ey]
    const real2* pulse                      // This step's source [Ex, Ey], shared by all cases
) {
    // Shared boundary arrays
    __shared__ real2 warp_edges[WARPS_PER_BLOCK];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    bool valid = (idx < N);

    // Case offsets
    const int c = blockIdx.y;
    int material_size = material_end - material_start;
    E += c * 2 * N;
    H += c * 3 * N;
    M += c * 3 * material_size;
    abc_left += c * 4;
    abc_right += c * 4;

    // Base Pointers
    real2* Ex = E;
    real2* Ey = E + N;
    real2* Hx = H;
    real2* Hy = H + N;
    real2* Hz = H + 2 * N;

    // GET H FIELD
    real2 Hx_curr = valid ? Hx[idx] : ZERO;
    real2 Hy_curr = valid ? Hy[idx] : ZERO;
    real2 Hz_curr = valid ? Hz[idx] : ZERO;

    // GET E FIELD
    real2 Ex_curr = valid ? Ex[idx] : ZERO;
    real2 Ey_curr = valid ? Ey[idx] : ZERO;

#ifndef HOST_BOUNDARY
    // Soft source: add the CPU-precomputed pulse to the source cell
    if (idx == PULSE_CELL) {
        real2 p_x = pulse[0], p_y = pulse[1];
        Ex_curr.x += p_x.x; Ex_curr.y += p_x.y;
        Ey_curr.x += p_y.x; Ey_curr.y += p_y.y;
        Ex[idx] = Ex_curr;
        Ey[idx] = Ey_curr;
    }

    // ABC boundary condition [Sullivan 1.3 pg. 4]: E^n(0) = E^(n-2)(1) at Courant 0.5.
    // Only the owner thread of each edge cell touches its buffer, so any block count is safe.
    if (idx == 0) {
        Ex_curr = abc_left[0];
        Ey_curr = abc_left[1];
        abc_left[0] = abc_left[2];
        abc_left[1] = abc_left[3];
        Ex[idx] = Ex_curr;
        Ey[idx] = Ey_curr;
    } else if (idx == N - 1) {
        Ex_curr = abc_right[0];
        Ey_curr = abc_right[1];
        abc_right[0] = abc_right[2];
        abc_right[1] = abc_right[3];
        Ex[idx] = Ex_curr;
        Ey[idx] = Ey_curr;
    }
#endif

    real2 Ex_left, Ey_left; 
    
    get_left_stencil(Ex_curr, Ex_left, warp_edges, Ex, idx);
    get_left_stencil(Ey_curr, Ey_left, warp_edges, Ey, idx);
    if (!valid) return;

    // GET M FIELD
    real *Mx = M + idx - material_start;
    real *My = Mx + material_size;
    real *Mz = My + material_size;
    bool in_material = (idx >= material_start && idx < material_end);
    real Mx_curr, My_curr, Mz_curr;
    real Mx_next, My_next, Mz_next;

    // LLG interaction between H and M, driven by H^n (before Faraday)
    if (in_material) {
        Mx_curr = *Mx; My_curr = *My; Mz_curr = *Mz;
        real M_mag = sqrt(Mx_curr*Mx_curr + My_curr*My_curr + Mz_curr*Mz_curr);
        real neg_coeff_damp = neg_lambda_mu[c] / M_mag;
        LLG_RK4_calculation(
            Mx_curr, My_curr, Mz_curr, 
            Hx_curr.x, Hy_curr.x, Hz_curr.x, 
            dt, neg_gamma_LL[c], neg_coeff_damp, 
            Mx_next, My_next, Mz_next
        );
    }

    // Faraday's law FDTD - Update H [Sullivan eq. 1.9b], H[i] sits at i-1/2
    // H[0] has no left neighbour and stays untouched
    if (idx > 0) {
        Hx_curr.x += coeff * (Ey_curr.x - Ey_left.x); // Hx affected by Ey
        Hx_curr.y += coeff * (Ey_curr.y - Ey_left.y);
        Hy_curr.x -= coeff * (Ex_curr.x - Ex_left.x); // Hy affected by Ex
        Hy_curr.y -= coeff * (Ex_curr.y - Ex_left.y);
    }

    if (in_material) {
        // Coupling - add M to H
        if (closed_system[c]) {
            Hx_curr.x += Mx_curr - Mx_next;
            Hy_curr.x += My_curr - My_next;
            Hz_curr.x += Mz_curr - Mz_next;
        }

        // Write back values
        *Mx = Mx_next;
        *My = My_next;
        *Mz = Mz_next;
    }
    Hx[idx] = Hx_curr;
    Hy[idx] = Hy_curr;
    Hz[idx] = Hz_curr;
}

// Main Electric Kernel - Ampere & ABC buffer fill
extern "C" __global__ void electric_kernel(
    real2* E, real2* H, int N,
    const real* coeff_E, const real* coeff_H, // Per case, per-cell ca, cb [Sullivan eq. 3.8]
    real2* abc_left, real2* abc_right
) {
    // Shared boundary arrays
    __shared__ real2 warp_edges[WARPS_PER_BLOCK];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    bool valid = (idx < N);

    // Case offsets
    const int c = blockIdx.y;
    E += c * 2 * N;
    H += c * 3 * N;
    coeff_E += c * N;
    coeff_H += c * N;
    abc_left += c * 4;
    abc_right += c * 4;

    // Base Pointers
    real2* Ex = E;
    real2* Ey = E + N;
    real2* Hx = H;
    real2* Hy = H + N;

    // GET E FIELD
    real2* Ex_ptr = Ex + idx;
    real2* Ey_ptr = Ey + idx;
    real2 Ex_curr = valid ? *Ex_ptr : ZERO;
    real2 Ey_curr = valid ? *Ey_ptr : ZERO;

    // GET H FIELD
    real2 Hx_curr = valid ? Hx[idx] : ZERO;
    real2 Hy_curr = valid ? Hy[idx] : ZERO;
    real2 Hx_right, Hy_right; 
    
    get_right_stencil(Hx_curr, Hx_right, warp_edges, Hx, N, idx);
    get_right_stencil(Hy_curr, Hy_right, warp_edges, Hy, N, idx);

#ifndef HOST_BOUNDARY
    // ABC: push E^n of the cells next to the edges into the second buffer slot
    if (idx == 1) {
        abc_left[2] = Ex_curr;
        abc_left[3] = Ey_curr;
    } else if (idx == N - 2) {
        abc_right[2] = Ex_curr;
        abc_right[3] = Ey_curr;
    }
#endif

    // Ampere's law FDTD - Update E [Sullivan eq. 1.9a]
    // E[N-1] has no right neighbour and is set by the ABC
    if (idx < N - 1) {
        real ca = coeff_E[idx], cb = coeff_H[idx];
        Ex_curr.x = ca * Ex_curr.x - cb * (Hy_right.x - Hy_curr.x); // Ex affected by Hy
        Ex_curr.y = ca * Ex_curr.y - cb * (Hy_right.y - Hy_curr.y);
        Ey_curr.x = ca * Ey_curr.x + cb * (Hx_right.x - Hx_curr.x); // Ey affected by Hx
        Ey_curr.y = ca * Ey_curr.y + cb * (Hx_right.y - Hx_curr.y);
    }

    // Write back values
    if (valid) {
        *Ex_ptr = Ex_curr;
        *Ey_ptr = Ey_curr;
    }
}
