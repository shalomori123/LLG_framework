#include <cuda_runtime.h>
#include <cuComplex.h>

#include "llg_kernel.cu"

#define WARPS_PER_BLOCK 32
#define ZERO make_float2(0.0f, 0.0f)
#define CLOSED_SYSTEM 1

__device__ __forceinline__ 
void get_left_stencil(float2 curr, float2& left, float2* warp_edges, 
                         float2* global_data, int idx) {
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

__device__ __forceinline__ 
void get_right_stencil(float2 curr, float2& right, float2* warp_edges, 
                          float2* global_data, int N, int idx) {
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

__global__ void magnetic_kernel(
    float2* E, float2* H, 
    float* M_curr, float* M_next,
    int N, float coeff,                     // Maxwell params
    int material_start, int material_end,   // LLG spatial bounds
    float dt, float neg_gamma_LL, float neg_coeff_damp // LLG constants
) {
    // Shared boundary arrays
    __shared__ float2 warp_edges[WARPS_PER_BLOCK];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    bool valid = (idx < N);
    int y_idx = N + idx;

    // GET H FIELD
    float2* Hx_curr = H + idx;
    float2* Hy_curr = H + y_idx;

    // calculate E field (Safe fetch for out-of-bounds threads)
    float2 Ex_curr = valid ? E[idx] : ZERO;
    float2 Ey_curr = valid ? E[y_idx] : ZERO;
    float2 Ex_right, Ey_right; // Register neighbors
    
    get_right_stencil(Ex_curr, Ex_right, warp_edges, E, N, idx);
    get_right_stencil(Ey_curr, Ey_right, warp_edges, E, 2*N, y_idx);

    if (valid) {
        // Faraday's law
        Hx_curr->x += coeff * (Ey_right.x - Ey_curr.x); // Hx affected by Ey
        Hx_curr->y += coeff * (Ey_right.y - Ey_curr.y);
        Hy_curr->x += coeff * (Ex_right.x - Ex_curr.x); // Hy affected by Ex
        Hy_curr->y += coeff * (Ex_right.y - Ex_curr.y);

        // LLG interaction between H and M
        LLG_RK4_calculation(H, M_curr, M_next, idx, 
                            material_start, material_end, N, 
                            dt, neg_gamma_LL, neg_coeff_damp);

        // Coupling - add M to H
        #ifdef CLOSED_SYSTEM
            // H_next[material_begin:material_end, :] += M_current - M_next
        #endif
    }
}

__global__ void electric_kernel(float2* E, float2* H, int N, float coeff,
                                float2* abc_left, float2* abc_right) {
    // Shared boundary arrays
    __shared__ float2 warp_edges[WARPS_PER_BLOCK];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    bool valid = (idx < N);
    int y_idx = N + idx;

    // GET E FIELD
    float2* Ex_curr = E + idx;
    float2* Ey_curr = E + y_idx;

    // calculate H field (Safe fetch for out-of-bounds threads)
    float2 Hx_curr = valid ? H[idx] : ZERO;
    float2 Hy_curr = valid ? H[y_idx] : ZERO;
    float2 Hx_left, Hy_left; 
    
    get_left_stencil(Hx_curr, Hx_left, warp_edges, H, idx);
    get_left_stencil(Hy_curr, Hy_left, warp_edges, H, y_idx);

    if (valid) {
        // Ampere's law
        Ex_curr->x += coeff * (Hy_left.x - Hy_curr.x); // Ex affected by Hy
        Ex_curr->y += coeff * (Hy_left.y - Hy_curr.y);
        Ey_curr->x += coeff * (Hx_left.x - Hx_curr.x); // Ey affected by Hx
        Ey_curr->y += coeff * (Hx_left.y - Hx_curr.y);
    }

    // ABC boundry condition logic [Sullivan 1.3 pg. 4]
    // Read buffer value and then write the next.
    if (idx == 0) {
        *Ex_curr = abc_left[0];
        *Ey_curr = abc_left[1];
    } else if (idx == N - 1) {
        *Ex_curr = abc_right[0];
        *Ey_curr = abc_right[1];
    }
    __syncthreads();
    if (idx == 1) {
        abc_left[0] = *Ex_curr;
        abc_left[1] = *Ey_curr;
    } else if (idx == N - 2) {
        abc_right[0] = *Ex_curr;
        abc_right[1] = *Ey_curr;
    }
}