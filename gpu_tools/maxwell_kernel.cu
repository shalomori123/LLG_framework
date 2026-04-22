#include <cuda_runtime.h>

#define WARPS_PER_BLOCK 32
#define BOUND 0.0f // TODO: ABC logic to the bounderies
#define CLOSED_SYSTEM 1

__device__ __forceinline__ 
void get_left_stencil(float curr, float& left, float* warp_edges, 
                         float* global_data, int idx) {
    const int lane_id = threadIdx.x % warpSize;
    const int warp_id = threadIdx.x / warpSize;

    // Inside a warp: shuffle
    left = __shfl_up_sync(0xffffffff, curr, 1);

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
            left = BOUND;                   // Padding
        }
    }

    // Sync to validate reading before future calls
    __syncthreads();
}

__device__ __forceinline__ 
void get_right_stencil(float curr, float& right, float* warp_edges, 
                          float* global_data, int N, int idx) {
    const int lane_id = threadIdx.x % warpSize;
    const int warp_id = threadIdx.x / warpSize;

    // Inside a warp: shuffle
    right = __shfl_down_sync(0xffffffff, curr, 1);

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
            right = BOUND;                   // Padding
        }
    }

    // Sync to validate reading before future calls
    __syncthreads();
}

__global__ void magnetic_kernel(float* E, float* H, int N, float coeff) {
    // Shared boundary arrays
    __shared__ float warp_edges[WARPS_PER_BLOCK];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // GET H FIELD
    float* Hx_curr = H + 3*idx;
    float* Hy_curr = H + 3*idx + 1;
    float* Hz_curr = H + 3*idx + 2;

    // calculate E field
    float Ex_curr = E[3*idx];
    float Ey_curr = E[3*idx + 1];
    float Ex_right, Ey_right; // Register neighbors
    get_right_stencil(Ex_curr, Ex_right, warp_edges, E, N, idx);
    get_right_stencil(Ey_curr, Ey_right, warp_edges, E, N, idx);

    // Faraday's law
    *Hx_curr += coeff * (Ey_right - Ey_curr); // Hx affected by Ey
    *Hy_curr += coeff * (Ex_right - Ex_curr); // Hy affected by Ex

    // LLG interaction between H and M
    // TODO: call LLG_RK4_kernel (which shouldn't be seperated kernel) and adapt parameters

    // Coupling - add M to H
    #ifdef CLOSED_SYSTEM
        // H_next[material_begin:material_end, :] += M_current - M_next
    #endif
}

__global__ void electric_kernel(float* E, float* H, int N, float coeff) {
    // Shared boundary arrays
    __shared__ float warp_edges[WARPS_PER_BLOCK];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // GET E FIELD
    float* Ex_curr = E + 3*idx;
    float* Ey_curr = E + 3*idx + 1;

    // calculate H field
    float Hx_curr = H[3*idx];
    float Hy_curr = H[3*idx + 1];
    float Hx_left, Hy_left; // Register neighbors
    get_left_stencil(Hx_curr, Hx_left, warp_edges, H, idx);
    get_left_stencil(Hy_curr, Hy_left, warp_edges, H, idx);

    // Faraday's law
    *Ex_curr += coeff * (Hy_left - Hy_curr); // Ex affected by Hy
    *Ey_curr += coeff * (Hx_left - Hx_curr); // Ey affected by Hx

    // TODO: add ABC logic
}