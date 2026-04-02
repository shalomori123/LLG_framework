#include <cuda_runtime.h>

#define WARPS_PER_BLOCK 32
#define BOUND 0.0f // TODO: ABC logic to the bounderies

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

__global__ void maxwell_kernel(float* E, float* H, int N, float coeff) {
    // Shared boundary arrays
    __shared__ float warp_edges[WARPS_PER_BLOCK];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // GET FIELDS
    float E_curr = E[idx];
    float H_curr = H[idx];
    float H_left, E_right; // Register neighbors
    get_left_stencil(H_curr, H_left, warp_edges, H, idx);
    get_right_stencil(E_curr, E_right, warp_edges, E, N, idx);

    // Physics
    E[idx] += coeff * (H_left - H_curr);
    H[idx] += coeff * (E_right - E_curr);
    // ...
}