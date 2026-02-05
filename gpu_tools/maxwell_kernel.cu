#include <cuda_runtime.h>
#define WARPS_PER_BLOCK 32


__device__ __forceinline__ 
void get_stencil(float local_val, float& left, float& right, 
                    float* left_edges, float* right_edges, 
                    float* global_data, int N) {
    const int lane_id = threadIdx.x % warpSize;
    const int warp_id = threadIdx.x / warpSize;
    const int global_idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Inside warp shuffle
    left = __shfl_up_sync(0xffffffff, local_val, 1);
    right = __shfl_down_sync(0xffffffff, local_val, 1);

    // Export boundaries to Shared MEM
    if (lane_id == 0)            left_edges[warp_id]  = local_val;
    if (lane_id == warpSize - 1) right_edges[warp_id] = local_val;
    __syncthreads();

    // Warp boundaries logic
    if (lane_id == 0) {
        if (warp_id > 0) {
            left = right_edges[warp_id - 1]; // From previous warp
        } else {
            left = (global_idx > 0) ? global_data[global_idx - 1] : 0.0f;
        }
    }

    if (lane_id == warpSize - 1) {
        if (warp_id < (blockDim.x / warpSize) - 1) {
            right = left_edges[warp_id + 1]; // From next warp
        } else {
            right = (global_idx < N - 1) ? global_data[global_idx + 1] : 0.0f;
        }
    }

    // Sync for reading validation before next calls
    __syncthreads();
}

__global__ void maxwell_kernel(float* E, float* H, int N, float coeff) {
    // Shared boundary arrays
    __shared__ float left_edges[WARPS_PER_BLOCK];
    __shared__ float right_edges[WARPS_PER_BLOCK];

    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (global_idx >= N) return;

    // GET FIELDS
    float E_curr = E[global_idx];
    float E_left, E_right; // Register neighbors
    get_stencil(E_curr, E_left, E_right, left_edges, right_edges, E, N);

    // Physics
    H[global_idx] -= coeff * (E_right - E_left);
    // ...
}