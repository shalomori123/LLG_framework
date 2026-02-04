#include <cuda_runtime.h>
#define WARPS_PER_BLOCK 32

__global__ void maxwell_kernel(float* global_E, float* output_H, int N) {
    __shared__ float left_edges[WARPS_PER_BLOCK];
    __shared__ float right_edges[WARPS_PER_BLOCK];

    int lane_id = threadIdx.x % warpSize;
    int warp_id = threadIdx.x / warpSize;
    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (global_idx >= N) return;

    float E_local = global_E[global_idx];

    float E_left = __shfl_up_sync(0xffffffff, E_local, 1);
    float E_right = __shfl_down_sync(0xffffffff, E_local, 1);

    // share edges to Shared Memory
    if (lane_id == 0)            left_edges[warp_id]  = E_local;
    if (lane_id == warpSize - 1) right_edges[warp_id] = E_local;
    __syncthreads();

    // boundaries
    if (lane_id == 0) {
        if (warp_id > 0) E_left = right_edges[warp_id - 1];
        else  E_left = (global_idx > 0) ? global_E[global_idx - 1] : 0;
    }

    if (lane_id == warpSize - 1) {
        if (warp_id < (blockDim.x / warpSize) - 1) E_right = left_edges[warp_id + 1];
        else   E_right = global_E[global_idx + 1];
    }

    // physics
    output_H[global_idx] = (E_right - E_left) * 0.5f; 
}