__global__ void maxwell_kernel(float* global_E, float* shared_boundary, float* output_H) {
    int lane_id = threadIdx.x % 32;
    int warp_id = threadIdx.x / 32;
    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;

    float E_local = global_E[global_idx];
    float E_left, E_right;

    E_left = __shfl_up_sync(0xffffffff, E_local, 1); // get E from index-1 neighbor

    // from (index + 1)
    E_right = __shfl_down_sync(0xffffffff, E_local, 1);

    // Boundary Conditions
    if (lane_id == 0) {
        E_left = (global_idx > 0) ? global_E[global_idx - 1] : 0.0f; 
    }
    
    // thread 31
    if (lane_id == 31) {
        E_right = global_E[global_idx + 1]; 
    }

    // physics
    output_H[global_idx] = (E_right - E_left) * 0.5f; 
}