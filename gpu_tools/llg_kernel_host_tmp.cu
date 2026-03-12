#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

#include <cuComplex.h>
#include <cuda_runtime.h>

#include "llg_kernel.cu"

namespace {

void check_cuda(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(status));
    }
}

void cross_product_host(
    double ax, double ay, double az,
    double bx, double by, double bz,
    double& rx, double& ry, double& rz) {
    rx = ay * bz - az * by;
    ry = az * bx - ax * bz;
    rz = ax * by - ay * bx;
}

void compute_llg_derivative_host(
    double mx, double my, double mz,
    double hx, double hy, double hz,
    double neg_gamma_ll, double neg_coeff_damp,
    double& kx, double& ky, double& kz) {
    double c1x, c1y, c1z;
    cross_product_host(mx, my, mz, hx, hy, hz, c1x, c1y, c1z);

    double c2x, c2y, c2z;
    cross_product_host(mx, my, mz, c1x, c1y, c1z, c2x, c2y, c2z);

    kx = neg_gamma_ll * c1x + neg_coeff_damp * c2x;
    ky = neg_gamma_ll * c1y + neg_coeff_damp * c2y;
    kz = neg_gamma_ll * c1z + neg_coeff_damp * c2z;
}

void llg_rk4_host(
    const std::vector<double>& m_curr,
    const std::vector<cuDoubleComplex>& h_curr,
    std::vector<double>& m_next,
    int material_size,
    double dt,
    double neg_gamma_ll,
    double neg_coeff_damp) {
    for (int idx = 0; idx < material_size; ++idx) {
        const int ptr_m = idx * 3;

        const double mx = m_curr[ptr_m];
        const double my = m_curr[ptr_m + 1];
        const double mz = m_curr[ptr_m + 2];

        const double hx = cuCreal(h_curr[idx]);
        const double hy = cuCreal(h_curr[material_size + idx]);
        const double hz = cuCreal(h_curr[2 * material_size + idx]);

        double k1x, k1y, k1z;
        double k2x, k2y, k2z;
        double k3x, k3y, k3z;
        double k4x, k4y, k4z;

        compute_llg_derivative_host(mx, my, mz, hx, hy, hz, neg_gamma_ll, neg_coeff_damp, k1x, k1y, k1z);
        compute_llg_derivative_host(
            mx + 0.5 * dt * k1x, my + 0.5 * dt * k1y, mz + 0.5 * dt * k1z,
            hx, hy, hz, neg_gamma_ll, neg_coeff_damp, k2x, k2y, k2z);
        compute_llg_derivative_host(
            mx + 0.5 * dt * k2x, my + 0.5 * dt * k2y, mz + 0.5 * dt * k2z,
            hx, hy, hz, neg_gamma_ll, neg_coeff_damp, k3x, k3y, k3z);
        compute_llg_derivative_host(
            mx + dt * k3x, my + dt * k3y, mz + dt * k3z,
            hx, hy, hz, neg_gamma_ll, neg_coeff_damp, k4x, k4y, k4z);

        m_next[ptr_m] = mx + (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        m_next[ptr_m + 1] = my + (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
        m_next[ptr_m + 2] = mz + (dt / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);
    }
}

__device__ __forceinline__ unsigned int get_smid() {
    unsigned int smid;
    asm volatile("mov.u32 %0, %smid;" : "=r"(smid));
    return smid;
}

__global__ void debug_thread_location_kernel(int* sm_ids, int material_size) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= material_size) return;

    const unsigned int smid = get_smid();
    sm_ids[idx] = static_cast<int>(smid);
    printf("GPU thread idx=%d block=%d thread=%d sm=%u\n", idx, blockIdx.x, threadIdx.x, smid);
}

}  // namespace

int main() {
    try {
        constexpr int material_size = 4;
        constexpr int threads_per_block = 128;
        const int blocks = (material_size + threads_per_block - 1) / threads_per_block;

        const double dt = 1e-15;
        const double neg_gamma_ll = -2.5e10;
        const double neg_coeff_damp = -4.0e8;

        std::vector<double> m_curr = {
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.2, 0.3, 0.9,
            -0.5, 0.4, 0.1,
        };

        // Component-major complex128 layout with shape (3, material_size).
        std::vector<cuDoubleComplex> h_curr = {
            make_cuDoubleComplex(0.0, 0.0),
            make_cuDoubleComplex(0.3, 0.0),
            make_cuDoubleComplex(0.8, 0.0),
            make_cuDoubleComplex(0.2, 0.0),
            make_cuDoubleComplex(0.5, 0.0),
            make_cuDoubleComplex(0.0, 0.0),
            make_cuDoubleComplex(0.1, 0.0),
            make_cuDoubleComplex(0.7, 0.0),
            make_cuDoubleComplex(1.0, 0.0),
            make_cuDoubleComplex(0.9, 0.0),
            make_cuDoubleComplex(0.4, 0.0),
            make_cuDoubleComplex(0.3, 0.0),
        };

        std::vector<double> m_next_gpu(m_curr.size(), 0.0);
        std::vector<double> m_next_cpu(m_curr.size(), 0.0);
        std::vector<int> sm_ids(material_size, -1);

        double* d_m_curr = nullptr;
        cuDoubleComplex* d_h_curr = nullptr;
        double* d_m_next = nullptr;
        int* d_sm_ids = nullptr;

        int device = -1;
        check_cuda(cudaGetDevice(&device), "cudaGetDevice");
        cudaDeviceProp props{};
        check_cuda(cudaGetDeviceProperties(&props, device), "cudaGetDeviceProperties");
        std::cout << "CUDA device " << device << ": " << props.name
                  << " (SM count=" << props.multiProcessorCount << ")\n";

        check_cuda(cudaMalloc(&d_m_curr, m_curr.size() * sizeof(double)), "cudaMalloc d_m_curr");
        check_cuda(cudaMalloc(&d_h_curr, h_curr.size() * sizeof(cuDoubleComplex)), "cudaMalloc d_h_curr");
        check_cuda(cudaMalloc(&d_m_next, m_next_gpu.size() * sizeof(double)), "cudaMalloc d_m_next");
        check_cuda(cudaMalloc(&d_sm_ids, sm_ids.size() * sizeof(int)), "cudaMalloc d_sm_ids");

        check_cuda(cudaMemcpy(d_m_curr, m_curr.data(), m_curr.size() * sizeof(double), cudaMemcpyHostToDevice), "copy M_curr");
        check_cuda(cudaMemcpy(d_h_curr, h_curr.data(), h_curr.size() * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice), "copy H_curr");

        debug_thread_location_kernel<<<blocks, threads_per_block>>>(d_sm_ids, material_size);
        check_cuda(cudaGetLastError(), "launch debug_thread_location_kernel");
        check_cuda(cudaDeviceSynchronize(), "sync debug_thread_location_kernel");
        check_cuda(cudaMemcpy(sm_ids.data(), d_sm_ids, sm_ids.size() * sizeof(int), cudaMemcpyDeviceToHost), "copy sm_ids");

        std::set<int> unique_sms(sm_ids.begin(), sm_ids.end());
        std::cout << "Host view: thread -> SM mapping\n";
        for (int idx = 0; idx < material_size; ++idx) {
            std::cout << "  thread " << idx << " -> sm " << sm_ids[idx] << '\n';
        }
        std::cout << "Unique SMs used: " << unique_sms.size() << '\n';

        LLG_RK4_kernel<<<blocks, threads_per_block>>>(
            d_m_curr, d_h_curr, d_m_next, material_size, dt, neg_gamma_ll, neg_coeff_damp);
        check_cuda(cudaGetLastError(), "launch LLG_RK4_kernel");
        check_cuda(cudaDeviceSynchronize(), "sync LLG_RK4_kernel");

        check_cuda(cudaMemcpy(m_next_gpu.data(), d_m_next, m_next_gpu.size() * sizeof(double), cudaMemcpyDeviceToHost), "copy M_next");

        llg_rk4_host(m_curr, h_curr, m_next_cpu, material_size, dt, neg_gamma_ll, neg_coeff_damp);

        double max_abs_err = 0.0;
        std::cout << std::setprecision(16);
        for (int idx = 0; idx < material_size; ++idx) {
            const int ptr = idx * 3;
            std::cout << "cell " << idx << '\n';
            std::cout << "  gpu: [" << m_next_gpu[ptr] << ", " << m_next_gpu[ptr + 1] << ", " << m_next_gpu[ptr + 2] << "]\n";
            std::cout << "  cpu: [" << m_next_cpu[ptr] << ", " << m_next_cpu[ptr + 1] << ", " << m_next_cpu[ptr + 2] << "]\n";

            for (int comp = 0; comp < 3; ++comp) {
                const double err = std::abs(m_next_gpu[ptr + comp] - m_next_cpu[ptr + comp]);
                if (err > max_abs_err) {
                    max_abs_err = err;
                }
            }
        }

        std::cout << "max_abs_err=" << max_abs_err << '\n';

        cudaFree(d_m_curr);
        cudaFree(d_h_curr);
        cudaFree(d_m_next);
        cudaFree(d_sm_ids);

        return max_abs_err < 1e-12 ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
