#!/usr/bin/env python3
"""
Temporary Python host for llg_kernel:
- Launches a debug kernel to print SM ID per thread (proves GPU execution).
- Runs the LLG RK4 kernel.
- Compares GPU result against CPU reference.

Expected layout:
- M_curr: float64, shape (N, 3)
- H_curr: complex128, shape (3, N)  # component-major (Hx, Hy, Hz)
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import sys
from typing import Tuple

import numpy as np

try:
    import cupy as cp
except Exception as exc:  # pragma: no cover - runtime guard
    raise SystemExit(
        "CuPy is required for gpu_tools/llg_kernel_host_tmp.py.\n"
        f"Import error: {exc}\n"
        "Install a CUDA-matching CuPy build, e.g. pip install cupy-cuda12x"
    )


DEBUG_KERNEL_SRC = r"""
__device__ __forceinline__ unsigned int get_smid() {
    unsigned int smid;
    asm volatile("mov.u32 %0, %smid;" : "=r"(smid));
    return smid;
}

extern "C" __global__ void debug_thread_location_kernel(int* sm_ids, int material_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= material_size) return;
    unsigned int smid = get_smid();
    sm_ids[idx] = (int)smid;
    printf("GPU thread idx=%d block=%d thread=%d sm=%u\n", idx, blockIdx.x, threadIdx.x, smid);
}
"""


def build_cuda_source() -> str:
    kernel_path = Path(__file__).with_name("llg_kernel.cu")
    if not kernel_path.exists():
        raise SystemExit(f"Missing CUDA source: {kernel_path}")
    return kernel_path.read_text(encoding="utf-8") + "\n\n" + DEBUG_KERNEL_SRC


def detect_cuda_include_dirs() -> list[str]:
    candidates: list[Path] = []

    for env_key in ("CUDA_PATH", "CUDA_HOME"):
        value = os.environ.get(env_key)
        if value:
            candidates.append(Path(value) / "include")

    nvcc_path = shutil.which("nvcc")
    if nvcc_path:
        nvcc = Path(nvcc_path).resolve()
        # Typical Linux toolkit layout:
        #   /usr/local/cuda/bin/nvcc -> /usr/local/cuda/include
        #   /usr/bin/nvcc (alternatives) -> /usr/include or /usr/local/cuda/include
        candidates.append(nvcc.parent.parent / "include")
        candidates.append(nvcc.parent.parent / "targets" / "x86_64-linux" / "include")

    candidates.extend(
        [
            Path("/usr/local/cuda/include"),
            Path("/usr/include"),
            Path("/usr/include/x86_64-linux-gnu"),
        ]
    )

    include_dirs: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        if (candidate / "cuComplex.h").exists():
            candidate_str = str(candidate)
            if candidate_str not in seen:
                include_dirs.append(candidate_str)
                seen.add(candidate_str)
    return include_dirs


def cross_product_host(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.empty_like(a)
    out[:, 0] = a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1]
    out[:, 1] = a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2]
    out[:, 2] = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    return out


def llg_rk4_host(
    m_curr: np.ndarray,
    h_curr: np.ndarray,
    dt: float,
    neg_gamma_ll: float,
    neg_coeff_damp: float,
) -> np.ndarray:
    h_real = np.real(h_curr).T  # (N, 3)
    mx = m_curr

    def deriv(m_guess: np.ndarray) -> np.ndarray:
        c1 = cross_product_host(m_guess, h_real)
        c2 = cross_product_host(m_guess, c1)
        return neg_gamma_ll * c1 + neg_coeff_damp * c2

    k1 = deriv(mx)
    k2 = deriv(mx + 0.5 * dt * k1)
    k3 = deriv(mx + 0.5 * dt * k2)
    k4 = deriv(mx + dt * k3)
    return mx + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def decode_device_name(name_value) -> str:
    if isinstance(name_value, bytes):
        return name_value.split(b"\x00", 1)[0].decode(errors="replace")
    if hasattr(name_value, "tobytes"):
        return name_value.tobytes().split(b"\x00", 1)[0].decode(errors="replace")
    return str(name_value)


def make_inputs(n: int, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    m = rng.normal(size=(n, 3)).astype(np.float64)
    h = (rng.normal(size=(3, n)) + 1j * rng.normal(size=(3, n))).astype(np.complex128)
    return m, h


def main() -> int:
    parser = argparse.ArgumentParser(description="Temporary Python host for LLG CUDA kernel")
    parser.add_argument("--n", type=int, default=256, help="material_size (number of cells)")
    parser.add_argument("--threads", type=int, default=128, help="threads per block")
    parser.add_argument("--seed", type=int, default=7, help="random seed for synthetic inputs")
    parser.add_argument("--tol", type=float, default=1e-12, help="max abs error tolerance")
    args = parser.parse_args()

    if args.n <= 0:
        raise SystemExit("--n must be > 0")
    if args.threads <= 0:
        raise SystemExit("--threads must be > 0")

    dt = 1e-15
    neg_gamma_ll = -2.5e10
    neg_coeff_damp = -4.0e8

    try:
        props = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
    except cp.cuda.runtime.CUDARuntimeError as exc:
        if "cudaErrorNoDevice" in str(exc):
            print("No CUDA device visible to CuPy (cudaErrorNoDevice).", file=sys.stderr)
            print(f"LD_LIBRARY_PATH={os.environ.get('LD_LIBRARY_PATH', '')}", file=sys.stderr)
            print(f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES', '')}", file=sys.stderr)
            print("If you are on WSL, try:", file=sys.stderr)
            print("  export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}", file=sys.stderr)
            print("  unset CUDA_VISIBLE_DEVICES", file=sys.stderr)
            print(
                "  python3 -c \"import cupy as cp; print(cp.cuda.runtime.getDeviceCount())\"",
                file=sys.stderr,
            )
            return 2
        raise

    dev_name = decode_device_name(props["name"])
    sm_count = int(props["multiProcessorCount"])
    print(f"CUDA device: {dev_name} (SM count={sm_count})")

    m_curr, h_curr = make_inputs(args.n, args.seed)
    m_next_cpu = llg_rk4_host(m_curr, h_curr, dt, neg_gamma_ll, neg_coeff_damp)

    include_dirs = detect_cuda_include_dirs()
    if not include_dirs:
        raise SystemExit(
            "Could not find CUDA headers (cuComplex.h). "
            "Set CUDA_PATH/CUDA_HOME or install CUDA toolkit headers."
        )
    compile_options = ["-std=c++14"] + [f"-I{path}" for path in include_dirs]
    print("Using CUDA include dirs:")
    for path in include_dirs:
        print(f"  {path}")

    module = cp.RawModule(
        code=build_cuda_source(),
        options=tuple(compile_options),
        name_expressions=("LLG_RK4_kernel", "debug_thread_location_kernel"),
    )
    llg_kernel = module.get_function("LLG_RK4_kernel")
    debug_kernel = module.get_function("debug_thread_location_kernel")

    m_curr_d = cp.asarray(m_curr.ravel(), dtype=cp.float64)
    h_curr_d = cp.asarray(h_curr.ravel(), dtype=cp.complex128)
    m_next_d = cp.zeros_like(m_curr_d)
    sm_ids_d = cp.full((args.n,), -1, dtype=cp.int32)

    blocks = (args.n + args.threads - 1) // args.threads
    debug_kernel((blocks,), (args.threads,), (sm_ids_d, np.int32(args.n)))
    llg_kernel(
        (blocks,),
        (args.threads,),
        (
            m_curr_d,
            h_curr_d,
            m_next_d,
            np.int32(args.n),
            np.float64(dt),
            np.float64(neg_gamma_ll),
            np.float64(neg_coeff_damp),
        ),
    )
    cp.cuda.runtime.deviceSynchronize()

    sm_ids = cp.asnumpy(sm_ids_d)
    unique_sm = np.unique(sm_ids[sm_ids >= 0])
    print("Host view: first 32 thread -> SM IDs")
    preview = min(32, args.n)
    for idx in range(preview):
        print(f"  thread {idx:4d} -> sm {sm_ids[idx]}")
    print(f"Unique SMs used: {len(unique_sm)}")

    m_next_gpu = cp.asnumpy(m_next_d).reshape(args.n, 3)
    max_abs_err = float(np.max(np.abs(m_next_gpu - m_next_cpu)))
    print(f"max_abs_err={max_abs_err:.6e}")

    if max_abs_err >= args.tol:
        print(f"FAILED: error >= tolerance ({args.tol})", file=sys.stderr)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
