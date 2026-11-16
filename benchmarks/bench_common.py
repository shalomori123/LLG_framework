"""
Shared helpers for profile_modes.py and scaling_sweep.py: paths, environment
capture, CUDA-graph device timing, kernel occupancy, the analytic traffic/FLOP
model and host RAM/CPU sampling.
"""
import json
import math
import os
import platform
import re
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GPU_TOOLS = ROOT / "gpu_tools"
PYTHON_CORE = ROOT / "python_core"
RESULTS = Path(__file__).resolve().parent / "results"
for _p in (GPU_TOOLS, PYTHON_CORE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np

# Ideal DRAM traffic of one time step, every value read/written once (see maxwell_kernel.cu):
#   magnetic: read Ex Ey Hx Hy Hz, write Hx Hy Hz   -> 8 complex per cell
#   electric: read Ex Ey Hx Hy, ca cb, write Ex Ey  -> 6 complex + 2 real per cell
#   LLG:      read + write Mx My Mz                 -> 6 real per material cell
COMPLEX_PER_CELL = 14
REAL_PER_CELL = 2
REAL_PER_MATERIAL_CELL = 6
# Estimated FLOPs (hand count, verify with ncu sm__sass_thread_inst_executed_op_*):
#   Faraday 12 + Ampere 16 per cell; RK4 LL step + |M| + coupling per material cell
FLOPS_PER_CELL = 28
FLOPS_PER_MATERIAL_CELL = 190

THREADS_PER_BLOCK = 256


def results_dir(kind: str, gpu_name: str) -> Path:
    slug = re.sub(r"[^A-Za-z0-9]+", "-", gpu_name).strip("-") or "gpu"
    path = RESULTS / f"{datetime.now():%Y%m%d-%H%M%S}_{slug}_{kind}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _decode(value):
    return value.decode(errors="replace") if isinstance(value, bytes) else value


def environment() -> dict:
    import cupy as cp
    props = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
    get = lambda key, default=None: props.get(key, default)
    peak_bw = None
    try:  # memoryClockRate [kHz] x bus width [bit] x 2 (DDR) / 8
        peak_bw = 2 * get("memoryClockRate") * 1e3 * get("memoryBusWidth") / 8
    except TypeError:
        pass
    env = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "gpu": _decode(get("name")),
        "compute_capability": f"{get('major')}.{get('minor')}",
        "sm_count": get("multiProcessorCount"),
        "vram_total_gb": get("totalGlobalMem", 0) / 2**30,
        "l2_cache_mb": get("l2CacheSize", 0) / 2**20,
        "regs_per_sm": get("regsPerMultiprocessor", 65536),
        "max_threads_per_sm": get("maxThreadsPerMultiProcessor", 2048),
        "max_blocks_per_sm": get("maxBlocksPerMultiProcessor", 32),
        "sm_clock_mhz": (get("clockRate") or 0) / 1e3,
        "peak_dram_bandwidth_gbs": peak_bw / 1e9 if peak_bw else None,
        "cuda_driver": cp.cuda.runtime.driverGetVersion(),
        "cuda_runtime": cp.cuda.runtime.runtimeGetVersion(),
        "cupy": cp.__version__,
        "numpy": np.__version__,
        "python": platform.python_version(),
        "os": platform.platform(),
        "cpu": platform.processor() or platform.machine(),
        "cpu_count": os.cpu_count(),
    }
    try:
        import psutil
        env["ram_total_gb"] = psutil.virtual_memory().total / 2**30
    except ImportError:
        pass
    return env


def kernel_attributes(kernel) -> dict:
    """num_regs / local_size_bytes / shared_size_bytes of a RawKernel or a
    RawModule.get_function() Function (which has no .attributes)."""
    if hasattr(kernel, "attributes"):
        return kernel.attributes
    import cupy as cp
    # CUfunction_attribute: 1 shared, 3 local, 4 num regs
    get = lambda attr: cp.cuda.driver.funcGetAttribute(attr, kernel.ptr)
    return {"num_regs": get(4), "local_size_bytes": get(3), "shared_size_bytes": get(1)}


def kernel_info(kernel, env: dict) -> dict:
    """Registers, spills, shared memory and theoretical occupancy at 256 threads/block."""
    attrs = kernel_attributes(kernel)
    regs = attrs["num_regs"]
    warps_per_block = THREADS_PER_BLOCK // 32
    regs_per_warp = math.ceil(regs * 32 / 256) * 256          # allocation granularity 256
    blocks_by_regs = (env["regs_per_sm"] // regs_per_warp) // warps_per_block if regs else 10**9
    blocks_by_threads = env["max_threads_per_sm"] // THREADS_PER_BLOCK
    blocks_per_sm = max(0, min(blocks_by_regs, blocks_by_threads, env["max_blocks_per_sm"]))
    return {
        "num_regs": regs,
        "local_bytes_spill": attrs["local_size_bytes"],
        "shared_bytes": attrs["shared_size_bytes"],
        "blocks_per_sm": blocks_per_sm,
        "occupancy": blocks_per_sm * warps_per_block * 32 / env["max_threads_per_sm"],
    }


def traffic_model(grid_size: int, material_size: int, cases: int, precision: str) -> dict:
    real_b = 8 if precision == "fp64" else 4
    bytes_per_step = cases * (grid_size * (COMPLEX_PER_CELL * 2 * real_b + REAL_PER_CELL * real_b)
                              + material_size * REAL_PER_MATERIAL_CELL * real_b)
    flops_per_step = cases * (grid_size * FLOPS_PER_CELL + material_size * FLOPS_PER_MATERIAL_CELL)
    return {
        "bytes_per_step": bytes_per_step,
        "flops_per_step": flops_per_step,
        "arithmetic_intensity": flops_per_step / bytes_per_step,
    }


def device_time_per_step(state, pulses, which: str = "step", graph_steps: int = 500,
                         min_seconds: float = 0.5):
    """Pure GPU time per step [s]: capture graph_steps launches into a CUDA graph and
    time graph replays with events, so host launch cost is excluded."""
    import cupy as cp
    launch = {
        "step": lambda k: state.step(pulses[k]),
        "magnetic": lambda k: state.launch_magnetic(pulses[k]),
        "electric": lambda k: state.launch_electric(),
    }[which]
    stream = cp.cuda.Stream(non_blocking=True)
    try:
        with stream:
            stream.begin_capture()
            for k in range(graph_steps):
                launch(k)
            graph = stream.end_capture()
    except Exception as exc:  # old CuPy / driver without graph capture
        return None, f"graph capture failed: {exc}"

    def replay(n):
        start, end = cp.cuda.Event(), cp.cuda.Event()
        start.record(stream)
        for _ in range(n):
            graph.launch(stream)
        end.record(stream)
        end.synchronize()
        return cp.cuda.get_elapsed_time(start, end) / 1e3

    replay(1)  # warm-up
    once = max(replay(1), 1e-6)
    repeats = int(min(200, max(1, math.ceil(min_seconds / once))))
    return replay(repeats) / (repeats * graph_steps), None


class HostMonitor:
    """Peak RSS (sampled in a thread, needs psutil) and process CPU time during a block."""

    def __init__(self, interval: float = 0.02):
        self.interval = interval
        self.peak_rss = self.base_rss = None
        self._stop = threading.Event()
        try:
            import psutil
            self._proc = psutil.Process()
        except ImportError:
            self._proc = None

    def _sample(self):
        while not self._stop.wait(self.interval):
            self.peak_rss = max(self.peak_rss, self._proc.memory_info().rss)

    def __enter__(self):
        if self._proc is not None:
            self.base_rss = self.peak_rss = self._proc.memory_info().rss
            self._thread = threading.Thread(target=self._sample, daemon=True)
            self._thread.start()
        self._cpu0, self._wall0 = time.process_time(), time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.wall = time.perf_counter() - self._wall0
        self.cpu = time.process_time() - self._cpu0
        if self._proc is not None:
            self._stop.set()
            self._thread.join()
            self.peak_rss = max(self.peak_rss, self._proc.memory_info().rss)

    def summary(self) -> dict:
        out = {"host_cpu_pct": 100 * self.cpu / self.wall if self.wall else None}
        if self.base_rss is not None:
            out["ram_peak_mb"] = self.peak_rss / 2**20
            out["ram_growth_mb"] = (self.peak_rss - self.base_rss) / 2**20
        return out


def all_cases():
    """The 12 production cases: every material, closed and open."""
    import sim_config as cfg
    return [dict(cfg.material_case(m), closed_system=cfg.STATES[s])
            for m in cfg.ALL_MATERIALS for s in ("close", "open")]


def cases_for_batch(C: int):
    combos = all_cases()
    return [combos[i % len(combos)] for i in range(C)]


def run_production(preset: dict, cases, precision: str, steps: int, save_time=(), output_paths=None) -> float:
    """Wall time [s] of sim.simulation_cases exactly as sim_run.py calls it (incl. setup)."""
    import cupy as cp
    import simulation_gpu as sim
    cp.cuda.Device().synchronize()
    t0 = time.perf_counter()
    sim.simulation_cases(
        cases, preset["z_indexes"], preset["ms0"], preset["max_E"], preset["phase_diff"], preset["fwhm"],
        preset["save_locations"], list(save_time), preset["dt"], steps,
        wavelength=preset["wavelength"], safety_factor=preset["safety_start"],
        time_interval=preset["time_interval"], precision=precision, output_paths=output_paths,
    )
    cp.cuda.Device().synchronize()
    return time.perf_counter() - t0


def write_json(path: Path, data):
    def default(value):
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
        return str(value)
    path.write_text(json.dumps(data, indent=2, default=default))


def write_csv(path: Path, rows: list[dict]):
    import csv
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
