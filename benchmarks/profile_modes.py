"""
Performance profile of every run mode (preset x precision x cases per launch).

For each mode it measures, on the production code path of sim_run.py:
  time      production us/step (frames + snapshot I/O to disk, like a real run),
            bare launch loop us/step (kernels only), pure GPU us/step from a
            CUDA graph (no host cost), and per-kernel GPU time
  host      GPU-busy %, host overhead per step, host CPU %, a cProfile of the
            loop grouped into launch / frame capture / flush / pulse / tqdm
  memory    VRAM held by the CuPy pool and model, host RAM peak, bytes written
  hardware  registers, spills, occupancy, SM fill, effective DRAM bandwidth vs
            peak, estimated GFLOP/s and arithmetic intensity
and a one-line bottleneck verdict. Optional --nsys / --ncu wrap each mode in
Nsight Systems / Nsight Compute and save their reports next to the results.

Output: benchmarks/results/<date-time>_<gpu>_modes/
    environment.json  summary.csv  summary.json  report.md
    modes/<mode>/metrics.json  cprofile.txt  cprofile.prof
    nsys/<mode>.nsys-rep (+ _*.csv)   ncu/<mode>.ncu-rep (+ .csv)

Examples:
    python benchmarks/profile_modes.py
    python benchmarks/profile_modes.py --presets rt --precisions fp64 --batches 1 12
    python benchmarks/profile_modes.py --steps 2000 --nsys --ncu
"""
import argparse
import cProfile
import io
import os
import pstats
import shutil
import subprocess
import sys
import time
from pathlib import Path

import bench_common as bc

DEFAULT_STEPS = {"full_width": 20000, "rt": 5000}
LAUNCH_LOOP_STEPS = 2000
GRAPH_STEPS = 500


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--presets", nargs="+", default=["full_width", "rt"])
    parser.add_argument("--precisions", nargs="+", choices=["fp64", "fp32"], default=["fp64", "fp32"])
    parser.add_argument("--batches", nargs="+", type=int, default=[1, 12], help="cases per launch")
    parser.add_argument("--steps", type=int, default=None,
                        help=f"time steps of the production run (default {DEFAULT_STEPS})")
    parser.add_argument("--profile-steps", type=int, default=3000, help="steps of the cProfile run")
    parser.add_argument("--nsys", action="store_true", help="also record an Nsight Systems timeline per mode")
    parser.add_argument("--ncu", action="store_true", help="also record Nsight Compute kernel reports per mode")
    parser.add_argument("--keep-bundles", action="store_true", help="keep the .simdata written by the run")
    parser.add_argument("--quiet", action="store_true", help="hide tqdm bars (also removes their cost)")
    parser.add_argument("--child", metavar="MODE", help=argparse.SUPPRESS)
    return parser.parse_args()


def mode_name(preset, precision, C):
    return f"{preset}-{precision}-C{C}"


def parse_mode(name):
    preset, precision, c = name.rsplit("-", 2)
    return preset, precision, int(c[1:])


def scaled_snapshots(preset_name, steps):
    """Keep the preset's snapshots-per-step rate when the run is truncated."""
    import numpy as np
    import sim_config as cfg
    full = cfg.build_preset(preset_name)
    count = int(round(len(full["save_time"]) * steps / full["simulation_time"]))
    return np.linspace(0, steps - 1, count).astype(int) if count else np.array([], dtype=int)


def dir_bytes(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.exists() else 0


def host_breakdown(prof: cProfile.Profile) -> dict:
    """Seconds per host activity. Cython calls (kernel __call__, cp.take) are only
    visible if CuPy exposes them to the profiler; otherwise they stay in loop self time."""
    stats = pstats.Stats(prof).stats
    groups = {"launch": 0.0, "frame_capture": 0.0, "flush": 0.0, "pulse": 0.0, "snapshot_io": 0.0,
              "tqdm": 0.0, "loop_self": 0.0}
    for (file, _, func), (_, _, tt, ct, _) in stats.items():
        if func in ("launch_magnetic", "launch_electric"):
            groups["launch"] += ct
        elif func == "flush":
            groups["flush"] += ct
        elif func in ("out_pulse",):
            groups["pulse"] += ct
        elif func == "simulation_cases":
            groups["loop_self"] += tt
        elif "take" in func:
            groups["frame_capture"] += tt
        elif "memmap" in file or "format.py" in file:
            groups["snapshot_io"] += tt
        elif "tqdm" in file:
            groups["tqdm"] += tt
    groups["total"] = sum(tt for (_, _, tt, _, _) in stats.values())
    return groups


def verdict(m) -> str:
    if m["gpu_busy_pct"] is not None and m["gpu_busy_pct"] < 50:
        if m["io_extraction_us"] > 0.3 * m["production_us_per_step"]:
            return "host-bound: frame capture / output I/O"
        return "host-bound: kernel launch (Python + driver)"
    if m["bandwidth_pct_peak"] and m["bandwidth_pct_peak"] > 50:
        return "GPU DRAM-bandwidth-bound"
    if m["sm_fill"] < 1:
        return "GPU under-filled (fewer blocks than SMs): latency-bound"
    return "GPU latency/compute mix: inspect ncu report"


def profile_mode(preset_name, precision, C, steps, profile_steps, mode_dir: Path, env, keep_bundles):
    import numpy as np
    import cupy as cp
    from cupyx.profiler import time_range
    import sim_config as cfg
    import simulation_gpu as sim

    preset = cfg.build_preset(preset_name, steps)
    save_time = scaled_snapshots(preset_name, steps)
    cases = bc.cases_for_batch(C)
    N = preset["grid_size"]
    mat = preset["z_indexes"][-1] - preset["z_indexes"][0]
    m = {"mode": mode_name(preset_name, precision, C), "preset": preset_name, "precision": precision,
         "cases": C, "grid_cells": N, "material_cells": mat, "steps": steps, "snapshots": len(save_time)}

    # Kernels: compile/load time (CuPy's disk cache makes repeat runs fast) and attributes
    t0 = time.perf_counter()
    magnetic, electric = sim.load_kernels(precision)
    m["kernel_load_s"] = time.perf_counter() - t0
    for name, kernel in (("magnetic", magnetic), ("electric", electric)):
        for key, value in bc.kernel_info(kernel, env).items():
            m[f"{key}_{name}"] = value
    m["blocks_total"] = -(-N // bc.THREADS_PER_BLOCK) * C
    m["sm_fill"] = m["blocks_total"] / env["sm_count"]
    m["waves"] = m["blocks_total"] / (env["sm_count"] * max(1, m["blocks_per_sm_magnetic"]))

    # Warm-up (allocations, first launches)
    bc.run_production(preset, cases, precision, min(200, steps))

    # 1. Production run: frames + scaled snapshots streamed to .simdata like sim_run.py
    pool, pinned = cp.get_default_memory_pool(), cp.get_default_pinned_memory_pool()
    pool.free_all_blocks(); pinned.free_all_blocks()
    scratch = mode_dir / "bundles"
    outputs = [scratch / f"case{c:02d}" for c in range(C)]
    with time_range(m["mode"]), bc.HostMonitor() as mon:
        elapsed = bc.run_production(preset, cases, precision, steps, save_time, outputs)
    free_b, total_b = cp.cuda.Device().mem_info
    m.update(mon.summary())
    m["vram_pool_mb"] = pool.total_bytes() / 2**20
    m["pinned_host_mb"] = pinned.total_bytes() / 2**20 if hasattr(pinned, "total_bytes") else None
    m["vram_device_used_mb"] = (total_b - free_b) / 2**20
    m["vram_model_mb"] = (C * sim.case_device_bytes(N, mat, precision) + sim.FRAME_BUFFER_BYTES) / 2**20
    m["io_mb_written"] = dir_bytes(scratch) / 2**20
    if not keep_bundles:
        shutil.rmtree(scratch, ignore_errors=True)

    m["production_us_per_step"] = 1e6 * elapsed / steps
    m["us_per_case_step"] = m["production_us_per_step"] / C
    m["ns_per_cell_step"] = 1e9 * elapsed / (steps * N * C)
    m["cell_updates_per_s"] = steps * N * C / elapsed

    # 2. Bare launch loop and pure device time on a fresh state
    state = sim.DeviceState(cases, preset["z_indexes"], preset["ms0"], preset["dt"],
                            preset["wavelength"], precision)
    pulse_func = sim.pulse_gen(preset["max_E"], preset["phase_diff"], preset["fwhm"], preset["dt"],
                               preset["wavelength"], 0, preset["safety_start"])
    k_steps = max(LAUNCH_LOOP_STEPS, GRAPH_STEPS)
    pulses = cp.asarray(pulse_func(np.arange(k_steps)).astype(state.complex_t))
    cp.cuda.Device().synchronize()
    t0 = time.perf_counter()
    for k in range(LAUNCH_LOOP_STEPS):
        state.step(pulses[k])
    cp.cuda.Device().synchronize()
    m["launch_loop_us_per_step"] = 1e6 * (time.perf_counter() - t0) / LAUNCH_LOOP_STEPS

    errors = []
    for which, key in (("step", "device_us_per_step"), ("magnetic", "magnetic_us"), ("electric", "electric_us")):
        seconds, err = bc.device_time_per_step(state, pulses, which, GRAPH_STEPS)
        m[key] = 1e6 * seconds if seconds else None
        if err:
            errors.append(err)
    if errors:
        m["notes"] = "; ".join(sorted(set(errors)))

    dev = m["device_us_per_step"]
    m["gpu_busy_pct"] = 100 * dev / m["production_us_per_step"] if dev else None
    m["host_overhead_us"] = m["production_us_per_step"] - dev if dev else None
    m["launch_overhead_us"] = max(0.0, m["launch_loop_us_per_step"] - dev) if dev else None
    m["io_extraction_us"] = m["production_us_per_step"] - m["launch_loop_us_per_step"]

    m.update(bc.traffic_model(N, mat, C, precision))
    if dev:
        m["bandwidth_gbs"] = m["bytes_per_step"] / (dev * 1e-6) / 1e9
        peak = env.get("peak_dram_bandwidth_gbs")
        m["bandwidth_pct_peak"] = 100 * m["bandwidth_gbs"] / peak if peak else None
        m["gflops"] = m["flops_per_step"] / (dev * 1e-6) / 1e9
    else:
        m["bandwidth_gbs"] = m["bandwidth_pct_peak"] = m["gflops"] = None
    del state, pulses

    # 3. Host-side function profile
    prof = cProfile.Profile()
    profile_steps = min(profile_steps, steps)
    prof.enable()
    bc.run_production(preset, cases, precision, profile_steps,
                      scaled_snapshots(preset_name, profile_steps), outputs)
    prof.disable()
    shutil.rmtree(scratch, ignore_errors=True)
    prof.dump_stats(mode_dir / "cprofile.prof")
    text = io.StringIO()
    pstats.Stats(prof, stream=text).sort_stats("tottime").print_stats(30)
    (mode_dir / "cprofile.txt").write_text(text.getvalue())
    breakdown = host_breakdown(prof)
    for key, seconds in breakdown.items():
        m[f"host_{key}_us_per_step"] = 1e6 * seconds / profile_steps

    m["bottleneck"] = verdict(m)
    return m


def run_child(mode, steps):
    """One production run, the target of nsys / ncu."""
    import sim_config as cfg
    from cupyx.profiler import time_range
    preset_name, precision, C = parse_mode(mode)
    preset = cfg.build_preset(preset_name, steps)
    cases = bc.cases_for_batch(C)
    bc.run_production(preset, cases, precision, min(200, steps))
    with time_range(mode):
        bc.run_production(preset, cases, precision, steps, scaled_snapshots(preset_name, steps))


def run_tool(cmd, log: Path):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        log.write_text(" ".join(map(str, cmd)) + "\n\n" + result.stdout + "\n" + result.stderr)
        return result.returncode == 0
    except OSError as exc:
        log.write_text(f"{cmd}\n{exc}")
        return False


def nsight(mode, steps, out: Path, use_nsys, use_ncu):
    child = [sys.executable, str(Path(__file__).resolve()), "--child", mode, "--steps", str(steps), "--quiet"]
    if use_nsys and shutil.which("nsys"):
        rep = out / "nsys" / mode
        rep.parent.mkdir(exist_ok=True)
        if run_tool(["nsys", "profile", "--trace=cuda,nvtx,osrt", "--sample=cpu", "--force-overwrite=true",
                     "-o", str(rep), *child], rep.with_suffix(".log")):
            for reports in ("cuda_gpu_kern_sum,cuda_api_sum,nvtx_sum", "gpukernsum,cudaapisum,nvtxsum"):
                if run_tool(["nsys", "stats", "--force-export=true", "--report", reports, "--format", "csv",
                             "--output", str(rep), str(rep) + ".nsys-rep"], rep.with_suffix(".stats.log")):
                    break
    elif use_nsys:
        print("  nsys not on PATH, skipped")
    if use_ncu and shutil.which("ncu"):
        rep = out / "ncu" / mode
        rep.parent.mkdir(exist_ok=True)
        ncu_child = child[:]
        ncu_child[ncu_child.index("--steps") + 1] = "300"
        if run_tool(["ncu", "--set", "full", "--kernel-name", "regex:magnetic_kernel|electric_kernel",
                     "--launch-skip", "400", "--launch-count", "2", "-f", "-o", str(rep), *ncu_child],
                    rep.with_suffix(".log")):
            run_tool(["ncu", "--import", str(rep) + ".ncu-rep", "--csv", "--page", "details"],
                     rep.with_suffix(".csv"))
    elif use_ncu:
        print("  ncu not on PATH, skipped")


REPORT_COLUMNS = [
    ("mode", "mode", "{}"), ("grid_cells", "cells", "{:,}"),
    ("production_us_per_step", "prod us/step", "{:.2f}"), ("device_us_per_step", "GPU us/step", "{:.2f}"),
    ("launch_loop_us_per_step", "launch loop us/step", "{:.2f}"), ("gpu_busy_pct", "GPU busy %", "{:.0f}"),
    ("ns_per_cell_step", "ns/cell-step", "{:.3f}"), ("bandwidth_pct_peak", "DRAM % peak", "{:.0f}"),
    ("sm_fill", "SM fill", "{:.2f}"), ("num_regs_magnetic", "regs mag", "{}"),
    ("occupancy_magnetic", "occ mag", "{:.2f}"), ("vram_pool_mb", "VRAM MB", "{:.0f}"),
    ("ram_peak_mb", "RAM MB", "{:.0f}"), ("bottleneck", "bottleneck", "{}"),
]


def fmt(value, spec):
    if value is None or value == "":
        return "-"
    try:
        return spec.format(value)
    except (ValueError, TypeError):
        return str(value)


def write_report(out: Path, env, rows):
    lines = [f"# Mode profile - {env['gpu']}", "",
             f"{env['timestamp']}  |  CC {env['compute_capability']}, {env['sm_count']} SMs, "
             f"{env['vram_total_gb']:.0f} GB, peak DRAM {fmt(env.get('peak_dram_bandwidth_gbs'), '{:.0f}')} GB/s  |  "
             f"CUDA runtime {env['cuda_runtime']}, CuPy {env['cupy']}, {env['cpu']} x{env['cpu_count']}", "",
             "| " + " | ".join(h for _, h, _ in REPORT_COLUMNS) + " |",
             "|" + "---|" * len(REPORT_COLUMNS)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(k), s) for k, _, s in REPORT_COLUMNS) + " |")
    lines += ["", "## Host time per step (cProfile, us)", "",
              "| mode | launch | frame capture | flush | pulse | snapshot I/O | tqdm | loop self |",
              "|---|---|---|---|---|---|---|---|"]
    for row in rows:
        keys = ["launch", "frame_capture", "flush", "pulse", "snapshot_io", "tqdm", "loop_self"]
        lines.append(f"| {row['mode']} | " + " | ".join(fmt(row.get(f"host_{k}_us_per_step"), "{:.2f}")
                                                     for k in keys) + " |")
    lines += ["", "## Columns", "",
              "- **prod us/step**: wall time of the production loop (frames + snapshots to disk) per step.",
              "- **GPU us/step**: pure device time from a CUDA-graph replay, no host launch cost.",
              "- **launch loop us/step**: Python loop of the two kernel launches only. If it is well above "
              "GPU us/step, the host cannot issue launches fast enough (launch-bound).",
              "- **GPU busy %** = GPU / prod. Below 50 % the run is host-bound.",
              "- **DRAM % peak**: analytic bytes per step (bench_common.traffic_model) / GPU time, vs "
              "memoryClockRate x bus width. Small grids live in L2, so this is only meaningful for rt.",
              "- **SM fill** = thread blocks / SMs. Below 1 most of the GPU is idle whatever the code does.",
              "- Host breakdown: Cython calls (kernel __call__, cp.take) may not be visible to cProfile and "
              "then show up as loop self time. The nsys timeline is the ground truth.", ""]
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    if args.quiet:
        os.environ["TQDM_DISABLE"] = "1"
    if args.child:
        run_child(args.child, args.steps or 2000)
        return

    env = bc.environment()
    out = bc.results_dir("modes", env["gpu"])
    bc.write_json(out / "environment.json", env)
    print(f"{env['gpu']}  ->  {out}")

    rows = []
    for preset in args.presets:
        steps = args.steps or DEFAULT_STEPS.get(preset, 5000)
        for precision in args.precisions:
            for C in args.batches:
                name = mode_name(preset, precision, C)
                print(f"\n=== {name}  ({steps} steps)")
                mode_dir = out / "modes" / name
                mode_dir.mkdir(parents=True, exist_ok=True)
                try:
                    m = profile_mode(preset, precision, C, steps, args.profile_steps, mode_dir, env,
                                     args.keep_bundles)
                except Exception as exc:  # e.g. out of memory: record and continue
                    m = {"mode": name, "preset": preset, "precision": precision, "cases": C,
                         "bottleneck": f"failed: {type(exc).__name__}: {exc}"}
                bc.write_json(mode_dir / "metrics.json", m)
                rows.append(m)
                print(f"  prod {fmt(m.get('production_us_per_step'), '{:.2f}')} us/step | "
                      f"GPU {fmt(m.get('device_us_per_step'), '{:.2f}')} us/step | "
                      f"busy {fmt(m.get('gpu_busy_pct'), '{:.0f}')}% | {m['bottleneck']}")
                if args.nsys or args.ncu:
                    nsight(name, min(steps, 5000), out, args.nsys, args.ncu)
                bc.write_csv(out / "summary.csv", rows)
                bc.write_json(out / "summary.json", {"environment": env, "modes": rows})
                write_report(out, env, rows)
    print(f"\nreport: {out / 'report.md'}")


if __name__ == "__main__":
    main()
