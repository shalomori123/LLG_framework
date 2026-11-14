"""
Space scaling: time per step vs number of grid cells N, CPU reference vs GPU.

The grid is the full_width preset (dz = 0.25 nm, same 47-cell material stack)
with only the vacuum gap z1 widened, so N runs from 78 cells (= full_width)
through 240,889 (= the rt grid size) up to ~1e6. The real rt preset (dz = 1 nm)
is measured as an extra point to check that the cell count, not dz, sets the cost.

Series per N:
  cpu                python_core/simulation_cpu.py (NumPy reference, 1 case)
  gpu_<p>            production loop of simulation_gpu (1 case), p = fp64/fp32
  gpu_<p>_device     pure GPU time per step (CUDA-graph replay, no host cost)
  gpu_fp64_C<C>      production loop with C cases per launch, time per case

Output: benchmarks/results/<date-time>_<gpu>_scaling/
    scaling.csv  scaling.json  environment.json
    fig_scaling_time.png        time per step and case [us] vs N
    fig_scaling_throughput.png  cell updates per second vs N
    fig_scaling_speedup.png     speedup over the CPU reference vs N

Examples:
    python benchmarks/scaling_sweep.py
    python benchmarks/scaling_sweep.py --sizes 78 1000 240889 --target-seconds 1
    python benchmarks/scaling_sweep.py --plot-only benchmarks/results/<dir> --figures-dir book/figures
"""
import argparse
import csv
import os
import shutil
import time
from pathlib import Path

import numpy as np

import bench_common as bc

FULL_WIDTH_CELLS = 78
RT_CELLS = 240_889


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sizes", nargs="+", type=int, default=None,
                        help="grid sizes N (default: 13 log-spaced from 78 to 1e6, plus 240,889)")
    parser.add_argument("--max-cells", type=int, default=1_000_000)
    parser.add_argument("--precisions", nargs="+", choices=["fp64", "fp32"], default=["fp64", "fp32"])
    parser.add_argument("--batch", type=int, default=12, help="cases per launch for the batched series, 0 = skip")
    parser.add_argument("--no-cpu", action="store_true", help="skip the CPU reference")
    parser.add_argument("--cpu-max-cells", type=int, default=None, help="skip CPU above this N")
    parser.add_argument("--target-seconds", type=float, default=2.0, help="GPU run length per point")
    parser.add_argument("--cpu-target-seconds", type=float, default=3.0, help="CPU run length per point")
    parser.add_argument("--no-rt-preset", action="store_true", help="skip the real rt preset point")
    parser.add_argument("--plot-only", type=Path, default=None, metavar="DIR",
                        help="re-plot an existing results dir (reads scaling.csv)")
    parser.add_argument("--figures-dir", type=Path, default=None, help="also copy the figures here")
    parser.add_argument("--quiet", action="store_true", help="hide tqdm bars")
    return parser.parse_args()


def widened_grid(N: int) -> dict:
    """full_width preset with the vacuum gap set so the grid has exactly N cells."""
    import sim_config as cfg
    grid = cfg.build_preset("full_width")
    z = grid["z_indexes"]
    mag, cond = z[1] - z[0], z[2] - z[1]
    z1 = N - cfg.grid_size([0, mag, mag + cond], grid["dt"], grid["wavelength"])
    if z1 < 8:
        raise ValueError(f"N={N} is smaller than the full_width grid ({FULL_WIDTH_CELLS} cells)")
    grid["z_indexes"] = [z1, z1 + mag, z1 + mag + cond]
    grid["save_locations"] = np.array([4] + grid["z_indexes"])
    grid["grid_size"] = cfg.grid_size(grid["z_indexes"], grid["dt"], grid["wavelength"])
    assert grid["grid_size"] == N
    return grid


def gpu_point(grid, precision, C, target_s):
    """Production us per step (all C cases) and the number of steps timed."""
    cases = bc.cases_for_batch(C)
    calib = 200
    t = bc.run_production(grid, cases, precision, calib)  # also warms up
    steps = int(np.clip(target_s * calib / t, 1000, 200_000))
    t = bc.run_production(grid, cases, precision, steps)
    return 1e6 * t / steps, steps


def gpu_device_point(grid, precision):
    import cupy as cp
    import simulation_gpu as sim
    state = sim.DeviceState(bc.cases_for_batch(1), grid["z_indexes"], grid["ms0"], grid["dt"],
                            grid["wavelength"], precision)
    pulse_func = sim.pulse_gen(grid["max_E"], grid["phase_diff"], grid["fwhm"], grid["dt"],
                               grid["wavelength"], 0, grid["safety_start"])
    pulses = cp.asarray(pulse_func(np.arange(500)).astype(state.complex_t))
    seconds, _ = bc.device_time_per_step(state, pulses, "step", 500)
    return 1e6 * seconds if seconds else None


def cpu_point(grid, target_s):
    import simulation_cpu as cpu
    case = bc.cases_for_batch(1)[0]

    def run(steps):
        t0 = time.perf_counter()
        cpu.simulation(grid["z_indexes"], case["eps_r"], case["conductivity"], case["damping"], grid["ms0"],
                       grid["max_E"], grid["phase_diff"], grid["fwhm"], grid["save_locations"], [],
                       case["closed_system"], grid["dt"], steps, wavelength=grid["wavelength"],
                       safety_factor=grid["safety_start"], time_interval=grid["time_interval"])
        return time.perf_counter() - t0

    calib = 5
    t = run(calib)
    steps = int(np.clip(target_s * calib / t, calib, 20_000))
    return 1e6 * run(steps) / steps, steps


def measure(kind, grid, args):
    N = grid["grid_size"]
    row = {"kind": kind, "grid_cells": N, "material_cells": grid["z_indexes"][-1] - grid["z_indexes"][0],
           "dz_nm": 2 * 299_792_458 * grid["dt"] * 1e9}
    for p in args.precisions:
        row[f"gpu_{p}_us"], row[f"gpu_{p}_steps"] = gpu_point(grid, p, 1, args.target_seconds)
        row[f"gpu_{p}_device_us"] = gpu_device_point(grid, p)
    if args.batch > 1:
        p = args.precisions[0]
        us, row[f"gpu_{p}_C{args.batch}_steps"] = gpu_point(grid, p, args.batch, args.target_seconds)
        row[f"gpu_{p}_C{args.batch}_us"] = us / args.batch
    if not args.no_cpu and (args.cpu_max_cells is None or N <= args.cpu_max_cells):
        row["cpu_us"], row["cpu_steps"] = cpu_point(grid, args.cpu_target_seconds)
        for key in [k for k in row if k.startswith("gpu_") and k.endswith("_us")]:
            if row[key]:
                row["speedup_" + key[4:-3]] = row["cpu_us"] / row[key]
    return row


# ----------------------------------------------------------------------------- plots
SERIES = [  # column, label, style
    ("cpu_us", "CPU (NumPy reference)", dict(color="0.35", marker="s", linestyle="-")),
    ("gpu_fp64_us", "GPU fp64", dict(color="C0", marker="o", linestyle="-")),
    ("gpu_fp64_device_us", "GPU fp64, device time only", dict(color="C0", marker="", linestyle="--")),
    ("gpu_fp32_us", "GPU fp32", dict(color="C1", marker="^", linestyle="-")),
    ("gpu_fp32_device_us", "GPU fp32, device time only", dict(color="C1", marker="", linestyle="--")),
]


def batch_series(columns):
    out = []
    for col in columns:
        if col.startswith("gpu_") and "_C" in col and col.endswith("_us"):
            p, c = col[4:-3].split("_C")
            out.append((col, f"GPU {p}, {c} cases per launch (per case)",
                        dict(color="C2", marker="D", linestyle="-")))
    return out


def read_rows(csv_path: Path):
    with csv_path.open() as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key, value in row.items():
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                row[key] = value or None
    return [r for r in rows if not r.get("error")]


def mark_presets(ax):
    for n, name in ((FULL_WIDTH_CELLS, "full_width"), (RT_CELLS, "rt")):
        ax.axvline(n, color="0.5", linestyle=":", linewidth=1)
        ax.text(n * 1.08, 0.97, f"{name}\n{n:,} cells", transform=ax.get_xaxis_transform(),
                va="top", ha="left", fontsize=8, color="0.3")


def plot(out_dir: Path, figures_dir: Path = None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = read_rows(out_dir / "scaling.csv")
    sweep = sorted((r for r in rows if r["kind"] == "sweep"), key=lambda r: r["grid_cells"])
    presets = [r for r in rows if r["kind"] != "sweep"]
    columns = list(rows[0].keys()) if rows else []
    series = [s for s in SERIES + batch_series(columns) if s[0] in columns]

    def draw(transform, ylabel, name, speedup=False):
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        for col, label, style in series:
            if speedup and col == "cpu_us":
                continue
            ycol = "speedup_" + col[4:-3] if speedup else col
            pts = [(r["grid_cells"], transform(r, r.get(ycol))) for r in sweep if r.get(ycol)]
            if pts:
                ax.plot(*zip(*pts), label=label, markersize=4, linewidth=1.3, **style)
            for r in presets:  # real rt preset, same series
                if r.get(ycol):
                    ax.plot(r["grid_cells"], transform(r, r[ycol]), linestyle="none", marker="*",
                            markersize=11, color=style["color"], markeredgecolor="black", markeredgewidth=0.5)
        if presets:
            ax.plot([], [], linestyle="none", marker="*", markersize=11, color="white",
                    markeredgecolor="black", label="real rt preset (dz = 1 nm)")
        if speedup:
            ax.axhline(1, color="black", linewidth=0.8)
            for r in sweep:
                if r["grid_cells"] in (FULL_WIDTH_CELLS, RT_CELLS) and r.get("speedup_fp64"):
                    ax.annotate(f"x{r['speedup_fp64']:.0f}", (r["grid_cells"], r["speedup_fp64"]),
                                textcoords="offset points", xytext=(6, -12), fontsize=9, color="C0")
        mark_presets(ax)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("grid cells N")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout()
        path = out_dir / name
        fig.savefig(path, dpi=200)
        plt.close(fig)
        if figures_dir:
            figures_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, figures_dir / name)
        print(f"saved {path}")

    draw(lambda r, us: us, "time per step and case [µs]", "fig_scaling_time.png")
    draw(lambda r, us: r["grid_cells"] / (us * 1e-6), "cell updates per second [cells·steps/s]",
         "fig_scaling_throughput.png")
    draw(lambda r, s: s, "speedup over CPU reference", "fig_scaling_speedup.png", speedup=True)


def main():
    args = parse_args()
    if args.plot_only:
        plot(args.plot_only, args.figures_dir)
        return
    if args.quiet:
        os.environ["TQDM_DISABLE"] = "1"

    env = bc.environment()
    out = bc.results_dir("scaling", env["gpu"])
    bc.write_json(out / "environment.json", env)
    print(f"{env['gpu']}  ->  {out}")

    sizes = args.sizes or sorted({int(round(n)) for n in np.geomspace(FULL_WIDTH_CELLS, 1e6, 13)} | {RT_CELLS})
    sizes = [n for n in sizes if FULL_WIDTH_CELLS <= n <= args.max_cells]

    points = [("sweep", n) for n in sizes]
    if not args.no_rt_preset:
        points.append(("preset_rt", None))

    import sim_config as cfg
    rows = []
    for kind, n in points:
        grid = widened_grid(n) if n else cfg.build_preset("rt")
        print(f"\n=== {kind}  N = {grid['grid_size']:,}")
        try:
            row = measure(kind, grid, args)
        except Exception as exc:  # e.g. out of memory at the largest N
            row = {"kind": kind, "grid_cells": grid["grid_size"], "error": f"{type(exc).__name__}: {exc}"}
        rows.append(row)
        print("  " + "  ".join(f"{k}={v:.3g}" for k, v in row.items() if k.endswith("_us") and v))
        bc.write_csv(out / "scaling.csv", rows)
        bc.write_json(out / "scaling.json", {"environment": env, "points": rows})
    plot(out, args.figures_dir)


if __name__ == "__main__":
    main()
