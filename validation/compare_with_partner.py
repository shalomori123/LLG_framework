"""
Compare this repo's GPU simulation against the partner repo's validated GPU backends.

Runs the Co_Pt case (open and close) on a truncated number of time steps with
identical inputs through both implementations and prints the max relative
difference of every output array.

Usage (from the LLG_framework root, on a CUDA machine):
    python validation/compare_with_partner.py --partner-repo ../llg_framework_cpu_version_updated
    python validation/compare_with_partner.py --preset rt --steps 20000
    python validation/compare_with_partner.py --precision fp32

Partner backends used:
    fp32 full_width -> simulation_gpu.py        fp32 rt -> simulation_gpu_rt.py
    fp64 full_width -> simulation_gpu_fp64.py   fp64 rt -> simulation_gpu_rt_fp64.py

Known, intentional difference: the partner backends double-buffer E, so the
saved value of the right ABC cell (E[N-1]) lags one step behind ours. That cell
is excluded from the "E times picture" comparison.
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
GPU_TOOLS = REPO_ROOT / "gpu_tools"

SPEED_OF_LIGHT = 1 / (8.854e-12 * 4 * np.pi * 1e-7) ** 0.5

ARRAY_KEYS = [
    "E in locations", "H in locations", "ms in locations",
    "E times picture", "H times picture", "ms time picture",
]
TOLERANCE = {"fp32": 1e-4, "fp64": 1e-9}


def co_pt_inputs(preset: str, steps: int) -> dict:
    """Co_Pt inputs, same formulas as the original run_material.py presets."""
    wavelength = 784e-9
    safety_start = 3
    if preset == "full_width":
        fwhm, dz, time_interval = 1.1e-12, 2e-9 / 8, 7
        ms0 = 3e5 * np.array([1, 0, 0, 0, 0, -0.0001]).reshape((2, 3))
    else:
        fwhm, dz, time_interval = 0.22e-12, 2e-9 / 2, 6
        ms0 = 3e5 * np.array([1, 0, 0, 0, 0, -0.00001]).reshape((2, 3))
    dt = dz / (2 * SPEED_OF_LIGHT)
    magnetic_steps = 10e-9 // dz
    conductive_steps = 2e-9 // dz
    sigma_pulse = ((0.5 / np.log(2)) ** 0.5) * fwhm

    if preset == "full_width":
        z1 = 8
        z_end = int(z1 + magnetic_steps + conductive_steps)
        z_indexes = [z1, int(z1 + magnetic_steps), z_end]
        save_loc = np.array([4] + z_indexes)
    else:
        z1 = int(2.15 * sigma_pulse / dt)
        z_end = int(z1 + magnetic_steps + conductive_steps)
        z_indexes = [z1, int(z1 + magnetic_steps), z_end]
        save_loc = np.array([4] + z_indexes + [0, z_end + 1])
        save_loc[3:5] = save_loc[2:4]
        save_loc[2] = int(z1 + magnetic_steps - 1)

    save_time = np.linspace(steps * 0.5, steps * 0.9, 5).astype(int)
    return dict(
        z_indices=z_indexes, eps_r=[1, 1], conductivity=[1.7e7, 7e6],
        damping=0.025, M0=ms0, max_E=2.4e5, phase_diff=np.pi / 2,
        pulse_width=fwhm, save_locations=save_loc, save_time_indices=save_time,
        dt=dt, total_time_steps=steps, wavelength=wavelength,
        safety_factor=safety_start, time_interval=time_interval,
    )


def load_module(name: str, path: Path):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def relative_diff(mine, ref) -> float:
    mine = np.asarray(mine)
    ref = np.asarray(ref)
    if mine.shape != ref.shape:
        raise ValueError(f"shape mismatch {mine.shape} vs {ref.shape}")
    scale = np.max(np.abs(ref))
    return float(np.max(np.abs(mine - ref)) / scale) if scale > 0 else float(np.max(np.abs(mine)))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--partner-repo", type=Path,
                        default=REPO_ROOT.parent / "llg_framework_cpu_version_updated")
    parser.add_argument("--preset", choices=["full_width", "rt"], default="full_width")
    parser.add_argument("--precision", choices=["fp32", "fp64"], default="fp64")
    parser.add_argument("--steps", type=int, default=30000)
    parser.add_argument("--state", choices=["open", "close", "both"], default="both")
    args = parser.parse_args()

    partner_file = {
        ("fp32", "full_width"): "simulation_gpu.py",
        ("fp32", "rt"): "simulation_gpu_rt.py",
        ("fp64", "full_width"): "simulation_gpu_fp64.py",
        ("fp64", "rt"): "simulation_gpu_rt_fp64.py",
    }[(args.precision, args.preset)]

    mine = load_module("mine_simulation_gpu", GPU_TOOLS / "simulation_gpu.py")
    partner = load_module("partner_backend", args.partner_repo.resolve() / partner_file)

    inputs = co_pt_inputs(args.preset, args.steps)
    states = ["open", "close"] if args.state == "both" else [args.state]
    tolerance = TOLERANCE[args.precision]
    all_ok = True

    for state in states:
        kwargs = dict(inputs, closed_system=(state == "close"))
        mine_kwargs = dict(kwargs, precision=args.precision)

        t0 = time.perf_counter()
        ours = mine.simulation(**mine_kwargs)
        t_ours = time.perf_counter() - t0
        t0 = time.perf_counter()
        ref = partner.simulation(**kwargs)
        t_ref = time.perf_counter() - t0

        print(f"\n=== Co_Pt_{state}  preset={args.preset}  precision={args.precision}  "
              f"steps={args.steps}  (ours {t_ours:.1f}s, partner {t_ref:.1f}s)")
        for key in ARRAY_KEYS:
            a, b = np.asarray(ours[key]), np.asarray(ref[key])
            if key == "E times picture":
                a, b = a[:, :-1], b[:, :-1]
            finite = bool(np.all(np.isfinite(a)))
            rel = relative_diff(a, b)
            ok = finite and rel <= tolerance
            all_ok &= ok
            print(f"  {key:<18} rel diff {rel:10.3e}   {'ok' if ok else 'FAIL'}"
                  f"{'' if finite else '  (non-finite values!)'}")

    print("\nALL OK" if all_ok else "\nSOME ARRAYS DIFFER")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
