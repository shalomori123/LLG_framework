"""
Run Maxwell-LLG GPU simulations.

Examples:
    python sim_run.py                                   # full_width, all materials, open+close
    python sim_run.py --preset rt --materials Co_Pt
    python sim_run.py --materials Co_Pt --state close --steps 30000 --output-dir quick_check
"""
import argparse
from pathlib import Path

from tqdm import tqdm

import sim_config as cfg
import simulation_gpu as sim


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--preset", choices=list(cfg.PRESETS), default="full_width")
    parser.add_argument("--materials", nargs="+", choices=cfg.ALL_MATERIALS, default=cfg.ALL_MATERIALS)
    parser.add_argument("--state", choices=["open", "close", "both"], default="both")
    parser.add_argument("--precision", choices=["fp64", "fp32"], default="fp64",
                        help="fp32 is faster on consumer GPUs but can't resolve the open/close difference")
    parser.add_argument("--steps", type=int, default=None,
                        help="truncate the run to this many time steps (quick checks)")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="default: simulation_result/<preset>")
    parser.add_argument("--overwrite", action="store_true", help="rerun cases whose bundle already exists")
    return parser.parse_args()


def main():
    args = parse_args()
    preset = cfg.build_preset(args.preset, args.steps)
    output_dir = args.output_dir or Path("simulation_result") / args.preset
    output_dir.mkdir(parents=True, exist_ok=True)

    states = ["open", "close"] if args.state == "both" else [args.state]
    runs = [(material, state) for material in args.materials for state in states]
    print(f"preset={args.preset}  precision={args.precision}  grid={preset['grid_size']} cells  "
          f"steps={preset['simulation_time']}  "
          f"runs={len(runs)}  -> {output_dir}")

    for material, state in tqdm(runs, desc="cases", unit="case"):
        bundle_path = output_dir / f"{material}_{state}.simdata"
        if sim.bundle_complete(bundle_path) and not args.overwrite:
            continue
        case = cfg.material_case(material)
        sim.simulation(
            preset["z_indexes"], case["eps_r"], case["conductivity"], case["damping"], preset["ms0"],
            preset["max_E"], preset["phase_diff"], preset["fwhm"],
            preset["save_locations"], preset["save_time"], cfg.STATES[state],
            preset["dt"], preset["simulation_time"],
            wavelength=preset["wavelength"], safety_factor=preset["safety_start"],
            time_interval=preset["time_interval"], precision=args.precision,
            output_path=bundle_path,
        )


if __name__ == "__main__":
    main()
