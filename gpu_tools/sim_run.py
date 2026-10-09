"""
Run Maxwell-LLG GPU simulations.

Examples:
    python sim_run.py                                   # full_width, all materials, open+close
    python sim_run.py --batch                           # same, all cases in the same kernel launches
    python sim_run.py --preset rt --materials Co_Pt
    python sim_run.py --materials Co_Pt --state close --steps 30000 --output-dir quick_check
"""
import argparse
from pathlib import Path

import sim_config as cfg
import simulation_gpu as sim


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--preset", choices=list(cfg.PRESETS), default="full_width")
    parser.add_argument("--materials", nargs="+", choices=cfg.ALL_MATERIALS, default=cfg.ALL_MATERIALS)
    parser.add_argument("--state", choices=["open", "close", "both"], default="both")
    parser.add_argument("--precision", choices=["fp64", "fp32"], default="fp64",
                        help="fp32 is faster on consumer GPUs but can't resolve the open/close difference")
    parser.add_argument("--batch", type=int, nargs="?", const=0, default=1, metavar="C",
                        help="run C cases per launch (blockIdx.y); without a value C is sized from free VRAM")
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
    runs = []
    for material in args.materials:
        for state in states:
            bundle_path = output_dir / f"{material}_{state}.simdata"
            if sim.bundle_complete(bundle_path) and not args.overwrite:
                continue
            case = dict(cfg.material_case(material), closed_system=cfg.STATES[state])
            runs.append((case, bundle_path))

    material_size = preset["z_indexes"][-1] - preset["z_indexes"][0]
    batch = args.batch or sim.max_batch(preset["grid_size"], material_size, args.precision)
    batch = max(1, min(batch, len(runs) or 1))
    print(f"preset={args.preset}  precision={args.precision}  grid={preset['grid_size']} cells  "
          f"steps={preset['simulation_time']}  runs={len(runs)}  batch={batch}  -> {output_dir}")

    for first in range(0, len(runs), batch):
        group = runs[first:first + batch]
        print("running: " + ", ".join(path.stem for _, path in group))
        sim.simulation_cases(
            [case for case, _ in group], preset["z_indexes"], preset["ms0"],
            preset["max_E"], preset["phase_diff"], preset["fwhm"],
            preset["save_locations"], preset["save_time"],
            preset["dt"], preset["simulation_time"],
            wavelength=preset["wavelength"], safety_factor=preset["safety_start"],
            time_interval=preset["time_interval"], precision=args.precision,
            output_paths=[path for _, path in group],
        )


if __name__ == "__main__":
    main()
