"""HGCC driver for the heavy Fig. 1 steps (clustering, DE, association, pathway association) on ALL kept
tumor cells of one dataset. Run from this folder:

    python 1_heavy.py --dataset Xenium_5K_LC            # one dataset
    python 1_heavy.py --dataset 3                       # by index into config.DATASETS (SLURM array)
    python 1_heavy.py --dataset Xenium_5K_LC --quick    # 20k-cell development run -> output/2_exploration/<ds>_quick/
    python 1_heavy.py --dataset all --steps 3,3b       # subset of steps; earlier outputs are loaded from cache

Outputs go to ../../output/2_exploration/<ds>/ and are consumed by fig1_exploration.ipynb after `make out-pull`.
"""

import argparse
import time
import warnings

import config as C
import pipeline as P

warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--dataset", default="all", help="dataset name, index into config.DATASETS, or 'all'")
ap.add_argument("--steps", default="3,3b,4,7,8", help="comma-separated subset of 3,3b,4,7,8")
ap.add_argument("--quick", action="store_true", help="subsample to config.QUICK_N_CELLS cells, fewer bootstraps")
ap.add_argument("--n-cells", type=int, default=None, help="override the quick-mode cell budget")
ap.add_argument("--force", action="store_true", help="recompute even if cached outputs exist")
args = ap.parse_args()

if args.dataset == "all":
    datasets = C.DATASETS
elif args.dataset.isdigit():
    datasets = [C.DATASETS[int(args.dataset)]]
else:
    datasets = [args.dataset]
steps = tuple(s.strip() for s in args.steps.split(","))

for ds in datasets:
    t0 = time.time()
    print(f"========== {ds} (steps {','.join(steps)}{', quick' if args.quick else ''}) ==========", flush=True)
    ctx = P.prepare(ds, quick=args.quick, n_cells=args.n_cells)
    P.step0_qc(ctx)
    P.run_heavy(ctx, steps=steps, force=args.force)
    print(f"========== {ds} done in {(time.time() - t0) / 60:.1f} min ==========", flush=True)
