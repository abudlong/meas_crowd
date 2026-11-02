#!/usr/bin/env python
# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Whole-detector meas_crowd runs from the command line.

    crowdDetector.py fit    VISIT DETECTOR [--out DIR] [--tile 400 --pad 30]
                            [--task-index K --n-tasks N]   (default: from SLURM_ARRAY_TASK_*)
                            [--region XMIN XMAX YMIN YMAX]  (test on part of the detector)
                            [--config stellar|crowded|FILE.json]
                            [--plots none|sample|all] [--all-source-pages]
    crowdDetector.py stitch [--out DIR]                     (mosaic + merged catalogue)
    crowdDetector.py plot   [--out DIR] [--view XMIN XMAX YMIN YMAX] [--png FILE]

DIR defaults to $MEAS_CROWD_OUT if set, else DEFAULT_OUT below.

--config picks the pipeline settings: 'stellar' (default) = CROWDED_STELLAR_FIELD_CONFIG,
'crowded' = CROWDED_FIELD_CONFIG (the settings of runs before November 2026), or a JSON file
(see makeConfig).  `fit` writes the settings it used to DIR/config.json.

--plots writes every diagnostic plot (PSF stars and stamps, PSF per pass, backgrounds, rounds,
source pages and cutouts, final comparison, ...) of the chosen tiles to DIR/plots/<tile>/:
'sample' (default) = 7 tiles (corners, centre, two edge midpoints), 'all' = every tile, 'none'.
Task 0 also writes the detector-level PSF plots to DIR/plots/detector/.

`fit` can be rerun: tiles already written are skipped.  Run `stitch` after all tasks finished
(or any time, to look at the tiles done so far).
"""

import argparse
import os
import sys

# where results go when neither --out nor $MEAS_CROWD_OUT is given
DEFAULT_OUT = ("/sdf/home/a/abudlong/project/repo-main-logs/crowdedFields/october2026/"
               "crowded_fields/8-detector144-stellarConfig")

# --config names -> preset in lsst.meas.crowd.session
CONFIGS = {"stellar": "CROWDED_STELLAR_FIELD_CONFIG", "crowded": "CROWDED_FIELD_CONFIG"}


def _envInt(name, default):
    v = os.environ.get(name)
    return int(v) if v not in (None, "") else default


def _plots(value):
    """--plots value -> fitTiles(plots=...)."""
    return value if value in ("none", "sample", "all") else [v.strip() for v in value.split(",")]


def main(argv=None):
    out = os.environ.get("MEAS_CROWD_OUT") or DEFAULT_OUT
    outHelp = f"output directory (default: {out})"
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fit", help="fit tiles (one share of them per Slurm array task)")
    f.add_argument("visit", type=int)
    f.add_argument("detector", type=int)
    f.add_argument("--out", default=out, help=outHelp)
    f.add_argument("--tile", type=int, default=400)
    f.add_argument("--pad", type=int, default=30)
    arrayMin = _envInt("SLURM_ARRAY_TASK_MIN", 0)
    f.add_argument("--task-index", type=int,
                   default=_envInt("SLURM_ARRAY_TASK_ID", arrayMin) - arrayMin)
    f.add_argument("--n-tasks", type=int, default=_envInt("SLURM_ARRAY_TASK_COUNT", 1))
    f.add_argument("--region", type=int, nargs=4, metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    f.add_argument("--repo", default="dp2_prep")
    f.add_argument("--collection", default="LSSTCam/runs/DRP/DP2")
    f.add_argument("--max-outer", type=int, default=3)
    f.add_argument("--max-inner", type=int, default=15)
    f.add_argument("--config", default="stellar",
                   help="stellar (default), crowded, or a JSON configuration file")
    f.add_argument("--plots", default="sample",
                   help="tiles to write all plots for: sample (default), all, none, or tile names "
                        "separated by commas")
    f.add_argument("--all-source-pages", action="store_true",
                   help="a per-source page for every source of a plotted tile (default 4 pages)")
    f.add_argument("--overwrite", action="store_true")

    s = sub.add_parser("stitch", help="assemble the tiles into a mosaic + catalogue")
    s.add_argument("--out", default=out, help=outHelp)

    p = sub.add_parser("plot", help="data | model | residual | residual/sigma of the mosaic")
    p.add_argument("--out", default=out, help=outHelp)
    p.add_argument("--view", type=int, nargs=4, metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    p.add_argument("--png", default=None, help="output file (default <out>/detector_view.png)")
    p.add_argument("--no-tiles", action="store_true")
    p.add_argument("--sources", action="store_true", help="circle the sources (zooms)")

    args = ap.parse_args(argv)
    print(f"output directory: {args.out}")
    import matplotlib
    matplotlib.use("Agg")
    from lsst.meas.crowd import detector

    if args.cmd == "fit":
        from lsst.meas.crowd import session
        base = getattr(session, CONFIGS[args.config]) if args.config in CONFIGS else args.config
        os.makedirs(args.out, exist_ok=True)
        cfg = session.makeConfig(base, saveTo=os.path.join(args.out, "config.json")
                                 if args.task_index == 0 else None)
        print(f"configuration: {CONFIGS.get(args.config, args.config)}")
        detector.fitTiles(args.visit, args.detector, cfg, outDir=args.out, tile=args.tile, pad=args.pad,
                          taskIndex=args.task_index, nTasks=args.n_tasks, region=args.region,
                          repo=args.repo, collection=args.collection, overwrite=args.overwrite,
                          maxOuter=args.max_outer, maxInner=args.max_inner,
                          plots=_plots(args.plots), allSourcePages=args.all_source_pages)
    elif args.cmd == "stitch":
        detector.stitchTiles(args.out)
    else:
        mosaic = detector.stitchTiles(args.out, verbose=False)
        png = args.png or os.path.join(args.out, "detector_view.png")
        st = detector.plotDetector(mosaic, view=args.view, showTiles=not args.no_tiles,
                                   sources=args.sources, out=png)
        print(f"residual/sigma robust sigma {st['residualSigma']:.3f}; wrote {png}")
        catPath = os.path.join(args.out, "detector_catalog.fits")
        if args.view is None and os.path.exists(catPath):     # photometry of the whole detector
            from astropy.table import Table
            from lsst.meas.crowd.plotting import plotGaiaComparison
            cat = Table.read(catPath)
            if "gaiaMatch" in cat.colnames:
                gpng = os.path.join(args.out, "detector_gaia_comparison.png")
                plotGaiaComparison(cat, out=gpng)
                print(f"wrote {gpng}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
