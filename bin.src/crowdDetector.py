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

    crowdDetector.py fit    VISIT DETECTOR --out DIR [--tile 400 --pad 30]
                            [--task-index K --n-tasks N]   (default: from SLURM_ARRAY_TASK_*)
                            [--region XMIN XMAX YMIN YMAX]  (test on part of the detector)
    crowdDetector.py stitch --out DIR                       (mosaic + merged catalogue)
    crowdDetector.py plot   --out DIR [--view XMIN XMAX YMIN YMAX] [--png FILE]

`fit` can be rerun: tiles already written are skipped.  Run `stitch` after all tasks finished
(or any time, to look at the tiles done so far).
"""

import argparse
import os
import sys


def _env_int(name, default):
    v = os.environ.get(name)
    return int(v) if v not in (None, "") else default


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fit", help="fit tiles (one share of them per Slurm array task)")
    f.add_argument("visit", type=int)
    f.add_argument("detector", type=int)
    f.add_argument("--out", required=True)
    f.add_argument("--tile", type=int, default=400)
    f.add_argument("--pad", type=int, default=30)
    array_min = _env_int("SLURM_ARRAY_TASK_MIN", 0)
    f.add_argument("--task-index", type=int,
                   default=_env_int("SLURM_ARRAY_TASK_ID", array_min) - array_min)
    f.add_argument("--n-tasks", type=int, default=_env_int("SLURM_ARRAY_TASK_COUNT", 1))
    f.add_argument("--region", type=int, nargs=4, metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    f.add_argument("--repo", default="dp2_prep")
    f.add_argument("--collection", default="LSSTCam/runs/DRP/DP2")
    f.add_argument("--max-outer", type=int, default=3)
    f.add_argument("--max-inner", type=int, default=15)
    f.add_argument("--overwrite", action="store_true")

    s = sub.add_parser("stitch", help="assemble the tiles into a mosaic + catalogue")
    s.add_argument("--out", required=True)

    p = sub.add_parser("plot", help="data | model | residual | residual/sigma of the mosaic")
    p.add_argument("--out", required=True)
    p.add_argument("--view", type=int, nargs=4, metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    p.add_argument("--png", default=None, help="output file (default <out>/detector_view.png)")
    p.add_argument("--no-tiles", action="store_true")
    p.add_argument("--sources", action="store_true", help="circle the sources (zooms)")

    args = ap.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    from lsst.meas.crowd import detector

    if args.cmd == "fit":
        detector.fit_tiles(args.visit, args.detector, out_dir=args.out, tile=args.tile, pad=args.pad,
                           task_index=args.task_index, n_tasks=args.n_tasks, region=args.region,
                           repo=args.repo, collection=args.collection, overwrite=args.overwrite,
                           max_outer=args.max_outer, max_inner=args.max_inner)
    elif args.cmd == "stitch":
        detector.stitch_tiles(args.out)
    else:
        mosaic = detector.stitch_tiles(args.out, verbose=False)
        png = args.png or os.path.join(args.out, "detector_view.png")
        st = detector.plot_detector(mosaic, view=args.view, show_tiles=not args.no_tiles,
                                    sources=args.sources, out=png)
        print(f"residual/sigma robust sigma {st['residual_sigma']:.3f}; wrote {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
