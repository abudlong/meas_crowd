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

"""Compare two sets of pipeline settings on tiles spread across one detector.

Each (settings, tile) pair is one task, so the fits run as a Slurm array
(``examples/usdf_validate.sh``); ``compare`` then writes a summary table and one PNG per fit::

    python examples/validate_settings.py list                     # the tasks
    python examples/validate_settings.py fit --task 3 --out DIR   # one task (Slurm sets --task)
    python examples/validate_settings.py compare --out DIR        # table + plots

Edit ``SETTINGS`` and ``POINTS`` below to test other settings or other parts of the detector.
"""

import argparse
import glob
import itertools
import os
import sys

import numpy as np

# "old" = CROWDED_FIELD_CONFIG as used in the whole-detector run; "new" = the corner-tile winner
SETTINGS = {
    "old": {},
    "new": dict(psfStackIsoFwhm=1.0, psfStackMinSnr=30.0, psfStackN=80,
                ellipTest=False, keepExtended=False),
}
# one tile around each point (detector px): corners, centre, an edge
POINTS = {
    "bottomLeft": (200, 200), "bottomRight": (3900, 150),
    "topLeft": (200, 3800), "topRight": (3900, 3800),
    "centre": (2040, 2000), "bottomMid": (2040, 150),
}
TASKS = list(itertools.product(SETTINGS, POINTS))


def _envInt(name, default):
    v = os.environ.get(name)
    return int(v) if v not in (None, "") else default


def fitOne(task, args):
    """Fit the tile of one (settings, point) task."""
    import lsst.meas.crowd as crowd
    sname, pname = TASKS[task]
    x, y = POINTS[pname]
    print(f"task {task}/{len(TASKS)}: settings {sname!r} {SETTINGS[sname]}, tile at {pname} ({x}, {y})")
    cfg = crowd.makeConfig(crowd.CROWDED_FIELD_CONFIG, **SETTINGS[sname])
    crowd.fitTiles(args.visit, args.detector, cfg, outDir=os.path.join(args.out, sname, pname),
                   region=(x, x + 1, y, y + 1), repo=args.repo, collection=args.collection,
                   psfSource="visit")


def metrics(path, png):
    """Return the comparison numbers of one tile file (and write its plot)."""
    import lsst.meas.crowd as crowd
    import matplotlib.pyplot as plt
    res = crowd.plotTile(path, out=png)
    plt.close("all")
    cat = crowd.loadTile(path)["catalog"]
    cat = cat[np.asarray(cat["owned"], bool)]
    snr = np.asarray(cat["psfSnr"], float)
    ext = np.asarray(cat["extended"], bool)
    good = ~ext & (snr > 20)
    apPsf = np.nanmedian(np.asarray(cat["apMag"], float)[good] - np.asarray(cat["psfMag"], float)[good])
    scatter = sep = np.nan
    if "gaiaMatch" in cat.colnames:
        m = np.asarray(cat["gaiaMatch"], bool) & (snr >= 50)
        dm = np.asarray(cat["gaiaDmagRel"], float)[m]
        dm = dm[np.isfinite(dm)]
        if dm.size:
            scatter = 1.4826 * np.median(np.abs(dm - np.median(dm)))
        if m.any() and "gaiaSep" in cat.colnames:
            sep = np.nanmedian(np.asarray(cat["gaiaSep"], float)[m])
    return dict(robustSigma=res["residualSigma"], nSources=len(cat), nExtended=int(ext.sum()),
                apMinusPsf=apPsf, gaiaScatter=scatter, gaiaSepMas=sep)


def compare(args):
    """Write the summary table (text and CSV) and one PNG per fitted tile."""
    import matplotlib
    matplotlib.use("Agg")
    figDir = os.path.join(args.out, "figures")
    os.makedirs(figDir, exist_ok=True)
    head = (f"{'tile':12s} {'set':4s} {'robust σ':>9s} {'sources':>8s} {'ext':>5s} "
            f"{'ap-PSF':>8s} {'Gaia sc.':>9s} {'Gaia sep':>9s}")
    lines, csv = [head], ["tile,settings,robustSigma,nSources,nExtended,apMinusPsf,gaiaScatter,gaiaSepMas"]
    for pname in POINTS:
        for sname in SETTINGS:
            files = sorted(glob.glob(os.path.join(args.out, sname, pname, "tile_*.fits")))
            if not files:
                lines.append(f"{pname:12s} {sname:4s}  -- no tile file (see logs/task_*.log)")
                continue
            r = metrics(files[0], os.path.join(figDir, f"{pname}_{sname}.png"))
            lines.append(f"{pname:12s} {sname:4s} {r['robustSigma']:9.3f} {r['nSources']:8d} "
                         f"{r['nExtended']:5d} {r['apMinusPsf']:+8.3f} {r['gaiaScatter']:9.3f} "
                         f"{r['gaiaSepMas']:9.1f}")
            csv.append(",".join([pname, sname] + [f"{r[k]}" for k in
                                                  ("robustSigma", "nSources", "nExtended",
                                                   "apMinusPsf", "gaiaScatter", "gaiaSepMas")]))
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(args.out, "validation_summary.txt"), "w") as fh:
        fh.write(text + "\n\nrobust σ: residual/σ over the tile (1 = noise); ap-PSF: median aperture - "
                 "PSF mag of point sources (S/N > 20); Gaia sc.: robust scatter of PSF mag - G "
                 "(S/N >= 50); Gaia sep: median distance to the Gaia position [mas]\n")
    with open(os.path.join(args.out, "validation_summary.csv"), "w") as fh:
        fh.write("\n".join(csv) + "\n")
    print(f"\nwrote {args.out}/validation_summary.txt, .csv and {figDir}/<tile>_<settings>.png")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="print the tasks")
    f = sub.add_parser("fit", help="fit one task's tile")
    f.add_argument("--task", type=int, default=_envInt("SLURM_ARRAY_TASK_ID", 0))
    f.add_argument("--visit", type=int, default=2025062100433)
    f.add_argument("--detector", type=int, default=144)
    f.add_argument("--repo", default="dp2_prep")
    f.add_argument("--collection", default="LSSTCam/runs/DRP/DP2")
    f.add_argument("--out", required=True)
    c = sub.add_parser("compare", help="summary table and plots")
    c.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "list":
        for i, (s, p) in enumerate(TASKS):
            print(f"{i:3d}  {s:4s}  {p:12s} {POINTS[p]}")
    elif args.cmd == "fit":
        import matplotlib
        matplotlib.use("Agg")
        fitOne(args.task, args)
    else:
        compare(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
