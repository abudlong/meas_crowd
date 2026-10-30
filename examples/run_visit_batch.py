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

"""Run meas_crowd on cutouts of one visit from a terminal (no notebook) and save everything.

    python examples/run_visit_batch.py --visit 2025062100433 --detector 144 \
        --cutouts 1000,1000,300 2500,2000,300 --out results_term

Every figure is written to <out>/figures/ as PNG (nothing is shown on screen); catalogues and
images go to <out>/ as FITS.  Use --repo to point at a Butler repository path or alias.
"""
import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import lsst.meas.crowd as crowd  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--visit", type=int, required=True)
    ap.add_argument("--detector", type=int, required=True)
    ap.add_argument("--cutouts", nargs="+", required=True, help="x0,y0,size (detector px)")
    ap.add_argument("--repo", default="dp2_prep", help="Butler repo path or alias")
    ap.add_argument("--collection", default="LSSTCam/runs/DRP/DP2")
    ap.add_argument("--out", default="results_term")
    ap.add_argument("--pad", type=int, default=20)
    ap.add_argument("--max-outer", type=int, default=3)
    ap.add_argument("--max-inner", type=int, default=15)
    ap.add_argument("--n-jobs", type=int, default=1, help="cutouts fitted in parallel")
    args = ap.parse_args()

    figdir = os.path.join(args.out, "figures")
    os.makedirs(figdir, exist_ok=True)
    counter = {"n": 0, "tag": "setup"}

    def saveAll(*a, **k):                      # replaces plt.show: save every open figure
        for num in plt.get_fignums():
            counter["n"] += 1
            plt.figure(num).savefig(os.path.join(figdir, f"{counter['n']:03d}_{counter['tag']}.png"),
                                    dpi=110, bbox_inches="tight")
        plt.close("all")
    plt.show = saveAll

    print(f"meas_crowd driver {crowd.pipeline.__version__}; C++ kernels: {crowd.HAVE_CXX}")
    cfg = crowd.makeConfig(crowd.CROWDED_FIELD_CONFIG)
    sess = crowd.prepareSession(visit=args.visit, detector=args.detector, cutouts=args.cutouts,
                                pad=args.pad, psfSource="visit", repo=args.repo,
                                collection=args.collection)
    crowd.attachGaia(sess, repo=args.repo)

    t0 = time.perf_counter()
    runs = crowd.runSessionGaia(sess, cfg, outDir=args.out, maxOuter=args.max_outer,
                                maxInner=args.max_inner, psfPolishCycles=3,
                                plotFinal=False, nJobs=args.n_jobs)
    print(f"fit: {time.perf_counter() - t0:.0f} s")

    for name, run in runs.items():
        counter["tag"] = name
        crowd.plotFinalComparison(run)
        plt.show()
        crowd.plotResidualGray(run)
        plt.show()
        phot = crowd.measurePhotometry(run)
        phot.write(os.path.join(args.out, f"{name}_photometry.fits"), overwrite=True)
        crowd.plotSourceCutouts(run, how="worst", n=8, phot=phot,
                                pdfPath=os.path.join(args.out, f"{name}_worst_sources.pdf"))
        crowd.diagnoseResidual(run)
        plt.show()
        crowd.starResidualTable(run)
        plt.show()
    print(f"figures in {figdir}/, catalogues and images in {args.out}/")


if __name__ == "__main__":
    main()
