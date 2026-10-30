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

"""Run a sample: one cutout plus a small detector region, with every plot of the user guide.

From a terminal, with the stack and meas_crowd set up (``setup -k -r .``)::

    python examples/run_sample.py --out sample_out
    python examples/run_sample.py --out sample_refit --psf-source refit
    python examples/run_sample.py --out sample_out --skip-detector      # cutout only, faster

Every figure is saved to ``<out>/figures/NNN_<step>.png`` (nothing is shown on screen); the
photometry, the source pages (PNG per page), the tile files and the mosaic go to ``<out>/``.  The same
steps as notebook cells are in ``examples/run_sample.ipynb``.

Steps (numbers as in the user guide):

1. run the pipeline on one cutout (and report what each pass did)
2. PSF candidate / used stars over the visit image, and the PSF model
3. initial background and the background of every outer pass
4. data | model | residual per outer pass and per inner iteration
5. per-source pages linked to coloured boxes on the data / model / residual
6. individual tiles of a small detector region
7. the stitched data | model | residual of that region (or the whole detector)
"""
import argparse
import glob
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import lsst.meas.crowd as crowd  # noqa: E402

# source catalogues that carry the visit PSF's calib_psf_candidate / calib_psf_used flags; the
# name depends on the pipeline version (list them with butler.registry.queryDatasetTypes("*star*"))
PSF_CATALOGS = ("single_visit_star_unstandardized", "single_visit_star_footprints", "src")


def load_visit(repo, collection, visit, detector, instrument="LSSTCam",
               visit_image_collections="LSSTCam/runs/*"):
    """Return the Butler, the visit image and its data id."""
    from lsst.daf.butler import Butler
    butler = Butler(repo, collections=collection)
    refs = butler.query_datasets("visit_image", collections=visit_image_collections, find_first=False,
                                 where=f"instrument='{instrument}' AND visit={visit} AND "
                                       f"detector={detector}")
    if not refs:
        raise RuntimeError(f"no visit_image for visit {visit}, detector {detector}")
    return butler, butler.get(refs[0]), refs[0].dataId


def visit_psf_catalog(butler, data_id, collections="LSSTCam/runs/*"):
    """Return the first source catalogue with PSF flags that exists, or None."""
    for name in PSF_CATALOGS:
        try:
            cat = butler.get(name, dataId=data_id, collections=collections)
            print(f"PSF stars of the visit PSF: from {name} ({len(cat)} rows)")
            return cat
        except Exception as e:
            print(f"  ({name}: {type(e).__name__})")
    print("no source catalogue with PSF flags found -- skipping the visit-PSF star plot")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--visit", type=int, default=2025062100433)
    ap.add_argument("--detector", type=int, default=144)
    ap.add_argument("--cutout", default="1000,1000,300", help="x0,y0,size (detector px)")
    ap.add_argument("--pad", type=int, default=20)
    ap.add_argument("--psf-source", choices=("visit", "refit"), default="visit")
    ap.add_argument("--repo", default="dp2_prep", help="Butler repo path or alias")
    ap.add_argument("--collection", default="LSSTCam/runs/DRP/DP2")
    ap.add_argument("--max-outer", type=int, default=3)
    ap.add_argument("--max-inner", type=int, default=15)
    ap.add_argument("--region", type=int, nargs=4, default=(800, 1599, 800, 1599),
                    metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
                    help="detector region for the tile run (default: 2 x 2 tiles of 400 px)")
    ap.add_argument("--whole-detector", action="store_true",
                    help="fit every tile (slow: ~100 tiles; use the Slurm array instead)")
    ap.add_argument("--skip-detector", action="store_true", help="only the cutout (steps 1-5)")
    ap.add_argument("--out", default="sample_out")
    args = ap.parse_args()

    figdir = os.path.join(args.out, "figures")
    os.makedirs(figdir, exist_ok=True)
    state = {"n": 0, "tag": "setup"}

    def save_all(*a, **k):                      # replaces plt.show: save every open figure
        for num in plt.get_fignums():
            state["n"] += 1
            plt.figure(num).savefig(os.path.join(figdir, f"{state['n']:03d}_{state['tag']}.png"),
                                    dpi=110, bbox_inches="tight")
        plt.close("all")
    plt.show = save_all

    def step(tag):
        plt.show()                              # flush the previous step's figures
        state["tag"] = tag
        print(f"\n==== {tag}")

    print(f"meas_crowd driver {crowd.pipeline.__version__}; C++ kernels: {crowd.HAVE_CXX}")

    # ---- 1. run the pipeline on one cutout
    step("1_run")
    butler, visit_image, data_id = load_visit(args.repo, args.collection, args.visit, args.detector)
    cfg = crowd.make_config(crowd.CROWDED_FIELD_CONFIG)
    sess = crowd.prepare_session(args.visit, args.detector, [args.cutout], pad=args.pad,
                                 psf_source=args.psf_source, visit_image=visit_image)
    crowd.attach_gaia(sess, visit_image=visit_image, repo=args.repo)
    t0 = time.perf_counter()
    runs = crowd.run_session_gaia(sess, cfg, out_dir=args.out, max_outer=args.max_outer,
                                  max_inner=args.max_inner, plot_final=False)
    if not runs:
        raise SystemExit("the cutout fit failed -- see the traceback printed above")
    name, run = next(iter(runs.items()))
    print(f"fit: {time.perf_counter() - t0:.0f} s")
    crowd.report(run)
    crowd.plot_final_comparison(run)
    crowd.plot_residual_gray(run)
    # made by the fit, with the Gaia comparison columns; saved as
    # <out>/visit<V>_det<D>/<cutout>/photometry.fits (and .csv)
    phot = run.photometry
    if getattr(sess, "gaia", None) is not None:
        crowd.plot_gaia_comparison(phot)

    # ---- 2. PSF stars and the PSF model
    step("2_psf")
    if args.psf_source == "refit":
        crowd.show_psf_stars(visit_image, crowd.psf_star_positions(session=sess), psf=sess.psf)
    else:
        cat = visit_psf_catalog(butler, data_id)
        if cat is not None:
            crowd.show_psf_stars(visit_image, crowd.psf_star_positions(catalog=cat), psf=sess.psf)
    hist = run.outer_history
    for k in range(len(hist) - 1):              # the in-cutout PSF update after each pass
        crowd.show_psf_stars(hist[k]["exposure"], crowd.psf_star_positions(run=run, outer=k),
                             psf=hist[k + 1]["psf"])
    crowd.fitLib.plot_psf_models_by_outer_pass(hist)

    # ---- 3. backgrounds
    step("3_background")
    crowd.show_backgrounds(run, initial=sess.bgImage)

    # ---- 4. model and residual per round
    step("4_rounds")
    crowd.show_rounds(run, level="outer")
    for k in range(len(hist)):
        crowd.show_rounds(run, level="inner", outer=k, max_rows=6)
    crowd.fitLib.plot_outer_loop_summary(hist)
    for k in range(len(hist) - 1):              # sources lost / gained between passes, and why
        crowd.compare_passes(run, a=k, b=k + 1)

    # ---- 5. source pages linked to boxes on the data / model / residual
    step("5_sources")
    crowd.show_linked_source_pages(run, source_view="all", sources_per_page=6, max_pages=1)
    crowd.show_linked_source_pages(run, max_pages=None,          # one PNG per page, memory flat
                                   out_dir=os.path.join(args.out, f"{name}_source_pages"))
    crowd.plot_source_cutouts(run, how="worst", n=8, phot=phot)

    if args.skip_detector:
        step("done")
        print(f"\nfigures in {figdir}/, files in {args.out}/")
        return

    # ---- 6. tiles of a detector region
    step("6_tiles")
    det_dir = os.path.join(args.out, "detector_run")
    if args.psf_source == "refit":              # the refit modified the visit image: reload it
        visit_image = None
    crowd.fit_tiles(args.visit, args.detector, cfg, out_dir=det_dir, tile=400, pad=30,
                    region=None if args.whole_detector else tuple(args.region),
                    repo=args.repo, collection=args.collection, psf_source=args.psf_source,
                    visit_image=visit_image)
    for path in sorted(glob.glob(os.path.join(det_dir, "tile_*.fits"))):
        crowd.plot_tile(path)

    # ---- 7. stitched data / model / residual
    step("7_detector")
    mosaic = crowd.stitch_tiles(det_dir)
    crowd.plot_detector(mosaic)
    if not args.whole_detector:
        x0, x1, y0, y1 = args.region
        crowd.plot_detector(mosaic, view=(x0, x0 + 199, y0, y0 + 199), binning=1, sources=True)
    step("done")
    print(f"\nfigures in {figdir}/, files in {args.out}/ and {det_dir}/")


if __name__ == "__main__":
    main()
