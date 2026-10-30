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

"""Draw the pipeline flow diagrams of doc/pipeline_guide.md as SVG.

Writes ``pipeline_flow.svg`` (the steps at a glance) and ``pipeline_flow_detailed.svg`` (each
step with what it does) next to this script.  Edit the text below and rerun::

    python doc/images/make_flow_diagrams.py
"""

import os
from xml.sax.saxutils import escape

HERE = os.path.dirname(os.path.abspath(__file__))
FONT = "Helvetica, Arial, sans-serif"
INK, QUIET, EDGE = "#202124", "#5f6368", "#9aa0a6"
ACCENT, ACCENT_FILL, BAND = "#3b6fc4", "#e8f0fe", "#f3f4f6"


class Svg:
    """A minimal SVG canvas."""

    def __init__(self, width, height, title):
        self.w, self.h, self.items = width, height, []
        self.title = title

    def text(self, x, y, s, size=12, weight="normal", color=INK, anchor="start", style=""):
        """Add one line of text (baseline at y)."""
        self.items.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}" '
                          f'fill="{color}" text-anchor="{anchor}"{style}>{escape(s)}</text>')

    def rect(self, x, y, w, h, fill="#ffffff", stroke=EDGE, width=1.25, rx=8):
        """Add a rounded rectangle."""
        self.items.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
                          f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')

    def path(self, d, dashed=False, arrow=True):
        """Add a connector (with an arrowhead at its end)."""
        dash = ' stroke-dasharray="5 4"' if dashed else ""
        end = ' marker-end="url(#arrow)"' if arrow else ""
        self.items.append(f'<path d="{d}" fill="none" stroke="{EDGE}" stroke-width="1.4"{dash}{end}/>')

    def diamond(self, cx, cy, hw, hh):
        """Add a decision diamond."""
        pts = f"{cx - hw},{cy} {cx},{cy - hh} {cx + hw},{cy} {cx},{cy + hh}"
        self.items.append(f'<polygon points="{pts}" fill="#ffffff" stroke="{EDGE}" stroke-width="1.25"/>')

    def box(self, x, y, w, title, lines=(), main=False, center=True, size=12, titleSize=14):
        """Add a box with a bold title and body lines; return its height."""
        h = 34 + 17 * len(lines) + (6 if lines else -6)
        self.rect(x, y, w, h, fill=ACCENT_FILL if main else "#ffffff",
                  stroke=ACCENT if main else EDGE, width=2 if main else 1.25)
        tx, anchor = (x + w / 2, "middle") if center else (x + 14, "start")
        self.text(tx, y + 23, title, size=titleSize, weight="bold", anchor=anchor)
        for i, s in enumerate(lines):
            self.text(tx, y + 43 + 17 * i, s, size=size, color=QUIET, anchor=anchor)
        return h

    def write(self, path):
        """Write the SVG file."""
        head = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" '
                f'width="{self.w}" height="{self.h}" font-family="{FONT}" role="img" '
                f'aria-label="{escape(self.title)}">\n<title>{escape(self.title)}</title>\n'
                '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
                f'markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" '
                f'fill="{EDGE}"/></marker></defs>\n'
                f'<rect width="{self.w}" height="{self.h}" fill="#ffffff"/>\n')
        with open(path, "w") as fh:
            fh.write(head + "\n".join(self.items) + "\n</svg>\n")


def minimal():
    """The steps at a glance."""
    s = Svg(760, 486, "meas_crowd pipeline flow")
    s.text(24, 32, "Each cutout loops fit, then refit, until sources and PSF stop changing",
           size=16, weight="bold")
    for x, t, b in ((24, "1. Setup", "cutouts, bad-pixel mask, PSF"),
                    (280, "2. Gaia catalogue", "starting positions from DR3"),
                    (536, "3. Preparation", "mask, PSF wings, Gaia offset")):
        s.box(x, 52, 200, t, [b])
    s.path("M224 80H280")
    s.path("M480 80H536")
    s.path("M636 108V122H290V136")
    # outer pass
    s.rect(24, 136, 532, 160, fill=BAND)
    s.text(40, 160, "4. Outer pass (PSF and background fixed)", size=14, weight="bold")
    s.box(40, 176, 156, "Inner loop", ["detect, fit, clean up;", "lower the threshold",
                                       "until the residual is clean"], main=True, size=11)
    s.box(212, 176, 156, "Clean-up", ["remove halo / ring fakes,", "polish positions,",
                                      "split pairs, fit galaxies"], size=11)
    s.box(384, 176, 156, "Noise + PSF check", ["add missing photon noise", "to the variance;",
                                               "measure the PSF error"], size=11)
    s.path("M196 222H212")
    s.path("M368 222H384")
    # decision
    s.path("M556 222H582")
    s.diamond(642, 222, 60, 42)
    s.text(642, 218, "5. Changed", size=12, anchor="middle")
    s.text(642, 234, "less than 1%?", size=12, anchor="middle")
    s.text(712, 212, "yes", size=12, color=QUIET, anchor="middle")
    s.path("M702 222H722V398")
    s.text(652, 290, "no", size=12, color=QUIET)
    s.path("M642 264V340H452")
    # refit
    s.box(172, 312, 280, "6. Refit for the next pass", ["background (tilted plane), PSF from stars"])
    s.path("M312 312V296")
    # outputs and detector
    s.box(536, 398, 200, "7. Outputs", ["photometry catalogue with", "Gaia comparison, files"])
    s.box(24, 398, 490, "8. Whole detector: fitTiles, then stitchTiles",
          ["400 px tiles with 30 px pad, then one mosaic + one catalogue"])
    s.path("M80 398V296", dashed=True)
    s.text(90, 344, "each tile", size=11, color=QUIET)
    s.text(90, 359, "runs 1 to 7", size=11, color=QUIET)
    s.write(os.path.join(HERE, "pipeline_flow.svg"))


def detailed():
    """Each step with what it does."""
    s = Svg(760, 1238, "meas_crowd pipeline flow, detailed")
    s.text(24, 32, "meas_crowd pipeline, step by step", size=16, weight="bold")
    s.text(24, 52, "Function names in brackets; configuration options in the text (e.g. bgMode).",
           size=12, color=QUIET)
    L, W = 24, 712
    y = 70
    steps = [
        ("1. Setup  (prepareSession)", [
            "Load the visit_image; cut out each box grown by pad px, before any background subtraction.",
            "Bad pixels: SAT, CR, BAD, INTRP, EDGE.  PSF: the visit's own (psfSource=\"visit\") or a",
            "new detector-level PCA PSF (\"refit\", which also subtracts a detector background)."]),
        ("2. Gaia catalogue  (attachGaia, selectGaiaSources)", [
            "Gaia DR3 over the detector (Butler reference catalogue or the ESA archive), moved to the",
            "visit epoch with proper motions, cached on disk.  Starting positions: the Gaia stars in",
            "the cutout and up to 1.5 FWHM outside it, brightest first, duplicates merged."]),
        ("3. Preparation  (start of runAdaptivePipeline)", [
            "Grow the saturated mask by 2 px; mark saturated-star zones.  Extend the PSF kernel from",
            "41 to 81 px with a Moffat wing and smooth its outer part.  Remove a global offset",
            "between the Gaia positions and the image."]),
    ]
    for title, lines in steps:
        h = s.box(L, y, W, title, lines, center=False)
        s.path(f"M380 {y + h}V{y + h + 22}")
        y += h + 22
    # outer pass container
    top = y
    idx = len(s.items)                       # the container is drawn first, sized at the end
    s.text(L + 16, top + 26, "4. Outer pass: the PSF and background stay fixed",
           size=14, weight="bold")
    y = top + 42
    inner = [
        "a. Detect peaks in the residual above the current threshold (5σ when starting from Gaia,",
        "    then lower).",
        "b. Add them, skipping edges, duplicates and positions already rejected as fakes.",
        "c. Joint non-negative fit of all sources; bright and saturated stars get correction terms",
        "    (shift, width, core and halo profile) so their residual is clean.",
        "d. Refine centroids; merge close pairs; drop sources below minSnr or with a bad shape",
        "    (too broad ones are kept as small galaxies).",
        "e. Residual check: are there significantly more positive than negative residual peaks?",
        "    Yes: lower the threshold (down to 3σ) and repeat.  No (or nothing changes): stop.",
    ]
    h = s.box(L + 16, y, 640,
              "Inner loop  (runInnerAdaptive): find and fit sources until the residual is clean",
              inner, main=True, center=False)
    s.path(f"M{L + 656} {y + h - 22}H{L + 672}V{y + 22}H{L + 657}")
    s.text(L + 684, y + h / 2, "repeat", size=11, color=QUIET, anchor="middle",
           style=f' transform="rotate(-90 {L + 684} {y + h / 2})"')
    s.path(f"M{L + 336} {y + h}V{y + h + 16}")
    y += h + 16
    h = s.box(L + 16, y, 640, "Clean-up at the end of the pass", [
        "Remove faint non-Gaia sources that are a bright star's halo or ring structure (halo and",
        "ring checks; these positions are not detected again).  Polish all positions together.",
        "Test unresolved pairs (add a second source) and small galaxies (elliptical model)."],
        center=False)
    s.path(f"M{L + 336} {y + h}V{y + h + 16}")
    y += h + 16
    h = s.box(L + 16, y, 640, "Noise and PSF-error check", [
        "Photon noise: if the bright stars' residual is pixel-scale noise, the variance plane is",
        "missing photon noise; add k × model to the variance for the later fits (autoPoisson).",
        "PSF error: measure ε, the PSF model error as a fraction of the model per pixel; plots show",
        "residual / sqrt(variance + (ε × model)²) (modelErrMeasure)."], center=False)
    y += h + 16
    s.rect(L, top, 690, y - top, fill=BAND)
    s.items.insert(idx, s.items.pop())
    # decision
    s.path(f"M380 {y}V{y + 24}")
    cy = y + 24 + 46
    s.diamond(380, cy, 125, 46)
    s.text(380, cy - 4, "5. Converged?", size=13, weight="bold", anchor="middle")
    s.text(380, cy + 13, "sources and PSF changed < 1%", size=11, color=QUIET, anchor="middle")
    s.text(248, cy + 4, "at least 2, at most maxOuter passes", size=11, color=QUIET, anchor="end")
    s.text(509, cy - 8, "no", size=12, color=QUIET)
    s.path(f"M505 {cy}H530")
    h6 = s.box(530, cy - 54, 206, "6. Refit for next pass", [
        "Background: tilted plane on", "the raw cutout, sources masked", "(bgMode).  PSF: residual",
        "stack of up to 25 clean stars,", "wings from bright stars."], center=False, size=11)
    s.path(f"M736 {cy - 54 + h6 / 2}H748V{top + 40}H{L + 690}")
    s.text(390, cy + 66, "yes", size=12, color=QUIET)
    s.path(f"M380 {cy + 46}V{cy + 82}")
    y = cy + 82
    h = s.box(L, y, W, "7. Outputs  (runCutoutGaia → CutoutRun)", [
        "run.photometry: PSF and aperture photometry with errors and the Gaia comparison",
        "(doc/catalog_columns.md).  run.outerHistory: every pass (images, model, PSF, background).",
        "With outDir: photometry.fits / .csv, catalog.fits, images.fits, psf_kernels.fits, summary."],
        center=False)
    y += h + 22
    s.path(f"M380 {y - 22}V{y}", dashed=True)
    s.text(390, y - 8, "repeated for every tile of a detector", size=11, color=QUIET)
    h = s.box(L, y, W, "8. Whole detector  (fitTiles → stitchTiles → plotDetector)", [
        "400 px tiles, each fitted with 30 px padding by steps 1 to 7 (a Slurm array shares the",
        "tiles; finished tiles are skipped on a restart) → tile_<name>.fits.  stitchTiles: one",
        "detector_mosaic.fits (data, model, residual, ...) and detector_catalog.fits (each source once)."],
        center=False)
    s.h = y + h + 24
    s.write(os.path.join(HERE, "pipeline_flow_detailed.svg"))


if __name__ == "__main__":
    minimal()
    detailed()
    print("wrote pipeline_flow.svg and pipeline_flow_detailed.svg")
