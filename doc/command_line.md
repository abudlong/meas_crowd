# meas_crowd: whole-detector runs from the command line (USDF)

How to fit a whole detector from a terminal at USDF, without a notebook: a quick test on a few
tiles, the full run as a Slurm array, how to follow it, and what it writes.  The Python side of
the same run (`fitTiles`, `stitchTiles`, `plotDetector`) is in the
[pipeline guide](pipeline_guide.md).

## What a detector run does

1. **Fit.** The detector is cut into tiles of at most 400 px whose interiors cover it exactly:
   11 × 10 = **110 tiles** of 370–371 × 400 px for an LSSTCam science detector (4072 × 4000 px),
   named by their corner and size (`tile_x0740_y0800_371x400.fits`).  Each tile is fitted on its
   interior grown by 30 px (so stars near a tile edge are fitted together with their neighbours)
   and written to `tile_<name>.fits`.  The Slurm array splits the tiles between its tasks (25
   tasks → 4–5 tiles each).
2. **Stitch.** The tile interiors are assembled into one detector image (`detector_mosaic.fits`)
   and the tile catalogues merged into one (`detector_catalog.fits`), each source kept once.
3. **Plot.** `detector_view.png`: data | model | residual | residual/σ of the whole detector.

The Slurm script runs the stitch and plot steps by itself once every fit task has ended.

## One-time setup

In a terminal on a USDF login node (use the same stack as for your notebooks):

```bash
source /sdf/group/rubin/sw/w_latest/loadLSST.bash
setup lsst_distrib
cd /sdf/data/rubin/user/$USER/code/build/meas_crowd     # your meas_crowd checkout
setup -k -r . && scons
```

Repeat the first three lines (and `setup -k -r .`) in every new terminal.  Run `scons` again
after updating the code.

## Where the results go

Every command takes `--out DIR`.  Without it the results go to `$MEAS_CROWD_OUT` if that is set,
else to `DEFAULT_OUT` at the top of `bin.src/crowdDetector.py`.  Each command prints the
directory it uses.  The Slurm script takes the directory from `OUT=`.

Set it once per terminal:

```bash
OUT=/sdf/home/a/abudlong/project/repo-main-logs/crowdedFields/october2026/crowded_fields/6-testScript/10.30.26_test2
```

`OUT=... command` on one line sets it for that command only; set it on its own line (as above)
to use `$OUT` in later commands.  A new folder per test (`..._test3`, ...) keeps runs apart.

## Step 1: a quick test on a few tiles

```bash
python -u bin.src/crowdDetector.py fit 2025062100433 144 --out $OUT/test_region --region 800 1599 800 1599
python bin.src/crowdDetector.py stitch --out $OUT/test_region
python bin.src/crowdDetector.py plot --out $OUT/test_region      # writes detector_view.png
```

`--region XMIN XMAX YMIN YMAX` fits only the tiles that overlap that area (here 2 × 2 tiles).
`python -u` prints every message as it happens, so you see how far it gets.

## Step 2: the whole detector as a Slurm array

```bash
OUT=$OUT VISIT=2025062100433 DETECTOR=144 PKG=$PWD bash examples/usdf_detector_array.sh
```

| Variable | Meaning | Default |
| --- | --- | --- |
| `VISIT`, `DETECTOR` | what to fit | 2025062100433, 144 |
| `NTASKS` | number of array tasks sharing the 110 tiles | 25 |
| `OUT` | output directory (created if needed) | the October 2026 test directory |
| `PKG` | your meas_crowd checkout | `$HOME/meas_crowd` |

Each fit task gets 1 CPU, 8 GB and 4 hours; the stitch job 16 GB and 30 minutes.  To change them,
or the partition / account (`milano`, `rubin:developers`; `sacctmgr show assoc user=$USER` lists
yours), edit the `sbatch` lines of `examples/usdf_detector_array.sh`.  Also check that its
`SETUP=` line sources the stack you use.

## Step 3: follow the run

```bash
squeue -u $USER                       # crowd_fit tasks; crowd_stitch waits (Dependency)
ls $OUT/tile_*.fits | wc -l           # finished tiles, out of 110
tail -f $OUT/logs/fit_0.log           # one task live (Ctrl-C stops watching, not the job)
sacct -u $USER -S today --format=JobID%20,JobName,State,Elapsed,MaxRSS | grep crowd
```

A task loads the visit image, the Gaia catalogue and the PSF before its first tile, so the first
tile files can take a while to appear.  To look at the tiles finished so far, stitch at any time
(missing tiles show in magenta):

```bash
python bin.src/crowdDetector.py stitch --out $OUT
python bin.src/crowdDetector.py plot --out $OUT --png $OUT/partial_view.png
```

## Step 4: restart or fill gaps

Run the same submit command again.  Tiles already on disk are skipped, so only missing ones
(after a `TIMEOUT`, `OUT_OF_MEMORY` or a failed task) are fitted, and a new stitch job follows.
`--overwrite` on `crowdDetector.py fit` refits everything.

## What you get

| File | Contents |
| --- | --- |
| `tile_<name>.fits` | one per tile, over its fit box: DATA, MODEL, RESIDUAL, RESIDUAL_SIGMA, VARIANCE, BAD and CATALOG (the tile's photometry; `owned` = inside this tile's interior) |
| `detector_mosaic.fits` | the whole detector: DATA, MODEL, RESIDUAL, RESIDUAL_SIGMA, VARIANCE, BAD, TILE (which tile each pixel came from, −1 = not fitted) |
| `detector_catalog.fits` | every source once, with the photometry and Gaia columns ([catalog_columns.md](catalog_columns.md)) |
| `detector_view.png` | data / model / residual / residual/σ of the whole detector |
| `logs/fit_<n>.log`, `logs/stitch.log` | one log per task and the stitch log |

DATA is the image the fit used (after the pipeline's background subtraction), RESIDUAL = DATA −
MODEL and RESIDUAL_SIGMA = RESIDUAL / √VARIANCE (NaN on bad pixels): about ±1 where the model
fits to the noise.  Tile files written before the residual images were added get them when you
next run `stitch` (nothing is refitted).

## Step 5: look at the results

```bash
# a zoom with the sources circled
python bin.src/crowdDetector.py plot --out $OUT --view 1000 1400 2000 2400 --sources --png $OUT/zoom.png
```

In Python (stack set up, or the notebook kernel):

```python
from astropy.table import Table
import lsst.meas.crowd as crowd

OUT = "/sdf/home/a/abudlong/project/repo-main-logs/crowdedFields/october2026/crowded_fields/6-testScript/10.30.26_test2"
cat = Table.read(f"{OUT}/detector_catalog.fits")
print(len(cat), "sources;", cat["gaiaMatch"].sum(), "matched to Gaia")
crowd.plotGaiaComparison(cat)

m = crowd.loadMosaic(OUT)                   # m["data"], m["model"], m["resid"], m["residSigma"]
crowd.plotDetector(m, view=(1000, 1400, 2000, 2400), binning=1, sources=True)
crowd.plotTile(f"{OUT}/tile_x0740_y0800_371x400.fits")      # any tile file (ls $OUT/tile_*)
```

`detector_mosaic.fits` opens in DS9.  `ds9 -multiframe $OUT/detector_mosaic.fits` loads every
extension (DATA, MODEL, RESIDUAL, RESIDUAL_SIGMA, ...) into its own frame; `ds9
"$OUT/detector_mosaic.fits[RESIDUAL]"` opens one.

## Command reference

```text
crowdDetector.py fit    VISIT DETECTOR [--out DIR] [--tile 400] [--pad 30]
                        [--region XMIN XMAX YMIN YMAX] [--task-index K --n-tasks N]
                        [--repo dp2_prep] [--collection LSSTCam/runs/DRP/DP2]
                        [--max-outer 3] [--max-inner 15] [--overwrite]
crowdDetector.py stitch [--out DIR]
crowdDetector.py plot   [--out DIR] [--view XMIN XMAX YMIN YMAX] [--png FILE]
                        [--no-tiles] [--sources]
```

- `--task-index` / `--n-tasks` default to the Slurm array variables, so inside an array each
  task picks its own share of the tiles.
- The command line always uses the visit's own PSF and `CROWDED_FIELD_CONFIG`.  For other
  settings, run `fitTiles` from Python with your configuration.

## Troubleshooting

| What you see | What it means / what to do |
| --- | --- |
| `mkdir: cannot create directory '.../lsst-scipipe-.../.cargo': Permission denied` at the top of every log | Printed by the stack's environment setup, not by meas_crowd; harmless. |
| A log shows nothing after that line | Python's output is buffered: the script sets `PYTHONUNBUFFERED=1` so lines appear as they are printed (add it to an older copy of the script's `SETUP=` line).  Check `squeue` / `sacct`: if the task is `RUNNING`, it is working. |
| `ls $OUT/tile_*.fits` finds nothing | Check `echo $OUT` first: set on the submit line only, `OUT` is empty in your shell afterwards.  Otherwise the first tiles may still be running. |
| Task state `TIMEOUT` or `OUT_OF_MEMORY` | Resubmit (finished tiles are kept); raise `--time` or `--mem` in the script if it happens again. |
| Task `FAILED`, or a Python traceback in a log | The last ~30 lines of that log say why. |
| A task runs for hours with 0 tiles | Possibly stuck on the Gaia catalogue: without a Butler refcat the code asks the ESA archive over the internet, which compute nodes may not reach.  Run the Step 1 test in a terminal to see where it stops. |
