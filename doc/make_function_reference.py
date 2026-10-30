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

"""Write doc/function_reference.md from the package's numpydoc docstrings.

Run from the repository root after changing a public function or its docstring::

    python doc/make_function_reference.py

Works with or without the Rubin stack (without it, the test stand-ins in ``tests/lsstStubs.py``
are used to import the package).  ``tests/test_docstrings.py`` fails when the page is out of
date.
"""

import importlib
import inspect
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "function_reference.md")

# (module, one-line role), in reading order: the driver and user-facing modules first, the
# original notebook libraries last
MODULES = [
    ("pipeline", "the driver: run the pipeline on cutouts from Gaia starting positions, inner / outer loops, "
                 "PSF refinement, diagnostics"),
    ("session", "load a visit, cut out regions, set up the PSF; configuration; saving runs"),
    ("gaiaCatalog", "the Gaia catalogue and the starting positions of each cutout"),
    ("photometry", "PSF and aperture photometry, errors, the Gaia comparison"),
    ("plotting", "final data / model / residual plots and residual diagnostics"),
    ("stagePlots", "plots of the intermediate stages and pass-to-pass comparison"),
    ("sourcePlots", "per-source cutouts, cross-sections and magnitudes"),
    ("detector", "whole-detector runs: tiles, stitching, mosaic and tile plots"),
    ("fastKernels", "fast PSF stamps (C++ with a NumPy fallback)"),
    ("fitLib", "the original cutout_fit_lib: joint fit, detection, PSF and background tools, "
               "plots"),
    ("adaptive", "the original cutout_adaptive: the adaptive inner-loop decisions"),
]

START_HERE = [
    ("session", "make_config"), ("session", "prepare_session"), ("gaiaCatalog", "attach_gaia"),
    ("pipeline", "run_session_gaia"), ("pipeline", "report"),
    ("photometry", "measure_photometry"), ("plotting", "plot_final_comparison"),
    ("plotting", "plot_gaia_comparison"), ("stagePlots", "show_rounds"),
    ("stagePlots", "show_linked_source_pages"), ("stagePlots", "compare_passes"),
    ("detector", "fit_tiles"), ("detector", "stitch_tiles"), ("detector", "plot_detector"),
]


def _import_package():
    """Import the package, with the test stand-ins when the Rubin stack is missing."""
    sys.path.insert(0, os.path.join(ROOT, "python"))
    try:
        import lsst.afw  # noqa: F401
    except ImportError:
        sys.path.insert(0, os.path.join(ROOT, "tests"))
        import lsstStubs
        lsstStubs.install()


def _anchor(text):
    """Return GitHub's anchor for a heading."""
    return re.sub(r"[^a-z0-9_\- ]", "", text.lower()).replace(" ", "-")


def _cell(lines):
    """Join docstring lines into one Markdown table cell."""
    text = " ".join(s.strip() for s in lines if s.strip())
    return text.replace("|", "\\|") or "–"


def _para(lines):
    """Join docstring lines into paragraphs (blank lines separate them)."""
    out, cur = [], []
    for s in list(lines) + [""]:
        if s.strip():
            cur.append(s.strip())
        elif cur:
            out.append(" ".join(cur))
            cur = []
    return "\n\n".join(out)


def _public(mod):
    """Return the module's public functions and classes in source order."""
    items = [(n, o) for n, o in vars(mod).items()
             if not n.startswith("_") and getattr(o, "__module__", None) == mod.__name__
             and (inspect.isfunction(o) or inspect.isclass(o))
             and o.__name__ == n]                  # skip aliases kept for old names
    return sorted(items, key=lambda t: inspect.getsourcelines(t[1])[1])


class _Name:
    """A default value shown by its name in the source (instead of ``<... at 0x...>``)."""

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name


def _signature(name, obj):
    """Return ``name(signature)``, with object defaults shown by their global name."""
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        return f"{name}(...)"
    glob = getattr(obj, "__globals__", {})
    params = []
    for p in sig.parameters.values():
        d = p.default
        if d is not p.empty and " at 0x" in repr(d):
            nm = next((k for k, v in glob.items() if v is d), type(d).__name__ + "()")
            p = p.replace(default=_Name(nm))
        params.append(p.replace(annotation=p.empty))     # annotations differ between versions
    return f"{name}{sig.replace(parameters=params, return_annotation=sig.empty)}"


def _doc_block(doc, lines):
    """Append the summary, parameters, returns and raises of a parsed docstring."""
    summary = _para(doc["Summary"])
    if summary:
        lines += [summary, ""]
    ext = _para(doc["Extended Summary"])
    if ext:
        lines += [ext, ""]
    for section, head in (("Parameters", "Parameter"), ("Attributes", "Attribute"),
                          ("Returns", "Returns"), ("Yields", "Yields")):
        rows = doc[section]
        if not rows:
            continue
        lines += [f"| {head} | Type | Description |", "| --- | --- | --- |"]
        for p in rows:
            nm = f"`{p.name}`" if p.name else "–"
            lines.append(f"| {nm} | {_cell([p.type]) if p.type else '–'} | {_cell(p.desc)} |")
        lines.append("")
    if doc["Raises"]:
        lines.append("Raises: " + "; ".join(
            f"`{p.type or p.name}` ({_cell(p.desc)})" for p in doc["Raises"]))
        lines.append("")


def build():
    """Return the text of the function reference page.

    Returns
    -------
    text : `str`
        The Markdown page.
    """
    from numpydoc.docscrape import ClassDoc, FunctionDoc
    _import_package()
    # the pipeline swaps some fitLib functions for faster ones while it runs; list the originals
    importlib.import_module("lsst.meas.crowd.pipeline").disable_speedups()
    mods = [(name, role, importlib.import_module(f"lsst.meas.crowd.{name}"))
            for name, role in MODULES]
    items = {name: _public(mod) for name, _, mod in mods}
    n_fun = sum(1 for v in items.values() for _, o in v if inspect.isfunction(o))
    n_cls = sum(1 for v in items.values() for _, o in v if inspect.isclass(o))

    lines = ["# meas_crowd: function reference", "",
             "Every public function and class of `lsst.meas.crowd`: what it does, its "
             "parameters (with types and defaults) and what it returns.  Generated from the "
             "docstrings by `doc/make_function_reference.py`; do not edit by hand (run "
             "`python doc/make_function_reference.py` after changing a docstring).  How the "
             "pipeline fits together is in [pipeline_guide.md](pipeline_guide.md); the output "
             "catalogue in [catalog_columns.md](catalog_columns.md).", "",
             f"{n_fun} functions and {n_cls} classes in {len(mods)} modules.  Everything listed "
             "here can be used as `crowd.<name>` after `import lsst.meas.crowd as crowd`, or "
             "from its module (`from lsst.meas.crowd import photometry`).  Names starting with "
             "`_` are internal helpers and are not listed.", "",
             "## Start here", "",
             "| Function | What it is for |", "| --- | --- |"]
    for m, n in START_HERE:
        obj = getattr(importlib.import_module(f"lsst.meas.crowd.{m}"), n)
        doc = FunctionDoc(obj)
        lines.append(f"| [`{m}.{n}`](#{_anchor(f'{m}.{n}')}) | {_cell(doc['Summary'])} |")
    lines += ["", "## Modules", "", "| Module | Role | Public functions / classes |",
              "| --- | --- | --- |"]
    for name, role, mod in mods:
        lines.append(f"| [`{name}`](#{_anchor(name)}) | {role} | {len(items[name])} |")
    lines.append("")

    for name, role, mod in mods:
        lines += [f"## {name}", "", f"`lsst.meas.crowd.{name}`: {role}.", ""]
        mdoc = inspect.getdoc(mod) or ""
        if mdoc:
            lines += [mdoc.split("\n\n")[0].replace("\n", " "), ""]
        lines += ["| Name | Summary |", "| --- | --- |"]
        for n, o in items[name]:
            d = ClassDoc(o) if inspect.isclass(o) else FunctionDoc(o)
            kind = " (class)" if inspect.isclass(o) else ""
            lines.append(f"| [`{n}`](#{_anchor(f'{name}.{n}')}){kind} | {_cell(d['Summary'])} |")
        lines.append("")
        for n, o in items[name]:
            if inspect.isclass(o):
                lines += [f"### {name}.{n}", "", f"Class.  `{_signature(n, o)}`", ""]
                _doc_block(ClassDoc(o, config={"show_class_members": False}), lines)
                for mn, mo in sorted(((k, v) for k, v in vars(o).items()
                                      if inspect.isfunction(v) and not k.startswith("_")),
                                     key=lambda t: inspect.getsourcelines(t[1])[1]):
                    lines += [f"#### {name}.{n}.{mn}", "", f"`{_signature(mn, mo)}`", ""]
                    _doc_block(FunctionDoc(mo), lines)
            else:
                lines += [f"### {name}.{n}", "", f"`{_signature(n, o)}`", ""]
                _doc_block(FunctionDoc(o), lines)
    return "\n".join(lines).rstrip() + "\n"


def main():
    """Write the page."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        text = build()
    with open(OUT, "w") as fh:
        fh.write(text)
    print(f"wrote {OUT} ({text.count(chr(10))} lines)")


if __name__ == "__main__":
    main()
