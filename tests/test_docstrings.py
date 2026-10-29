# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Every public function, class and method has an LSST-style (numpydoc) docstring.

Uses numpydoc's validator with the checks Rubin packages use (e.g. daf_butler): all except
See Also / Examples / extended summary, and named return values are allowed (DM style)."""

import importlib
import inspect
import unittest

import lsstStubs

lsstStubs.install()

try:
    from numpydoc.validate import validate
except ImportError:
    validate = None

MODULES = ["adaptive", "detector", "fastKernels", "fitLib", "gaiaSeeds", "photometry",
           "pipeline", "plotting", "session", "sourcePlots", "stagePlots"]
IGNORE = {"SA01", "EX01", "SS06", "GL01", "ES01", "RT02", "SS05", "GL02"}


def _public_objects(mod):
    for name, obj in vars(mod).items():
        if name.startswith("_") or getattr(obj, "__module__", None) != mod.__name__:
            continue
        if inspect.isfunction(obj):
            yield f"{mod.__name__}.{name}"
        elif inspect.isclass(obj):
            yield f"{mod.__name__}.{name}"
            for mname, mobj in vars(obj).items():
                if inspect.isfunction(mobj) and not mname.startswith("_"):
                    yield f"{mod.__name__}.{name}.{mname}"


@unittest.skipIf(validate is None, "numpydoc not installed")
class DocstringTestCase(unittest.TestCase):
    def testNumpydoc(self):
        problems = []
        for m in MODULES:
            mod = importlib.import_module(f"lsst.meas.crowd.{m}")
            for qual in _public_objects(mod):
                errs = [f"{c} {msg}" for c, msg in validate(qual)["errors"] if c not in IGNORE]
                if errs:
                    problems.append(f"{qual}: " + "; ".join(errs))
        self.assertEqual(problems, [], "\n" + "\n".join(problems))


if __name__ == "__main__":
    unittest.main()
