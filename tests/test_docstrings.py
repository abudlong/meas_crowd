# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Every public function, class and method has an LSST-style (numpydoc) docstring.

Uses numpydoc's validator with the checks Rubin packages use (e.g. daf_butler): all except
See Also / Examples / extended summary, and named return values are allowed (DM style)."""

import ast
import glob
import importlib
import inspect
import os
import unittest

import lsstStubs

lsstStubs.install()

try:
    from numpydoc.validate import validate
except ImportError:
    validate = None

MODULES = ["adaptive", "detector", "fastKernels", "fitLib", "gaiaCatalog", "photometry",
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


@unittest.skipIf(validate is None, "numpydoc not installed")
class FunctionReferenceTestCase(unittest.TestCase):
    def testFunctionReferenceUpToDate(self):
        """doc/function_reference.md matches the docstrings (regenerate it if this fails)."""
        import importlib.util
        import warnings
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(here, "..", "doc", "make_function_reference.py")
        spec = importlib.util.spec_from_file_location("make_function_reference", path)
        gen = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gen)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            text = gen.build()
        with open(gen.OUT) as fh:
            on_disk = fh.read()
        self.assertTrue(text == on_disk, "doc/function_reference.md is out of date: run "
                        "`python doc/make_function_reference.py` from the repository root")


class UnitTestDocTestCase(unittest.TestCase):
    def testUnitTestDocListsEveryTest(self):
        """doc/unit_tests.md describes every test method of tests/test_*.py."""
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "..", "doc", "unit_tests.md")) as fh:
            doc = fh.read()
        missing = []
        for path in sorted(glob.glob(os.path.join(here, "test_*.py"))):
            with open(path) as fh:
                tree = ast.parse(fh.read())
            for cls in (c for c in tree.body if isinstance(c, ast.ClassDef)):
                for m in cls.body:
                    if isinstance(m, ast.FunctionDef) and m.name.startswith("test"):
                        if f"`{cls.name}.{m.name}`" not in doc:
                            missing.append(f"{os.path.basename(path)}: {cls.name}.{m.name}")
        self.assertEqual(missing, [], "add these to doc/unit_tests.md:\n" + "\n".join(missing))


if __name__ == "__main__":
    unittest.main()
