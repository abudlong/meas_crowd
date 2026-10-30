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


class LicenceHeaderTestCase(unittest.TestCase):
    def testEverySourceFileHasTheLicenceHeader(self):
        """Every Python, C++ and shell file starts with the LSST GPLv3 header."""
        root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
        skip = {os.path.join("python", "lsst", "__init__.py"),
                os.path.join("python", "lsst", "meas", "__init__.py")}     # namespace stubs
        missing = []
        for sub in ("python", "tests", "examples", "bin.src", "doc", "src", "include", "."):
            pattern = "*" if sub == "." else os.path.join("**", "*")
            for path in glob.glob(os.path.join(root, sub, pattern), recursive=sub != "."):
                rel = os.path.relpath(path, root)
                if (rel in skip or not rel.endswith((".py", ".cc", ".h", ".sh"))
                        or ".tests" in rel or "build" + os.sep in rel):
                    continue
                with open(path) as fh:
                    head = fh.read(2000)
                if ("This file is part of meas_crowd." not in head
                        or "GNU General Public License" not in head):
                    missing.append(rel)
        self.assertEqual(sorted(set(missing)), [], "add the licence header to these files")


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
