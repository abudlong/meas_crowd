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

"""Tests that the old snake_case names still work after the switch to camelCase."""

import unittest
import warnings

import lsstStubs

lsstStubs.install()

import lsst.meas.crowd as crowd  # noqa: E402
from lsst.meas.crowd import pipeline, session  # noqa: E402
from lsst.meas.crowd._compat import camelKeys, toCamel  # noqa: E402


class SnakeCaseCompatTestCase(unittest.TestCase):
    def testToCamel(self):
        self.assertEqual(toCamel("run_session_gaia"), "runSessionGaia")
        self.assertEqual(toCamel("_fit_one"), "_fitOne")
        self.assertEqual(toCamel("PIPELINE_DEFAULTS"), "PIPELINE_DEFAULTS")
        self.assertEqual(toCamel("report"), "report")

    def testOldFunctionNames(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.assertIs(crowd.run_session_gaia, crowd.runSessionGaia)
            self.assertIs(pipeline.run_cutout_gaia, pipeline.runCutoutGaia)
        self.assertEqual(len(w), 2)
        self.assertTrue(all(issubclass(x.category, FutureWarning) for x in w))
        with self.assertRaises(AttributeError):
            pipeline.no_such_function

    def testOldConfigKeys(self):
        self.assertEqual(camelKeys({"stack_kwargs": {"n_use": 5}}), {"stackKwargs": {"nUse": 5}})
        cfg = session.makeConfig(deep=True, stack_kwargs=dict(n_use=5))
        self.assertEqual(cfg["stackKwargs"]["nUse"], 5)
        P = pipeline._params(dict(seed_threshold_start=9.0, n_centroid_rounds=2))
        self.assertEqual(P["gaiaThresholdStart"], 9.0)
        self.assertEqual(P["nCentroidRounds"], 2)


if __name__ == "__main__":
    unittest.main()
