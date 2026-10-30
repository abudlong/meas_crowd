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

"""Old snake_case names (used before the switch to Rubin camelCase) mapped to the current ones."""

__all__ = ["toCamel", "camelKeys", "moduleGetattr"]

import re
import warnings


def toCamel(name):
    """Return the camelCase form of a snake_case name.

    Parameters
    ----------
    name : `str`
        A name such as ``run_session_gaia``; leading and trailing underscores and UPPER_CASE
        constants are kept as they are.

    Returns
    -------
    camel : `str`
        The name in camelCase (``runSessionGaia``).
    """
    core = name.strip("_")
    if "_" not in core or core.isupper():
        return name
    lead = name[:len(name) - len(name.lstrip("_"))]
    trail = name[len(name.rstrip("_")):]
    camel = re.sub(r"_+([A-Za-z0-9])", lambda m: m.group(1).upper(), core)
    return lead + camel + trail


def camelKeys(cfg):
    """Return a copy of a configuration with old snake_case keys renamed to camelCase.

    Parameters
    ----------
    cfg : `dict`
        A configuration, possibly from a JSON file saved by an older version; nested dicts
        (``stackKwargs`` and the like) are converted too.

    Returns
    -------
    cfg : `dict`
        The configuration with camelCase keys.
    """
    if not isinstance(cfg, dict):
        return cfg
    return {(toCamel(k) if isinstance(k, str) else k): camelKeys(v) for k, v in cfg.items()}


def moduleGetattr(globs):
    """Make a module ``__getattr__`` that finds the camelCase name for an old snake_case one.

    Parameters
    ----------
    globs : `dict`
        The module's ``globals()``.

    Returns
    -------
    getattr : callable
        The module ``__getattr__`` (PEP 562); it warns with `FutureWarning` when an old name
        is used.
    """
    def __getattr__(name):
        camel = toCamel(name)
        if camel != name and camel in globs:
            warnings.warn(f"{name} has been renamed {camel}", FutureWarning, stacklevel=2)
            return globs[camel]
        raise AttributeError(f"module {globs['__name__']!r} has no attribute {name!r}")
    return __getattr__
