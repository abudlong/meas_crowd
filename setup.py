"""Stand-alone build of the C++ extension, for use WITHOUT the Rubin stack (laptops, CI).

Inside the stack use ``setup -r . && scons`` instead.  Here::

    pip install pybind11
    python setup.py build_ext --inplace     # puts _crowdLib*.so into python/lsst/meas/crowd
"""
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext = Pybind11Extension(
    "lsst.meas.crowd._crowdLib",
    ["python/lsst/meas/crowd/_crowdLib.cc", "src/stamps.cc"],
    include_dirs=["include"],
    cxx_std=17,
    extra_compile_args=["-O3"],
)

setup(
    name="meas_crowd",
    package_dir={"": "python"},
    ext_modules=[ext],
    cmdclass={"build_ext": build_ext},
)
