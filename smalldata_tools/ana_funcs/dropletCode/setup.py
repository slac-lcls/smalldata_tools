import numpy
from distutils.core import setup
from Cython.Build import cythonize

setup(
        ext_modules = cythonize("droplet.pyx"),
        include_dirs = [numpy.get_include()],
        #ext_modules = cythonize("convolve.pyx")
)
#python setup.py build_ext --inplace
