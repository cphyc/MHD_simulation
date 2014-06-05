from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'Simple MHD simulation code',
    ext_modules = cythonize("cviewer.pyx"),
)
