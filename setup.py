from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

# run with "export CFLAGS="-I /usr/local/lib/python2.7/site-packages/numpy/core/include $CFLAGS""
# and then "python setup.py build_ext --inplace"

# setup(
#     ext_modules=cythonize("vdP.pyx"),
#     include_dirs=[numpy.get_include()]
# )
# setup(
#     ext_modules=[
#         Extension("vdP_heinrich2.pyx", ["vdP_heinrich2.c"],
#                   include_dirs=[numpy.get_include()]),
#     ],
# )

extensions = [
    Extension("cython_func", ["cython_func.pyx"],
        include_dirs = [numpy.get_include()]),
]
setup(ext_modules = cythonize(extensions),
)
