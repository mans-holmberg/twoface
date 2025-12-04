from setuptools import setup, Extension
from Cython.Build import cythonize, build_ext
import numpy
import os

from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        "twoface._twoface",
        sources=[
            "c_src/_twoface.pyx",  # Cython wrapper
            "c_src/twoface.c",  # C code
        ],
        include_dirs=["c_src"],
        extra_compile_args=["-O3", "-ffast-math"],
    )
]

setup(name="twoface",
    version="0.1.0",
    author="Måns Holmberg",
    packages=["twoface"],
    description="Transit light-curve modeling for asymmetric planets transiting spotted stars.",
    ext_modules=cythonize(extensions),
    include_dirs=[numpy.get_include()],
    install_requires=["numpy"],
    zip_safe=False,
    cmdclass={'build_ext': build_ext},
)