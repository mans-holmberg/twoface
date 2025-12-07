from setuptools import setup, Extension
import numpy
from Cython.Build import cythonize

extensions = [
    Extension(
        "twoface._twoface",
        sources=["c_src/_twoface.pyx", "c_src/twoface.c"],
        include_dirs=["c_src", numpy.get_include()],
        extra_compile_args=["-O3", "-ffast-math"],
    )
]

setup(name="twoface",
    version="0.1.0",
    author="Måns Holmberg",
    packages=["twoface"],
    description="Transit lightcurves for asymmetric planets transiting spotted stars",
    ext_modules=cythonize(extensions),
    install_requires=["numpy"],
    zip_safe=False,
)
