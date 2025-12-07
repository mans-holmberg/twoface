from setuptools import setup, Extension
import numpy

try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

ext_src = "c_src/_twoface.pyx" if USE_CYTHON else "c_src/_twoface.c"

extensions = [
    Extension(
        "twoface._twoface",
        sources=[ext_src, "c_src/twoface.c"],
        include_dirs=["c_src", numpy.get_include()],
        extra_compile_args=["-O3", "-ffast-math"],
    )
]

if USE_CYTHON:
    extensions = cythonize(extensions)

setup(name="twoface",
    version="0.1.0",
    author="Måns Holmberg",
    packages=["twoface"],
    description="Transit lightcurves for asymmetric planets transiting spotted stars",
    ext_modules=extensions,
    install_requires=["numpy"],
    zip_safe=False,
)
