from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension(
        name="balltree.balltree",
        sources=["cython/wrapper.pyx", "src/point.c", "src/ballnode.c", "src/balltree.c"],
        extra_compile_args=["-O3", "-march=native", "-funroll-all-loops", "-ffast-math"]
    ),
]

setup(
    name="balltree",
    author="Jan Luca van den Busch",
    version="0.0.1",
    ext_modules=cythonize(ext_modules, language_level="3"),
    include_dirs=["include", np.get_include()],
)
