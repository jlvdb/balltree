from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension(
        name="balltree.balltree",
        sources=[
            "cython/wrapper.pyx",
            "src/point.c",
            "src/pointbuffers.c",
            "src/ballnode.c",
            "src/ballnode_query.c",
            "src/balltree.c",
            "src/balltree_serialize.c",
        ],
        extra_compile_args=[
            "-O3",
            "-ffast-math",
            "-march=native",
            "-funroll-all-loops",
        ]
    ),
]

setup(
    name="balltree",
    author="Jan Luca van den Busch",
    version="0.1.0",
    ext_modules=cythonize(ext_modules, language_level="3"),
    include_dirs=["include", np.get_include()],
)
