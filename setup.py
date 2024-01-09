from setuptools import setup, Extension
import numpy as np

ext_module = Extension(
    name="pyballtree",
    sources=[
        "python/balltree_wrapper.c",
        "src/point.c",
        "src/pointbuffers.c",
        "src/ballnode.c",
        "src/ballnode_stats.c",
        "src/ballnode_query.c",
        "src/balltree.c",
        "src/balltree_serialize.c",
    ],
    include_dirs=["include", np.get_include()],
    extra_compile_args=[
        "-O3",
        "-ffast-math",
        "-march=native",
        "-funroll-all-loops",
    ],
)

setup(
    name="pyballtree",
    author="Jan Luca van den Busch",
    version="0.1.0",
    description="Python wrapper for C BallTree",
    ext_modules=[ext_module],
)
