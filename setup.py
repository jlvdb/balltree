from setuptools import setup, Extension
import numpy as np

ext_module = Extension(
    name="balltree",
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
        "-DSET_PYERR_STRING",  # required to propagate C errors to python
    ],
)

if __name__ == "__main__":
    setup(
        ext_modules=[ext_module],
        url="https://github.com/jlvdb/balltree.git",
    )
