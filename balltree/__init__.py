__version__ = "0.4"

from .angulartree import AngularTree
from .balltree import BallTree, default_leafsize as _df

default_leafsize = _df
"""Default leaf size for trees, defined as ``DEFAULT_LEAFSIZE`` in `balltree.h`"""


__all__ = [
    "AngularTree",
    "BallTree",
    "default_leafsize",
]
