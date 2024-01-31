from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "sgn",
    "angular_to_cylinder",
    "cylinder_to_angular",
    "angular_to_euclidean",
    "euclidean_to_angular",
    "angle_to_chorddist",
    "chorddist_to_angle",
]


def require_dim(data: NDArray, ndim: int) -> None:
    data_dim = data.shape[1]
    if data_dim != ndim:
        raise ValueError(f"input must have {ndim} dimensions but got {data_dim}")


def sgn(val: ArrayLike) -> NDArray:
    return 2.0 * (val >= 0.0) - 1.0  # 1.0 if (val >= 0.0) else -1.0


def angular_to_cylinder(radec: NDArray) -> NDArray:
    """
    Convert angular coordinates in radian to cylindrical coordinates.

    ``radec`` can be of shape (2,) or (N, 2), result is guaranteed to be (N, 2).
    """
    radec = np.atleast_2d(radec)
    require_dim(radec, 2)
    xy = np.empty_like(radec, dtype=np.float64)
    xy[:, 0] = radec[:, 0]
    xy[:, 1] = np.sin(radec[:, 1])
    return xy


def cylinder_to_angular(xy: NDArray) -> NDArray:
    """
    Convert cylindrical coordinates to angular coordinates in radian.

    ``xy`` can be of shape (2,) or (N, 2), result is guaranteed to be (N, 2).
    """
    xy = np.atleast_2d(xy)
    require_dim(xy, 2)
    radec = np.empty_like(xy, dtype=np.float64)
    radec[:, 0] = xy[:, 0]
    radec[:, 1] = np.arcsin(xy[:, 1])
    return radec


def angular_to_euclidean(radec: NDArray) -> NDArray:
    """
    Convert angular coordinates in radian to Euclidean coordinates.

    ``radec`` can be of shape (2,) or (N, 2), result is guaranteed to be (N, 3).
    """
    radec = np.atleast_2d(radec)
    require_dim(radec, 2)
    ra = radec[:, 0]
    dec = radec[:, 1]
    cos_dec = np.cos(dec)

    xyz = np.empty((len(ra), 3), dtype=np.float64)
    xyz[:, 0] = np.cos(ra) * cos_dec
    xyz[:, 1] = np.sin(ra) * cos_dec
    xyz[:, 2] = np.sin(dec)
    return xyz


def euclidean_to_angular(xyz: NDArray) -> NDArray:
    """
    Convert Euclidean coordinates to angular coordinates in radian.

    ``xyz`` can be of shape (3,) or (N, 3), result is guaranteed to be (N, 2).
    """
    xyz = np.atleast_2d(xyz)
    require_dim(xyz, 3)
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    r_d2 = np.sqrt(x * x + y * y)
    r_d3 = np.sqrt(x * x + y * y + z * z)

    radec = np.empty((len(x), 2), dtype=np.float64)
    x_normed = np.ones_like(x)  # fallback for zero-division, arccos(1)=0.0
    np.divide(x, r_d2, where=r_d2 > 0.0, out=x_normed)
    radec[:, 0] = np.arccos(x_normed) * sgn(y) % (2.0 * np.pi)
    radec[:, 1] = np.arcsin(z / r_d3)
    return radec


def angle_to_chorddist(angle: ArrayLike) -> NDArray:
    """Convert great circle distance (in radian) to chord distance."""
    return 2.0 * np.sin(angle / 2.0, dtype=np.float64)


def chorddist_to_angle(chord: ArrayLike) -> NDArray:
    """Convert chord distance to great circle distance (in radian)."""
    return 2.0 * np.arcsin(chord / 2.0, dtype=np.float64)
