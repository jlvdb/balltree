from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._balltree import BallTree, default_leafsize


def sgn(val: ArrayLike) -> ArrayLike:
    return np.where(val == 0, 1.0, np.sign(val))


def radec_to_xy(radec: NDArray[np.float64]) -> NDArray[np.float64]:
    radec = np.atleast_2d(radec)
    xy = np.empty_like(radec)
    xy[:, 0] = radec[:, 0]
    xy[:, 1] = np.sin(radec[:, 1])
    return xy


def xy_to_radec(xy: NDArray[np.float64]) -> NDArray[np.float64]:
    radec = np.empty_like(xy)
    radec[:, 0] = xy[:, 0]
    radec[:, 1] = np.arcsin(xy[:, 1])
    return radec


def radec_to_xyz(radec: NDArray[np.float64]) -> NDArray[np.float64]:
    radec = np.atleast_2d(radec)
    ra = radec[:, 0]
    dec = radec[:, 1]
    cos_dec = np.cos(dec)

    xyz = np.empty((len(ra), 3))
    xyz[:, 0] = np.cos(ra) * cos_dec
    xyz[:, 1] = np.sin(ra) * cos_dec
    xyz[:, 2] = np.sin(dec)
    return xyz


def xyz_to_radec(xyz: NDArray[np.float64]) -> NDArray[np.float64]:
    xyz = np.atleast_2d(xyz)
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    r_d2 = np.sqrt(x * x + y * y)
    r_d3 = np.sqrt(x * x + y * y + z * z)

    radec = np.empty((len(x), 2))
    x_normed = np.ones_like(x)  # fallback for zero-division, arccos(1)=0.0
    np.divide(x, r_d2, where=r_d2 > 0.0, out=x_normed)
    radec[:, 0] = np.arccos(x_normed) * sgn(y) % (2.0 * np.pi)
    radec[:, 1] = np.arcsin(x / r_d3)
    return radec


def angle_to_radius(angle: float) -> float:
    return 2.0 * np.sin(angle / 2.0)


def radius_to_angle(radius: float) -> float:
    return 2.0 * np.arcsin(radius / 2.0)


class AngularTree:
    _tree: BallTree

    def __init__(
        self,
        radec: NDArray[np.float64],
        weight: NDArray[np.float64] | None = None,
        leafsize: int = default_leafsize,
    ) -> None:
        xyz = radec_to_xyz(radec)
        self._tree = BallTree(xyz, weight, leafsize=leafsize)

    @property
    def data(self) -> NDArray[np.float64]:
        data = self._tree.data
        radec = xy_to_radec(np.transpose([data["x"], data["y"], data["z"]]))

        dtype = [("ra", "f8"), ("dec", "f8"), ("weight", "f8")]
        array = np.empty(len(data), dtype=dtype)
        array["ra"] = radec[:, 0]
        array["dec"] = radec[:, 1]
        array["weight"] = data["weight"]
        return array

    @property
    def num_data(self) -> int:
        return self._tree.num_data

    @property
    def leafsize(self) -> int:
        return self._tree.leafsize

    @property
    def sum_weight(self) -> float:
        return self._tree.sum_weight

    @property
    def center(self) -> tuple(float, float, float):
        return tuple(xy_to_radec(self._tree.center)[0])

    @property
    def radius(self) -> float:
        return radius_to_angle(self._tree.radius)

    @classmethod
    def from_random(
        cls,
        ra_min: float,
        ra_max: float,
        dec_min: float,
        dec_max: float,
        size: int,
    ) -> AngularTree:
        x_min, y_min = radec_to_xy(ra_min, dec_min)
        x_max, y_max = radec_to_xy(ra_max, dec_max)
        x = np.random.uniform(x_min, x_max, size)
        y = np.random.uniform(y_min, y_max, size)
        radec = xy_to_radec(np.transpose([x, y]))
        return cls(radec)

    @classmethod
    def from_file(cls, fpath: str) -> AngularTree:
        new = AngularTree.__new__(AngularTree)
        new._tree = BallTree.from_file(fpath)
        return new

    def to_file(self, fpath: str) -> None:
        self._tree.to_file(fpath)

    def count_nodes(self) -> int:
        return self._tree.count_nodes()

    def get_node_data(self) -> NDArray:
        return self._tree.get_node_data()

    def count_radius(
        self,
        radec: NDArray[np.float64],
        angle: float,
        weight: NDArray[np.float64] | None = None,
    ) -> float:
        xyz = radec_to_xyz(radec)
        radius = angle_to_radius(angle)
        return self._tree.count_radius(xyz, radius, weight)

    def count_range(
        self,
        radec: NDArray[np.float64],
        ang_min: float,
        ang_max: float,
        weight: NDArray[np.float64] | None = None,
    ) -> float:
        xyz = radec_to_xyz(radec)
        rmin = angle_to_radius(ang_min)
        rmax = angle_to_radius(ang_max)
        return self._tree.count_range(xyz, rmin, rmax, weight)

    def dualcount_radius(
        self,
        other: AngularTree,
        angle: float,
    ) -> float:
        radius = angle_to_radius(angle)
        return self._tree.dualcount_radius(other, radius)

    def dualcount_range(
        self,
        other: AngularTree,
        ang_min: float,
        ang_max: float,
    ) -> float:
        rmin = angle_to_radius(ang_min)
        rmax = angle_to_radius(ang_max)
        return self._tree.dualcount_range(other, rmin, rmax)
