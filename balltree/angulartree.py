from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from balltree import coordinates as coord
from balltree.balltree import BallTree, default_leafsize

__all__ = [
    "AngularTree",
]


class AngularTree:
    _tree: BallTree

    def __init__(
        self,
        radec: NDArray,
        weight: NDArray | None = None,
        leafsize: int = default_leafsize,
    ) -> None:
        xyz = coord.angular_to_euclidean(radec)
        self._tree = BallTree(xyz, weight, leafsize=leafsize)

    @property
    def data(self) -> NDArray[np.float64]:
        data = self._tree.data
        radec = coord.euclidean_to_angular(
            np.transpose([data["x"], data["y"], data["z"]])
        )

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
        return tuple(coord.euclidean_to_angular(self._tree.center)[0])

    @property
    def radius(self) -> float:
        center_xyz = coord.angular_to_euclidean(self.center)[0]
        data_xyz = self._tree.data
        # compute the maximum distance from the center project one the sphere
        dx = data_xyz["x"] - center_xyz[0]
        dy = data_xyz["y"] - center_xyz[1]
        dz = data_xyz["z"] - center_xyz[2]
        dist = np.sqrt(dx * dx + dy * dy + dz * dz)
        return coord.chorddist_to_angle(dist.max())

    @classmethod
    def from_random(
        cls,
        ra_min: float,
        ra_max: float,
        dec_min: float,
        dec_max: float,
        size: int,
    ) -> AngularTree:
        ((x_min, y_min),) = coord.angular_to_cylinder([ra_min, dec_min])
        ((x_max, y_max),) = coord.angular_to_cylinder([ra_max, dec_max])
        x = np.random.uniform(x_min, x_max, size)
        y = np.random.uniform(y_min, y_max, size)
        radec = coord.cylinder_to_angular(np.transpose([x, y]))
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
        radec: NDArray,
        angle: float,
        weight: NDArray | None = None,
    ) -> float:
        xyz = coord.angular_to_euclidean(radec)
        radius = coord.angle_to_chorddist(angle)
        return self._tree.count_radius(xyz, radius, weight)

    def count_range(
        self,
        radec: NDArray,
        ang_min: float,
        ang_max: float,
        weight: NDArray | None = None,
    ) -> float:
        xyz = coord.angular_to_euclidean(radec)
        rmin = coord.angle_to_chorddist(ang_min)
        rmax = coord.angle_to_chorddist(ang_max)
        return self._tree.count_range(xyz, rmin, rmax, weight)

    def dualcount_radius(
        self,
        other: AngularTree,
        angle: float,
    ) -> float:
        if not isinstance(other, self.__class__):
            raise TypeError("'other' must be of type 'AngularTree'")
        radius = coord.angle_to_chorddist(angle)
        return self._tree.dualcount_radius(other._tree, radius)

    def dualcount_range(
        self,
        other: AngularTree,
        ang_min: float,
        ang_max: float,
    ) -> float:
        if not isinstance(other, self.__class__):
            raise TypeError("'other' must be of type 'AngularTree'")
        rmin = coord.angle_to_chorddist(ang_min)
        rmax = coord.angle_to_chorddist(ang_max)
        return self._tree.dualcount_range(other._tree, rmin, rmax)
