from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from balltree import coordinates as coord
from balltree._balltree import BallTree, default_leafsize


class AngularTree:
    _tree: BallTree

    def __init__(
        self,
        radec: NDArray,
        weight: NDArray | None = None,
        leafsize: int = default_leafsize,
    ) -> None:
        xyz = coord.radec_to_xyz(radec)
        self._tree = BallTree(xyz, weight, leafsize=leafsize)

    @property
    def data(self) -> NDArray[np.float64]:
        data = self._tree.data
        radec = coord.xy_to_radec(np.transpose([data["x"], data["y"], data["z"]]))

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
        return tuple(coord.xy_to_radec(self._tree.center)[0])

    @property
    def radius(self) -> float:
        return coord.radius_to_angle(self._tree.radius)

    @classmethod
    def from_random(
        cls,
        ra_min: float,
        ra_max: float,
        dec_min: float,
        dec_max: float,
        size: int,
    ) -> AngularTree:
        x_min, y_min = coord.radec_to_xy(ra_min, dec_min)
        x_max, y_max = coord.radec_to_xy(ra_max, dec_max)
        x = np.random.uniform(x_min, x_max, size)
        y = np.random.uniform(y_min, y_max, size)
        radec = coord.xy_to_radec(np.transpose([x, y]))
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
        xyz = coord.radec_to_xyz(radec)
        radius = coord.angle_to_radius(angle)
        return self._tree.count_radius(xyz, radius, weight)

    def count_range(
        self,
        radec: NDArray,
        ang_min: float,
        ang_max: float,
        weight: NDArray | None = None,
    ) -> float:
        xyz = coord.radec_to_xyz(radec)
        rmin = coord.angle_to_radius(ang_min)
        rmax = coord.angle_to_radius(ang_max)
        return self._tree.count_range(xyz, rmin, rmax, weight)

    def dualcount_radius(
        self,
        other: AngularTree,
        angle: float,
    ) -> float:
        radius = coord.angle_to_radius(angle)
        return self._tree.dualcount_radius(other, radius)

    def dualcount_range(
        self,
        other: AngularTree,
        ang_min: float,
        ang_max: float,
    ) -> float:
        rmin = coord.angle_to_radius(ang_min)
        rmax = coord.angle_to_radius(ang_max)
        return self._tree.dualcount_range(other, rmin, rmax)
