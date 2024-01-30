from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from balltree import coordinates as coord
from balltree.balltree import BallTree, default_leafsize

__all__ = [
    "AngularTree",
]


class AngularTree(BallTree):
    def __init__(
        self,
        radec: NDArray,
        weight: NDArray | None = None,
        leafsize: int = default_leafsize,
    ) -> None:
        xyz = coord.angular_to_euclidean(radec)
        super().__init__(xyz, weight, leafsize=leafsize)

    @property
    def data(self) -> NDArray[np.float64]:
        data = super().data
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
    def center(self) -> tuple(float, float, float):
        return tuple(coord.euclidean_to_angular(super().center)[0])

    @property
    def radius(self) -> float:
        center = coord.angular_to_euclidean(self.center)[0]
        radec_flat = self.data.view("f8")
        shape = (self.num_data, -1)
        xyz = coord.angular_to_euclidean(radec_flat.reshape(shape)[:, :-1])
        # compute the maximum distance from the center project one the sphere
        diff = xyz - center[np.newaxis, :]
        dist = np.sqrt(np.sum(diff**2, axis=1))
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

    def brute_radius(
        self,
        radec: NDArray,
        angle: float,
        weight: NDArray | None = None,
    ) -> float:
        xyz = coord.angular_to_euclidean(radec)
        radius = coord.angle_to_chorddist(angle)
        return super().brute_radius(xyz, radius, weight)

    def brute_range(
        self,
        radec: NDArray,
        ang_min: float,
        ang_max: float,
        weight: NDArray | None = None,
    ) -> float:
        xyz = coord.angular_to_euclidean(radec)
        rmin = coord.angle_to_chorddist(ang_min)
        rmax = coord.angle_to_chorddist(ang_max)
        return super().brute_range(xyz, rmin, rmax, weight)

    def count_radius(
        self,
        radec: NDArray,
        angle: float,
        weight: NDArray | None = None,
    ) -> float:
        xyz = coord.angular_to_euclidean(radec)
        radius = coord.angle_to_chorddist(angle)
        return super().count_radius(xyz, radius, weight)

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
        return super().count_range(xyz, rmin, rmax, weight)

    def dualcount_radius(
        self,
        other: AngularTree,
        angle: float,
    ) -> float:
        if not isinstance(other, self.__class__):
            raise TypeError("'other' must be of type 'AngularTree'")
        radius = coord.angle_to_chorddist(angle)
        return super().dualcount_radius(other, radius)

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
        return super().dualcount_range(other, rmin, rmax)
