from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from balltree import coordinates as coord
from balltree.balltree import BallTree, default_leafsize

__all__ = [
    "AngularTree",
]


class AngularTree(BallTree):
    """
    A wrapper around ``balltree.BallTree`` for using angular coordinates.

    The the method and attributes of this tree class work in angular units as
    they are commonly used in astronomy, i.e. `right ascension` and
    `declination` in radian for the first and second coordinates. Distances and
    radii are expressed as great circle distance in radians. Internally, data is
    still represented in Euclidean (x, y, z) coordinates.

    The data point(s) ``radec`` can be a numpy array of shape (2,) or (N, 2),
    or an equivalent python object. The optional ``weights`` can be a float
    or a 1-dim sequence of matching length, the optional ``leafsize``
    determines when the tree query algorithms switch from traversal to brute
    force.
    """

    def __init__(
        self,
        radec: ArrayLike,
        weight: ArrayLike | None = None,
        leafsize: int = default_leafsize,
    ) -> None:
        xyz = coord.angular_to_euclidean(radec)
        super().__init__(xyz, weight, leafsize=leafsize)

    @classmethod
    def from_random(
        cls,
        ra_min: float,
        ra_max: float,
        dec_min: float,
        dec_max: float,
        size: int,
        leafsize: int = default_leafsize
    ) -> AngularTree:
        """
        Build a new AngularTree instance from randomly generated points.
        
        The (ra, dec) coordinates are generated uniformly in the interval
        [`ra_min`, `ra_max`) and [`dec_min`, `dec_max`), respectively. `size`
        controlls the number of points generated. The optional `leafsize`
        determines when the tree query algorithms switch from traversal to brute
        force.
        """
        ((x_min, y_min),) = coord.angular_to_cylinder([ra_min, dec_min])
        ((x_max, y_max),) = coord.angular_to_cylinder([ra_max, dec_max])
        x = np.random.uniform(x_min, x_max, size)
        y = np.random.uniform(y_min, y_max, size)
        radec = coord.cylinder_to_angular(np.transpose([x, y]))
        return cls(radec, leafsize=leafsize)

    @classmethod
    def from_file(cls, path: str) -> AngularTree:
        return super().from_file(path)

    @property
    def data(self) -> NDArray:
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
    def num_points(self) -> int:
        return super().num_points

    @property
    def leafsize(self) -> int:
        return super().leafsize

    @property
    def sum_weight(self) -> float:
        return super().sum_weight

    @property
    def center(self) -> tuple(float, float):
        return tuple(coord.euclidean_to_angular(super().center)[0])

    @property
    def radius(self) -> float:
        center = coord.angular_to_euclidean(self.center)[0]
        radec_flat = self.data.view("f8")
        shape = (self.num_points, -1)
        xyz = coord.angular_to_euclidean(radec_flat.reshape(shape)[:, :-1])
        # compute the maximum distance from the center project one the sphere
        diff = xyz - center[np.newaxis, :]
        dist = np.sqrt(np.sum(diff**2, axis=1))
        return coord.chorddist_to_angle(dist.max())

    def to_file(self, path: str) -> None:
        """Store a representation of the tree instance in a binary file."""
        return super().to_file(path)

    def count_nodes(self) -> int:
        """Get a count of all nodes of the tree, including the root node."""
        return super().count_nodes()

    def get_node_data() -> NDArray:
        """
        Collect the meta data of all tree nodes in a numpy array.

        The array fields record `depth` (starting from the root node),
        `num_points`, `sum_weight`, `x`, `y`, `z` (node center) and node `radius`.

        .. Note::
            The node coordinates and radius are currently not converted to
            anlges.
        """
        return super().get_node_data()

    def brute_radius(
        self,
        radec: ArrayLike,
        angle: float,
        weight: ArrayLike | None = None,
    ) -> float:
        """
        Count neighbours within a given angle in radian using brute force.
        
        The query point(s) ``radec`` can be a numpy array of shape (2,) or (N, 2),
        or an equivalent python object. The optional ``weights`` can be a float
        or a 1-dim sequence of matching length.
        """
        xyz = coord.angular_to_euclidean(radec)
        radius = coord.angle_to_chorddist(angle)
        return super().brute_radius(xyz, radius, weight)

    def count_radius(
        self,
        radec: ArrayLike,
        angle: float,
        weight: ArrayLike | None = None,
    ) -> float:
        """
        Count neighbours within a given angle in radian using tree traversal.
        
        The query point(s) ``radec`` can be a numpy array of shape (2,) or (N, 2),
        or an equivalent python object. The optional ``weights`` can be a float
        or a 1-dim sequence of matching length.
        """
        xyz = coord.angular_to_euclidean(radec)
        radius = coord.angle_to_chorddist(angle)
        return super().count_radius(xyz, radius, weight)

    def dualcount_radius(
        self,
        other: AngularTree,
        angle: float,
    ) -> float:
        """
        Count all pairs within a given angle in radian from points in another tree.
        
        The pairs between the two trees are computed with the efficient dualtree
        algorithm.
        """
        if not isinstance(other, self.__class__):
            raise TypeError("'other' must be of type 'AngularTree'")
        radius = coord.angle_to_chorddist(angle)
        return super().dualcount_radius(other, radius)

    def brute_range(
        self,
        radec: ArrayLike,
        angles: ArrayLike,
        weight: ArrayLike | None = None,
    ) -> NDArray:
        """
        Count neighbours within a sequence of angles in radian using brute force.

        The query point(s) ``radec`` can be a numpy array of shape (2,) or (N, 2),
        or an equivalent python object. The ``angles`` must either be a float or
        monotonic sequence. The optional ``weights`` can be a float or a 1-dim
        sequence of matching length.

        Returns an array of counts. The first element contains the count of all
        neighbours ``0 <= r <= r_1``, the following values contain the
        incremental counts ``r_i-1 <= r <= r_i``.
        """
        xyz = coord.angular_to_euclidean(radec)
        radii = coord.angle_to_chorddist(angles)
        return super().brute_range(xyz, radii, weight)

    def count_range(
        self,
        radec: ArrayLike,
        angles: ArrayLike,
        weight: ArrayLike | None = None,
    ) -> NDArray:
        """
        Count neighbours within a sequence of angles in radian using tree traversal.

        The query point(s) ``radec`` can be a numpy array of shape (2,) or (N, 2),
        or an equivalent python object. The ``angles`` must either be a float or
        monotonic sequence. The optional ``weights`` can be a float or a 1-dim
        sequence of matching length.

        Returns an array of counts. The first element contains the count of all
        neighbours ``0 <= r <= r_1``, the following values contain the
        incremental counts ``r_i-1 <= r <= r_i``.
        """
        xyz = coord.angular_to_euclidean(radec)
        radii = coord.angle_to_chorddist(angles)
        return super().count_range(xyz, radii, weight)

    def dualcount_range(
        self,
        other: AngularTree,
        angles: ArrayLike,
    ) -> NDArray:
        """
        Count all pairs within a sequence of angles in radian from points in another
        tree.
        
        The pairs between the two trees are computed with the efficient dualtree
        algorithm. The ``radii`` must either be a float or monotonic sequence.
        The optional ``weights`` can be a float or a 1-dim.

        Returns an array of counts. The first element contains the count of all
        neighbours ``0 <= r <= r_1``, the following values contain the
        incremental counts ``r_i-1 <= r <= r_i``.
        """
        if not isinstance(other, self.__class__):
            raise TypeError("'other' must be of type 'AngularTree'")
        radii = coord.angle_to_chorddist(angles)
        return super().dualcount_range(other, radii)
