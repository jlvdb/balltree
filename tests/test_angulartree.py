import numpy as np
import numpy.testing as npt
import pytest

from balltree import coordinates, default_leafsize
from balltree.angulartree import AngularTree


@pytest.fixture
def mock_data():
    return np.deg2rad(
        [
            [0.0, 0.0],
            [90.0, 0.0],
            [180.0, 0.0],
            [270.0, 0.0],
            [0.0, 90.0],
            [0.0, -90.0],
        ]
    )


@pytest.fixture
def mock_data_small():
    return np.deg2rad(
        [
            [0.0, -45.0],
            [0.0, 0.0],
            [0.0, 45.0],
        ]
    )


@pytest.fixture
def mock_tree(mock_data):
    return AngularTree(mock_data)


def data_to_view(data, weight=True):
    dtype = [("ra", "f8"), ("dec", "f8"), ("weight", "f8")]
    data = np.atleast_2d(data)
    array = np.empty(len(data), dtype=dtype)
    array["ra"] = data[:, 0]
    array["dec"] = data[:, 1]
    array["weight"] = weight if weight is not None else 1.0
    return array


class TestAngularTree:
    def test_init(self, mock_data):
        tree = AngularTree(mock_data)
        xyz_tree = tree._tree.data
        xyz_data = coordinates.angular_to_euclidean(mock_data)
        npt.assert_almost_equal(xyz_data[:, 0], xyz_tree["x"])
        npt.assert_almost_equal(xyz_data[:, 1], xyz_tree["y"])
        npt.assert_almost_equal(xyz_data[:, 2], xyz_tree["z"])

    def test_data(self, mock_tree, mock_data):
        npt.assert_array_equal(mock_tree.data, data_to_view(mock_data))

    def test_num_data(self, mock_tree, mock_data):
        assert mock_tree.num_data == len(mock_data)

    def test_leafsize(self, mock_tree):
        assert mock_tree.leafsize == default_leafsize

    def test_sum_weight(self, mock_tree):
        assert mock_tree.sum_weight == np.ones(mock_tree.num_data).sum()

    def test_center(self, mock_data_small):
        tree = AngularTree(mock_data_small)
        npt.assert_array_almost_equal(tree.center, np.median(mock_data_small, axis=0))

    def test_radius(self, mock_data_small):
        tree = AngularTree(mock_data_small)
        npt.assert_almost_equal(tree.radius, np.deg2rad(45))

    def test_from_random(self):
        limit = 1.0
        size = 10000
        tree = AngularTree.from_random(0, limit, -limit, limit, size)
        data = tree.data
        assert tree.num_data == size
        assert data["ra"].min() >= 0.0
        assert data["ra"].max() <= limit
        assert data["dec"].min() >= -limit
        assert data["dec"].max() <= limit

    def test_to_from_file(self, mock_data, tmp_path):
        fpath = str(tmp_path / "tree.dump")
        orig = AngularTree(mock_data, leafsize=4)
        orig.to_file(fpath)

        restored = AngularTree.from_file(fpath)
        assert orig.leafsize == restored.leafsize
        assert orig.num_data == restored.num_data
        assert orig.count_nodes() == restored.count_nodes()
        npt.assert_array_equal(orig.data, restored.data)

    def test_count_nodes(self, mock_data):
        assert AngularTree(mock_data, leafsize=4).count_nodes() == 3

    def test_brute_radius(self, mock_tree):
        point = (0.0, 0.0, 0.0)
        eps = 1e-9
        radius = np.pi / 2.0
        assert mock_tree.brute_radius(point, radius - eps) == 1
        assert mock_tree.brute_radius(point, radius + eps) == 5
        assert mock_tree.brute_radius(point, 2.0 * radius + eps) == 6

    def test_brute_range(self, mock_tree):
        point = (0.0, 0.0, 0.0)
        eps = 1e-9
        radius = np.pi / 2.0
        assert mock_tree.brute_range(point, radius - eps, 2.0 * radius + eps) == 5

    def test_count_radius(self, mock_tree):
        point = (0.0, 0.0, 0.0)
        eps = 1e-9
        radius = np.pi / 2.0
        assert mock_tree.count_radius(point, radius - eps) == 1
        assert mock_tree.count_radius(point, radius + eps) == 5
        assert mock_tree.count_radius(point, 2.0 * radius + eps) == 6

    def test_count_range(self, mock_tree):
        point = (0.0, 0.0, 0.0)
        eps = 1e-9
        radius = np.pi / 2.0
        assert mock_tree.count_range(point, radius - eps, 2.0 * radius + eps) == 5

    def test_dualcount_radius(self, mock_tree):
        eps = 1e-9
        radius = np.pi / 2.0
        N = mock_tree.num_data
        assert mock_tree.dualcount_radius(mock_tree, radius - eps) == 1 * N
        assert mock_tree.dualcount_radius(mock_tree, radius + eps) == 5 * N
        assert mock_tree.dualcount_radius(mock_tree, 2.0 * radius + eps) == 6 * N

    def test_dualcount_range(self, mock_tree):
        eps = 1e-9
        radius = np.pi / 2.0
        N = mock_tree.num_data
        assert (
            mock_tree.dualcount_range(mock_tree, radius - eps, 2.0 * radius + eps)
            == 5 * N
        )
