import numpy as np
import numpy.testing as npt
import pytest

from balltree import default_leafsize
from balltree.angulartree import AngularTree

PI_2 = np.pi / 2.0


@pytest.fixture
def mock_data():
    # NOTE: if values are changed, test_angles, test_count, test_count_range change
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


@pytest.fixture
def test_angles():
    eps = 1e-9
    return np.array([PI_2 - eps, PI_2 + eps, np.pi])


@pytest.fixture
def test_count():
    return np.array([1.0, 5.0, 6.0])


@pytest.fixture
def test_count_range():
    return np.diff([1, 5, 6], prepend=0.0)  # prepend is for count with itself


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
        npt.assert_array_equal(tree.data, data_to_view(mock_data))

    def test_init_wrong_shape(self, mock_data):
        with pytest.raises(ValueError, match="dimensions"):
            AngularTree(mock_data[:, :-1])  # covers cases of method calls

    def test_num_points(self, mock_tree, mock_data):
        assert mock_tree.num_points == len(mock_data)

    def test_leafsize(self, mock_tree):
        assert mock_tree.leafsize == default_leafsize

    def test_sum_weight(self, mock_tree):
        assert mock_tree.sum_weight == np.ones(mock_tree.num_points).sum()

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
        assert tree.num_points == size
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
        assert orig.num_points == restored.num_points
        assert orig.count_nodes() == restored.count_nodes()
        npt.assert_array_equal(orig.data, restored.data)

    def test_count_nodes(self, mock_data):
        assert AngularTree(mock_data, leafsize=4).count_nodes() == 3

    def test_brute_radius(self, mock_tree, test_angles, test_count):
        point = (0.0, 0.0)
        for angle, count in zip(test_angles, test_count):
            assert mock_tree.brute_radius(point, angle) == count

    def test_count_radius(self, mock_tree, test_angles, test_count):
        point = (0.0, 0.0)
        for angle, count in zip(test_angles, test_count):
            assert mock_tree.count_radius(point, angle) == count

    def test_dualcount_radius(self, mock_tree, test_angles, test_count):
        N = mock_tree.num_points
        for angle, count in zip(test_angles, test_count):
            assert mock_tree.dualcount_radius(mock_tree, angle) == N * count

    def test_brute_range(self, mock_tree, test_angles, test_count_range):
        point = (0.0, 0.0)
        npt.assert_array_almost_equal(
            mock_tree.brute_range(point, test_angles), test_count_range
        )

    def test_count_range(self, mock_tree, test_angles, test_count_range):
        point = (0.0, 0.0)
        npt.assert_array_almost_equal(
            mock_tree.count_range(point, test_angles), test_count_range
        )

    def test_dualcount_range(self, mock_tree, test_angles, test_count_range):
        N = mock_tree.num_points
        npt.assert_array_almost_equal(
            mock_tree.dualcount_range(mock_tree, test_angles), test_count_range * N
        )
