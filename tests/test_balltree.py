import numpy as np
import numpy.testing as npt
import pytest

from balltree import BallTree, default_leafsize

radius_testvalues = [0.01, 0.1, 1.0, 2.0]
rminmax_testvalues = [(0.1, 0.2), (1.0, 2.0), (2.0, 3.0)]


@pytest.fixture
def mock_data():
    # NOTE: if values are changed, mock_radius changes
    return np.array(
        [
            [5.1, 5.2, 5.3],
            [6.1, 6.2, 6.3],
            [7.1, 7.2, 7.3],
            [3.1, 3.2, 3.3],
            [2.1, 2.2, 2.3],
            [4.1, 4.2, 4.3],  # median
            [1.1, 1.2, 1.3],
        ]
    )


@pytest.fixture
def mock_data_median(mock_data):
    return np.median(mock_data, axis=0)


@pytest.fixture
def mock_weight(mock_data):
    return np.arange(len(mock_data))


@pytest.fixture
def mock_center(mock_data_median):
    return mock_data_median


@pytest.fixture
def mock_radius():
    return np.sqrt(27)


@pytest.fixture
def mock_tree(mock_data, mock_weight):
    return BallTree(mock_data, weight=mock_weight)


@pytest.fixture
def rand_data_weight():
    rng = np.random.default_rng(12345)
    size = 1000
    return (rng.uniform(-1.0, 1.0, size=(size, 3)), rng.normal(1.0, 0.02, size))


def data_to_view(data, weight=True):
    dtype = [("x", "f8"), ("y", "f8"), ("z", "f8"), ("weight", "f8")]
    data = np.atleast_2d(data)
    array = np.empty(len(data), dtype=dtype)
    array["x"] = data[:, 0]
    array["y"] = data[:, 1]
    array["z"] = data[:, 2]
    array["weight"] = weight if weight is not None else 1.0
    return array


def euclidean_distance(points, tpoint):
    diff = points - tpoint[np.newaxis, :]
    diffsq = diff**2
    return np.sqrt(diffsq.sum(axis=1))


def brute_force(data_weight, point_weight, radius):
    data, weight = data_weight
    point, pweight = point_weight
    dist = euclidean_distance(data, point)
    mask = dist <= radius
    return pweight * np.sum(weight[mask])


class TestBallTree:
    def test_init_no_xyz(mock_data):
        with pytest.raises(TypeError):
            BallTree()
        with pytest.raises(TypeError):
            BallTree(weight=[0.5])
        with pytest.raises(TypeError):
            BallTree(leafsize=default_leafsize)

    def test_init_empty_xyz(mock_data):
        with pytest.raises(ValueError, match="at least one element"):
            BallTree(np.empty(shape=(0, 3)))
        with pytest.raises(ValueError, match="shape"):
            BallTree([])

    def test_init_single(self, mock_data_median):
        tree = BallTree(mock_data_median)
        assert tree.leafsize == default_leafsize
        assert tree.num_points == 1
        assert tree.count_nodes() == 1
        npt.assert_array_equal(tree.data, data_to_view(mock_data_median))

    def test_init_single_list(self, mock_data_median):
        tree = BallTree(mock_data_median.tolist())
        npt.assert_array_equal(tree.data, data_to_view(mock_data_median))

    def test_init_single_wrong_shape(self, mock_data_median):
        with pytest.raises(ValueError, match="shape"):
            BallTree(mock_data_median[:-1])

    def test_init_single_weight(self, mock_data_median):
        weight = np.array([0.5])
        tree = BallTree(mock_data_median, weight)
        npt.assert_array_equal(tree.data, data_to_view(mock_data_median, weight))

    def test_init_single_weight_kwarg(self, mock_data_median):
        weight = np.array([0.5])
        tree = BallTree(mock_data_median, weight=weight)
        npt.assert_array_equal(tree.data, data_to_view(mock_data_median, weight))

    def test_init_single_weight_list(self, mock_data_median):
        weight = [0.5]
        tree = BallTree(mock_data_median, weight)
        npt.assert_array_equal(tree.data, data_to_view(mock_data_median, weight))

    def test_init_single_weight_scalar(self, mock_data_median):
        weight = 0.5
        tree = BallTree(mock_data_median, weight)
        npt.assert_array_equal(tree.data, data_to_view(mock_data_median, [weight]))

    def test_init_single_weight_wrong_shape(self, mock_data_median):
        with pytest.raises(ValueError, match="same length"):
            BallTree(mock_data_median, [0.5, 0.5])
        with pytest.raises(ValueError, match="same length"):
            BallTree(mock_data_median, [])

    def test_init_multi_leaf_only(self, mock_data):
        tree = BallTree(mock_data)
        assert tree.leafsize == default_leafsize
        assert tree.num_points == len(mock_data)
        assert tree.count_nodes() == 1
        npt.assert_array_equal(tree.data, data_to_view(mock_data))

    def test_init_multi_two_childs(self, mock_data, mock_data_median):
        leafsize = 4
        tree = BallTree(mock_data, leafsize=leafsize)
        assert tree.leafsize == leafsize
        assert tree.num_points == len(mock_data)
        assert tree.count_nodes() == 3
        # check the partitioning
        pivot_x = mock_data_median[0]
        idx_pivot = tree.num_points // 2
        npt.assert_array_equal(tree.data[idx_pivot], data_to_view(mock_data_median))
        assert np.all(tree.data[:idx_pivot]["x"] <= pivot_x)
        assert np.all(tree.data[idx_pivot + 1 :]["x"] >= pivot_x)
        # check the data elements where swapped correctly
        npt.assert_array_equal(np.sort(tree.data), np.sort(data_to_view(mock_data)))

    def test_init_multi_wrong_shape(self, mock_data):
        with pytest.raises(ValueError, match="shape"):
            BallTree(mock_data[:, :-1])

    def test_init_multi_list(self, mock_data):
        mock_list = mock_data.tolist()
        npt.assert_array_equal(BallTree(mock_list).data, BallTree(mock_data).data)

    def test_init_multi_weight(self, mock_data, mock_weight):
        tree = BallTree(mock_data, mock_weight)
        tree_data_sorted = np.sort(tree.data, order="weight")
        npt.assert_array_equal(tree_data_sorted, data_to_view(mock_data, mock_weight))

    def test_init_wrong_types(self):
        with pytest.raises(ValueError):
            BallTree(None)
        with pytest.raises(ValueError):
            BallTree(1.0)
        with pytest.raises(TypeError):
            BallTree(["asldfkj"])

    def test_data(self, mock_data):
        # checks if the lifetime of the underlying buffer is correclty managed
        npt.assert_array_equal(data_to_view(mock_data), BallTree(mock_data).data)

    def test_sum_weight(self, mock_tree, mock_weight):
        assert mock_tree.sum_weight == mock_weight.sum()

    def test_center(self, mock_tree, mock_center):
        npt.assert_array_equal(mock_tree.center, mock_center)

    def test_radius(self, mock_tree, mock_radius):
        assert mock_tree.radius == mock_radius

    def test_from_random(self):
        low = -1.0
        high = 1.0
        size = 100000
        tree = BallTree.from_random(low, high, size)
        data = np.array(tree.data.tolist())
        assert len(data) == size
        assert data.min() >= low
        assert data.max() <= high
        assert np.all(data[:, -1] == 1.0)

    def test_to_from_file(self, mock_data, mock_weight, tmp_path):
        fpath = str(tmp_path / "tree.dump")
        orig = BallTree(mock_data, mock_weight, leafsize=4)
        orig.to_file(fpath)

        restored = BallTree.from_file(fpath)
        assert orig.leafsize == restored.leafsize
        assert orig.num_points == restored.num_points
        assert orig.count_nodes() == restored.count_nodes()
        npt.assert_array_equal(orig.data, restored.data)

    @pytest.mark.parametrize("radius", radius_testvalues)
    def test_brute_radius_multi(self, radius, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        count = 0.0
        for p, w in zip(data, weight):
            count += brute_force((data, weight), (p, w), radius)
        npt.assert_almost_equal(tree.brute_radius(data, radius, weight), count)

    def test_count_radius_no_radius(self, rand_data_weight):
        data, _ = rand_data_weight
        tree = BallTree(data)
        with pytest.raises(TypeError):
            tree.count_radius(data[0])

    @pytest.mark.parametrize("radius", radius_testvalues)
    def test_count_radius_single(self, radius, rand_data_weight):
        data, _ = rand_data_weight
        weight = np.ones(len(data))
        tree = BallTree(data)

        p = data[0]
        w = 1.0
        count = brute_force((data, weight), (p, w), radius)
        npt.assert_almost_equal(tree.count_radius(p, radius), count)

    @pytest.mark.parametrize("radius", radius_testvalues)
    def test_count_radius_single_weight(self, radius, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        p = data[0]
        w = weight[0]
        count = brute_force((data, weight), (p, w), radius)
        npt.assert_almost_equal(tree.count_radius(p, radius, weight=w), count)

    @pytest.mark.parametrize("radius", radius_testvalues)
    def test_count_radius_multi(self, radius, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        count = 0.0
        for p, w in zip(data, weight):
            count += brute_force((data, weight), (p, w), radius)
        npt.assert_almost_equal(tree.count_radius(data, radius, weight), count)

    @pytest.mark.parametrize("radius", radius_testvalues)
    def test_dualcount_radius(self, radius, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        count = 0.0
        for p, w in zip(data, weight):
            count += brute_force((data, weight), (p, w), radius)
        npt.assert_almost_equal(tree.dualcount_radius(tree, radius), count)

    @pytest.mark.parametrize("rmin,rmax", rminmax_testvalues)
    def test_brute_range_multi(self, rmin, rmax, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        count_min = 0.0
        count_max = 0.0
        for p, w in zip(data, weight):
            count_min += brute_force((data, weight), (p, w), rmin)
            count_max += brute_force((data, weight), (p, w), rmax)
        result = [count_min, count_max - count_min]
        npt.assert_almost_equal(tree.brute_range(data, (rmin, rmax), weight), result)

    def test_count_range_single_radius(self, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        p = data[0]
        expect = [tree.count_radius(p, 0.2)]
        npt.assert_almost_equal(tree.count_range(p, 0.2), expect)
        npt.assert_almost_equal(tree.count_range(p, [0.2]), expect)

    def test_count_range_many_radii(self, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        radii = [0.1, 0.2, 0.3, 0.4]
        p = data[0]
        count = [tree.count_radius(p, r) for r in radii]
        expect = np.diff(count, prepend=0.0)
        npt.assert_almost_equal(tree.count_range(p, radii), expect)

    @pytest.mark.parametrize("rmin,rmax", rminmax_testvalues)
    def test_count_range_single(self, rmin, rmax, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        p = data[0]
        w = weight[0]
        count_min = brute_force((data, weight), (p, w), rmin)
        count_max = brute_force((data, weight), (p, w), rmax)
        result = [count_min, count_max - count_min]
        npt.assert_almost_equal(tree.count_range(p, (rmin, rmax), weight=w), result)

    @pytest.mark.parametrize("rmin,rmax", rminmax_testvalues)
    def test_dualcount_range(self, rmin, rmax, rand_data_weight):
        data, weight = rand_data_weight
        tree = BallTree(data, weight)

        count_min = 0.0
        count_max = 0.0
        for p, w in zip(data, weight):
            count_min += brute_force((data, weight), (p, w), rmin)
            count_max += brute_force((data, weight), (p, w), rmax)
        result = [count_min, count_max - count_min]
        npt.assert_almost_equal(tree.dualcount_range(tree, (rmin, rmax)), result)
