import numpy as np
import numpy.testing as npt
import pytest

from balltree import coordinates


@pytest.fixture
def coords_angular_ra():
    return np.deg2rad(
        [
            [-90.0, 0.0],
            [0.0, 0.0],
            [90.0, 0.0],
            [180.0, 0.0],
            [270.0, 0.0],
            [360.0, 0.0],
            [450.0, 0.0],
        ]
    )


@pytest.fixture
def coords_cylinder_ra(coords_angular_ra):
    return coords_angular_ra


@pytest.fixture
def coords_euclidean_ra():
    return np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )


@pytest.fixture
def coords_angular_dec():
    return np.deg2rad(
        [
            [-90.0, 90.0],
            [0.0, 90.0],
            [90.0, 90.0],
            [-90.0, -90.0],
            [0.0, -90.0],
            [90.0, -90.0],
        ]
    )


@pytest.fixture
def coords_cylinder_dec(coords_angular_dec):
    values = coords_angular_dec.copy()
    values[:, 1] = coordinates.sgn(coords_angular_dec[:, 1])
    return values


@pytest.fixture
def coords_euclidean_dec():
    return np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, -1.0],
        ]
    )


@pytest.fixture
def chord_length():
    return np.array([0.0, np.sqrt(2.0), 2.0])


@pytest.fixture
def angle():
    return np.array([0.0, np.pi / 2.0, np.pi])


def test_sgn():
    values = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    result = np.array([-1.0, -1.0, 1.0, 1.0, 1.0])
    npt.assert_array_equal(coordinates.sgn(values), result)


def test_angular_to_cylinder_ra(coords_angular_ra, coords_cylinder_ra):
    result = coordinates.angular_to_cylinder(coords_angular_ra)
    expected = coords_cylinder_ra
    npt.assert_array_almost_equal(result, expected)


def test_angular_to_cylinder_dec(coords_angular_dec, coords_cylinder_dec):
    result = coordinates.angular_to_cylinder(coords_angular_dec)
    expected = coords_cylinder_dec
    npt.assert_array_almost_equal(result, expected)


def test_cylinder_to_angular_ra(coords_cylinder_ra, coords_angular_ra):
    result = coordinates.cylinder_to_angular(coords_cylinder_ra)
    expected = coords_angular_ra
    npt.assert_array_almost_equal(result, expected)


def test_cylinder_to_angular_dec(coords_cylinder_dec, coords_angular_dec):
    result = coordinates.cylinder_to_angular(coords_cylinder_dec)
    expected = coords_angular_dec
    npt.assert_array_almost_equal(result, expected)


def test_angular_to_euclidean_ra(coords_angular_ra, coords_euclidean_ra):
    result = coordinates.angular_to_euclidean(coords_angular_ra)
    expected = coords_euclidean_ra
    npt.assert_array_almost_equal(result, expected)


def test_angular_to_euclidean_dec(coords_angular_dec, coords_euclidean_dec):
    result = coordinates.angular_to_euclidean(coords_angular_dec)
    expected = coords_euclidean_dec
    npt.assert_array_almost_equal(result, expected)


def test_euclidean_to_angular_ra(coords_euclidean_ra, coords_angular_ra):
    result = coordinates.euclidean_to_angular(coords_euclidean_ra)
    # some points are outside of the range of [0, 2pi)
    expected = coords_angular_ra % (2.0 * np.pi)
    npt.assert_array_almost_equal(result, expected)


def test_euclidean_to_angular_dec(coords_euclidean_dec, coords_angular_dec):
    result = coordinates.euclidean_to_angular(coords_euclidean_dec)
    # all tested points are coordinate singularities, disregard RA component
    expected = coords_angular_dec.copy()
    expected[:, 0] = 0.0
    npt.assert_array_almost_equal(result, expected)


def test_angle_to_chorddist(angle, chord_length):
    result = coordinates.angle_to_chorddist(angle)
    expected = chord_length
    npt.assert_array_almost_equal(result, expected)


def test_chorddist_to_angle(chord_length, angle):
    result = coordinates.chorddist_to_angle(chord_length)
    expected = angle
    npt.assert_array_almost_equal(result, expected)
