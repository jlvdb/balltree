import numpy as np
import numpy.testing as npt
import pytest

from balltree import coordinates


def test_sgn():
    values = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    result = np.array([-1.0, -1.0, 1.0, 1.0, 1.0])
    npt.assert_array_equal(coordinates.sgn(values), result)


@pytest.mark.skip
def test_angular_to_cylinder():
    pass


@pytest.mark.skip
def test_cylinder_to_angular():
    pass


@pytest.mark.skip
def test_angular_to_euclidean():
    pass


@pytest.mark.skip
def test_euclidean_to_angular():
    pass


@pytest.mark.skip
def test_angle_to_chorddist():
    pass


@pytest.mark.skip
def test_chorddist_to_angle():
    pass
