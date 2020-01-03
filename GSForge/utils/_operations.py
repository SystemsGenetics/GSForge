import numpy as np


def shuffle_along_axis(a, axis):
    """
    Shuffle along an axis. `Source <https://stackoverflow.com/a/55317373>`_.
    """
    idx = np.random.rand(*a.shape).argsort(axis=axis)
    return np.take_along_axis(a, idx, axis=axis)
