import numpy as np


def shuffle_along_axis(array, axis):
    """
    Shuffle along an axis.

    :param array:
        A ``numpy.ndarray`` object.

    :param axis:
        The axis to be shuffled along.

    `Source <https://stackoverflow.com/a/55317373>`_.
    """
    idx = np.random.rand(*array.shape).argsort(axis=axis)
    return np.take_along_axis(array, idx, axis=axis)


def null_rank_distribution(real, shadow):
    """
    Returns the percent for which shadow values are ranked higher than a given real value.

    :param real: The 'real' scores.

    :param shadow: Scores from fake, or shadow, features.

    :return: An array of the same shape as ``real``, with the percent of times a shadow
      score had a higher value.
    """
    return np.sum(real[:, None] <= shadow[None, :], axis=-1) / shadow.shape[0]


