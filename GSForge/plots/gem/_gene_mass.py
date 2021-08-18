# import itertools
#
# import numpy as np
#
# import holoviews as hv
#
#
# def _create_paths(count_array: np.ndarray):
#     """
#     Returns an iterable containing coordinates for a single path, where the genes are integer-indexed
#     along the x-axis, and the mass of the gene (length * counts), is along the y-axis.
#
#     This function should not need to be called directly.
#
#     :param count_array:
#
#     :return:
#         An iterable containing coordinate tuples (x, y) for the gene mass path.
#     """
#     index_range = np.arange(count_array.shape[0])
#     starting_coords = list(zip(index_range, np.zeros(index_range.shape)))
#     ending_coords = list(zip(index_range, count_array))
#     return itertools.chain.from_iterable(zip(starting_coords, ending_coords, starting_coords))
#
#
# def plot_sample(sample: np.ndarray):
#     """
#     Plots the gene counts for the given sample array as a bar-plot like single-path element.
#
#     :param np.ndarray sample:
#         An array containing the count values.
#
#     :return:
#     """
#     paths = _create_paths(sample)
#     plot = hv.Path(paths)
#     return plot.opts(width=400, height=400, padding=0.05)
#
#
# def plot_sample_mass(sample: np.ndarray, lengths: np.ndarray):
#     """
#     Plots gene masses (counts * lengths) for the given sample array.
#
#     :param np.ndarray sample:
#         An array containing the count values.
#
#     :param np.ndarray lengths:
#         An array containing the length values.
#
#     :return:
#     """
#
#     masses = sample * lengths
#     paths = _create_paths(masses)
#     plot = hv.Path(paths)
#     return plot.opts(width=400, height=400, padding=0.05)
#
#
# def plot_count_sum_mass(counts: np.ndarray, lengths: np.ndarray):
#     """
#     Plots the sum of the gene masses (along the sample axis) for the given count array.
#
#     :param counts:
#     :param lengths:
#     :return:
#     """
#     sample_counts_sum = counts.sum(axis=0)
#     mass_sums = sample_counts_sum * lengths
#     paths = _create_paths(mass_sums)
#     plot = hv.Path(paths)
#     return plot.opts(width=400, height=400, padding=0.05)
