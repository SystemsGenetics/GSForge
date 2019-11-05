import numpy as np
import holoviews as hv
import matplotlib.pyplot as plt
from tqdm.autonotebook import tqdm
import seaborn as sns
import matplotlib.patches as mpatches

from ..models import OperationInterface


class SampleDistribution(OperationInterface):

    @staticmethod
    def np_sample_distributions(counts: np.ndarray, labels: np.ndarray = None):
        """
        Calculate and overlay kernel density estimates of the given count matrix on a per-sample basis.

        If this gives a ValueError you may need to install statsmodels for a more robust kernel estimate.

        :param numpy.ndarray counts:  The count matrix to be displayed.
        :param numpy.ndarray labels: If provided the output density estimates will be colored by their
            label membership.

        :returns: matplotlib.axes with overlayed kernel density estimates.
        """

        fig, ax = plt.subplots(1, figsize=(15, 8))

        if labels is not None:
            label_set = list(np.unique(labels))
            colors = {label: color for label, color in zip(label_set, sns.color_palette(n_colors=len(label_set)))}

        for index, sample_row in tqdm(enumerate(counts), total=counts.shape[0]):

            if labels is not None:
                label = labels[index]
                color = colors[label]
                sns.kdeplot(sample_row, ax=ax, legend=True, shade=False, gridsize=250, color=color)

            else:
                sns.kdeplot(sample_row, ax=ax, legend=False, shade=False, gridsize=250)
        if labels is not None:
            patches = [mpatches.Patch(color=color, label=label)
                       for label, color in colors.items()]
            plt.legend(handles=patches)

        return ax
