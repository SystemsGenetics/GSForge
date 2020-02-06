import param
import numpy as np
import holoviews as hv
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

from GSForge.models import OperationInterface


# TODO: Convert to a holoviews-based function.
class SampleWiseDistribution(OperationInterface):

    ax = param.Parameter()

    @staticmethod
    def np_sample_distributions(counts: np.ndarray, labels: np.ndarray = None, ax=None):
        """
        Calculate and overlay kernel density estimates of the given count matrix on a per-sample basis.

        If this gives a ValueError you may need to install `statsmodels` for a more robust kernel estimate.

        :param numpy.ndarray counts:
            The count matrix to be displayed.

        :param numpy.ndarray labels:
            If provided the output density estimates will be colored by their label membership.

        :returns:
            `matplotlib.axes` with per-sample kernel density estimates.
        """

        if ax is None:
            ax = plt.gca()

        if labels is not None:
            label_set = list(np.unique(labels))
            colors = {label: color for label, color in zip(label_set, sns.color_palette(n_colors=len(label_set)))}

        for index, sample_row in enumerate(counts):

            if labels is not None:
                label = labels[index]
                color = colors[label]
                sns.kdeplot(sample_row, ax=ax, legend=True, shade=False, gridsize=250, color=color)

            else:
                sns.kdeplot(sample_row, ax=ax, legend=False, shade=False, gridsize=250)
        if labels is not None:
            patches = [mpatches.Patch(color=color, label=label)
                       for label, color in colors.items()]
            ax.legend(handles=patches)

        return ax

    def process(self):
        """

        :return:
        """
        counts = self.x_count_data.values.copy()
        labels = self.y_annotation_data.values.copy() if self.y_annotation_data is not None else None
        ax = self.np_sample_distributions(counts, labels, ax=self.ax)
        ax.set_title(f"{self.active_count_variable} value distribution per-sample.");

        return ax
