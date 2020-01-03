"""
``Analytics`` are intended to more closely rank or compare a GEM subset, rather than the entire GEM.
These functions are intended for analyzing and comparing subsets generated by the functions found in ``prospectors``.


Methods and notation from [method_compare]_ used.

: :math:`LS` :
    Learning Sample, :math:`n` instances of input-output values.
: :math:`n` :
    Number of input-output value pairs in :math:`LS`.
: :math:`m` :
    Number of input variables (features or genes) in :math:`LS`.
: :math `X_i` :
    Input array of :math:`LS`. Ranges from :math:`i=1, ..., m`.
: :math:`LS` :
    An algorithm that outputs some relevance score, :math:`s_i`, for each input variable :math `X_i`.


.. [method_compare] `A comparison of per sample global scaling and per gene normalization methods for differential expression analysis of RNA-seq data <https://doi.org/10.1371/journal.pone.0176185>`_

"""

# import warnings

import param
import xarray as xr
import numpy as np
from sklearn.base import clone

from ..models import OperationInterface
from ..utils._operations import shuffle_along_axis, null_rank_distribution

__all__ = [
    "rank_genes_by_model",
    "nFDR",
    "mProbes",
]


class rank_genes_by_model(OperationInterface):
    """
    Given some machine learning model, this operation runs n_iterations
    and returns a summary dataset of the ranking results.
    """
    # TODO: Note that this uses the OperationInterface.

    model = param.Parameter()
    n_iterations = param.Integer(default=1)

    def process(self):
        if self.n_iterations == 1:
            return self._rank_genes_by_model()
        else:
            rankings = [self._rank_genes_by_model() for _ in range(self.n_iterations)]
            ranking_ds = xr.concat(rankings, "feature_importance_iter")
            ranking_ds["feature_importance_iter"] = (["feature_importance_iter", ], np.arange(self.n_iterations))
            ranking_ds["feature_importance_mean"] = ([self.gem.gene_index_name],
                                                     ranking_ds.mean(dim="feature_importance_iter"))
            ranking_ds["feature_importance_std"] = ([self.gem.gene_index_name],
                                                    ranking_ds.std(dim="feature_importance_iter"))
            return ranking_ds.to_dataset()

    def _rank_genes_by_model(self):
        x_data = self.x_count_data
        y_data = self.y_annotation_data

        if isinstance(self.annotation_variables, list):
            y_data = y_data.to_dataframe().values

        model = self.model.fit(x_data, y_data)

        attrs = {'Ranking Model': str(model),
                 "count_variable": self.count_variable,
                 "annotation_variables": self.annotation_variables}

        data = xr.DataArray(data=model.feature_importances_,
                            dims=[self.gem.gene_index_name],
                            coords=[self.get_gene_index()],
                            name="feature_importances",
                            attrs=attrs)
        return data


class nFDR(OperationInterface):
    """
    nFDR (False Discovery Rate) [method_compare]_.

    nFDR trains two models and compares their ``feature_importances_`` attributes to estimate
    the false discovery rate.

    The FDR estimated is the percent of instances a shuffled output feature has a higher feature
    importance score than the same non-shuffled feature score.

    This is repeated up to ``n_iterations``.
    """

    model = param.Parameter()
    n_iterations = param.Integer(default=1)

    @staticmethod
    def nFDR(x_data, y_data, model):
        y_shadowed = np.random.permutation(y_data)

        real_model = clone(model).fit(x_data, y_data)
        shadow_model = clone(model).fit(x_data, y_shadowed)

        real_scores = real_model.feature_importances_
        shadow_scores = shadow_model.feature_importances_

        null_rank_dist = null_rank_distribution(real_scores, shadow_scores)
        return null_rank_dist

    def _nFDR_to_xarray(self, nFRD_values):
        attrs = {'Ranking Model': str(self.model),
                 "count_variable": self.count_variable,
                 "annotation_variables": self.annotation_variables}
        return xr.DataArray(
            data=nFRD_values,
            coords=[self.get_gene_index()],
            dims=[self.gem.gene_index_name],
            attrs=attrs,
            name="nFDR")

    def process(self):

        x_data = self.x_count_data.values
        y_data = self.y_annotation_data.values

        if self.n_iterations == 1:
            nrd_values = self.nFDR(x_data, y_data, self.model)
            return self._nFDR_to_xarray(nrd_values)

        else:
            fdrs = [self.nFDR(x_data, y_data, self.model) for i in range(self.n_iterations)]
            fdrs = [self._nFDR_to_xarray(values) for values in fdrs]
            fdr_ds = xr.concat(fdrs, "nFDR_iter")
            fdr_ds.name = "nFDR"
            fdr_ds["nFDR_iter"] = (["nFDR_iter", ], np.arange(self.n_iterations))
            fdr_ds["nFDR_mean"] = ([self.gem.gene_index_name], fdr_ds.mean(dim="nFDR_iter"))
            fdr_ds["nFDR_std"] = ([self.gem.gene_index_name], fdr_ds.std(dim="nFDR_iter"))
            return fdr_ds.to_dataset()


class mProbes(OperationInterface):
    """
    mProbes [method_compare]_ works by randomly permuting the feature values in the supplied data.
    e.g. count values are shuffled within each samples feature (gene) array.

    It then ranks the real and shadowed features (for ``n_iterations``) with the supplied ``model``
    via a call to ``model.fit()``. It then examines ``model.feature_importances_`` for the feature
    importance values, and then calculates the null rank distribution.

    This is repeated upto ``n_iterations``.
    """

    model = param.Parameter()
    n_iterations = param.Integer(default=1)

    @staticmethod
    def mProbes(x_data, y_data, model):
        shadowed_array = shuffle_along_axis(x_data, 1)
        x_shadowed = np.hstack((x_data, shadowed_array))
        model = model.fit(x_shadowed, y_data)
        ranks = model.feature_importances_
        real, shadow = ranks.reshape((2, x_data.shape[1]))
        null_rank_dist = null_rank_distribution(real, shadow)
        return null_rank_dist

    def _mprobes_fdr_to_xarray(self, values):
        attrs = {'Ranking Model': str(self.model),
                 "count_variable": self.count_variable,
                 "annotation_variables": self.annotation_variables}
        return xr.DataArray(
            data=values,
            coords=[self.get_gene_index()],
            dims=[self.gem.gene_index_name],
            attrs=attrs,
            name="nFDR")

    def process(self):

        x_data = self.x_count_data.values
        y_data = self.y_annotation_data.values

        if self.n_iterations == 1:
            return self._mprobes_fdr_to_xarray(self.mProbes(x_data, y_data, self.model))
        else:
            fdrs = [self.mProbes(x_data, y_data, self.model) for i in range(self.n_iterations)]
            fdrs = [self._mprobes_fdr_to_xarray(values) for values in fdrs]
            ranking_ds = xr.concat(fdrs, "mProbes_iter")
            ranking_ds.name = "mProbes NRD"
            ranking_ds["mProbes_iter"] = (["mProbes_iter", ], np.arange(self.n_iterations))
            ranking_ds["mProbes_mean"] = ([self.gem.gene_index_name], ranking_ds.mean(dim="mProbes_iter"))
            ranking_ds["mProbes_std"] = ([self.gem.gene_index_name], ranking_ds.std(dim="mProbes_iter"))
            return ranking_ds.to_dataset()
