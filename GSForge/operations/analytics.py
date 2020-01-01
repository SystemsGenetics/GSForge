"""
`Analytics` are  functions  that returns only one set of results for each set of arguments provided.
"""

# import warnings

import param
import xarray as xr
import numpy as np

from sklearn.feature_selection import chi2, f_classif
from sklearn.model_selection import cross_val_score, train_test_split

from ..models import OperationInterface

__all__ = [
    "get_data",
    "train_test_split_wrapper",
    "rank_genes_by_model",
    "calculate_null_rank_distribution",
    "chi_squared_test",
    "f_classification_test",
]


# TODO: Add linter ignore for class capitaliation.
class get_data(OperationInterface):
    """
    Gets the GEM matrix and an optional annotation column.
    """
    # TODO: Expand comment, describe how the sample and gene indexes are built.

    def process(self):
        return self.x_count_data, self.y_annotation_data


class train_test_split_wrapper(OperationInterface):
    """
    Performs an `sklearn.preprocessing.train_test_split()` call on the subset of data
    specified by the interface options (the same options passed to `get_data()`.

    :returns: x_train, x_test, y_train, y_test
    """
    # TODO: Add links and reference to the sklearn function and docs.
    train_test_split_options = param.Parameter(default=dict())

    def process(self):
        # Get the subset of data selected by this operation.
        y_index = self.get_sample_index()

        # Get the sample index and make the train and test indexes.
        train_idx, test_idx = train_test_split(y_index, **self.train_test_split_options)

        x_train = self.x_count_data.sel({self.gem.sample_index_name: train_idx})
        x_test = self.x_count_data.sel({self.gem.sample_index_name: test_idx})

        y_train = self.y_annotation_data.sel({self.gem.sample_index_name: train_idx})
        y_test = self.y_annotation_data.sel({self.gem.sample_index_name: test_idx})

        return x_train, x_test, y_train, y_test


class chi_squared_test(OperationInterface):
    """
    Compute chi-squared stats between each non-negative feature and class.
    See the `Scikit-learn documentation <https://scikit-learn.org/>`_
    """
    # TODO: Note that this uses the OperationInterface.

    def process(self):
        x_data, y_data = self.x_count_data, self.y_annotation_data

        chi2_scores, chi2_pvals = chi2(np.nan_to_num(x_data), y_data)
        attrs = {"Method": "Chi-Squared",
                 "x_variable": self.x_variable,
                 "y_variables": self.y_variables}
        data = xr.Dataset({"chi2_scores": (["Gene"], chi2_scores),
                           "chi2_pvals": (["Gene"], chi2_pvals)},
                          coords={"Gene": x_data[self.gem.gene_index_name]},
                          attrs=attrs)
        return data


class f_classification_test(OperationInterface):
    """
    Compute the ANOVA F-value for the provided sample.
    See the `Scikit-learn documentation <https://scikit-learn.org/>`_
    """
    # TODO: Note that this uses the OperationInterface.

    def process(self):
        x_data, y_data = self.x_count_data, self.y_annotation_data

        f_scores, f_pvals = f_classif(np.nan_to_num(x_data), y_data)
        attrs = {"Method": "ANOVA F-value",
                 "x_variable": self.x_variable,
                 "y_variables": self.y_variables}
        data = xr.Dataset({"f_scores": (["Gene"], f_scores),
                           "f_pvals": (["Gene"], f_pvals)},
                          coords={"Gene": x_data[self.gem.gene_index_name]},
                          attrs=attrs)
        return data


class rank_genes_by_model(OperationInterface):
    """Given some machine learning model, this operation runs n_iterations
    and returns a summary dataset of the ranking results."""
    # TODO: Note that this uses the OperationInterface.

    model = param.Parameter()
    n_iterations = param.Integer(default=1)

    def process(self):
        if self.n_iterations == 1:
            return self._rank_genes_by_model()
        else:
            rankings = [self._rank_genes_by_model() for _ in range(self.n_iterations)]
            ranking_ds = xr.concat(rankings, "feature_importance_iter")
            # ranking_ds.name = "feature_importances"
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


def _null_rank_distribution(a, b):
    # TODO: Cite and add some description.
    return np.sum(a[:, None] <= b[None, :], axis=-1) / b.shape[0]


class calculate_null_rank_distribution(OperationInterface):
    """
    Probability for an irrelevant feature to be ranked above or at the same position as a given gene.
    """
    # TODO: Cite and update docstring.
    model = param.Parameter()
    n_iterations = param.Integer(default=1)

    def process(self):
        if self.n_iterations == 1:
            return self._calculate_null_rank_distribution()
        else:
            rankings = [self._calculate_null_rank_distribution() for i in range(self.n_iterations)]
            ranking_ds = xr.concat(rankings, "nrd_iter")
            ranking_ds.name = "null rank distribution"
            ranking_ds["nrd_iter"] = (["nrd_iter", ], np.arange(self.n_iterations))
            ranking_ds["nrd_mean"] = ([self.gem.gene_index_name], ranking_ds.mean(dim="nrd_iter"))
            ranking_ds["nrd_std"] = ([self.gem.gene_index_name], ranking_ds.std(dim="nrd_iter"))
            return ranking_ds.to_dataset()

    def _calculate_null_rank_distribution(self):
        # Get the genes that make up this lineament.
        selected_genes = self.get_gene_index()

        # Get the genes that are not part of that subset.
        # print(selected_genes)
        unselected_genes = self.gem.gene_index.drop(
            selected_genes, dim=self.gem.gene_index_name)[self.gem.gene_index_name].values.copy()

        # Randomly pick genes from theses (supposedly) irrelevant features.
        irrelevant_genes = np.random.choice(unselected_genes, selected_genes.shape[0])
        shadowed_gene_index = np.append(selected_genes, irrelevant_genes)

        # Construct the shadowed subset.
        shadowed_subset = self.gem.data.sel({self.gem.gene_index_name: shadowed_gene_index,
                                             self.gem.sample_index_name: self.get_sample_index()})
        x_values = shadowed_subset[self.count_variable].values
        y_values = shadowed_subset[self.annotation_variables].values

        if isinstance(self.annotation_variables, list):
            y_values = y_values.to_dataframe().values

        model = self.model.fit(x_values, y_values)
        ranks = model.feature_importances_
        real, shadow = ranks.reshape((2, selected_genes.shape[0]))
        null_rank_dist = _null_rank_distribution(real, shadow)

        attrs = {'Ranking Model': str(model),
                 "count_variable": self.count_variable,
                 "annotation_variables": self.annotation_variables}

        return xr.DataArray(data=null_rank_dist, coords=[self.get_gene_index()],
                            dims=[self.gem.gene_index_name], attrs=attrs,
                            name="null_rank_distribution")


class calculate_family_wise_error_rates(OperationInterface):
    # TODO: Cite and update docstring.

    """
    Document me!
    """
    model = param.Parameter()
    n_iterations = param.Integer(default=1)

    def process(self):
        if self.n_iterations == 1:
            return self._calculate_family_wise_error_rates()
        else:
            rankings = [self._calculate_family_wise_error_rates() for i in range(self.n_iterations)]
            ranking_ds = xr.concat(rankings, "fwe_iter")
            ranking_ds.name = "Family-wise error rate"
            ranking_ds["fwe_iter"] = (["fwe_iter", ], np.arange(self.n_iterations))
            ranking_ds["fwe_mean"] = ([self.gem.gene_index_name], ranking_ds.mean(dim="fwe_iter"))
            ranking_ds["fwe_std"] = ([self.gem.gene_index_name], ranking_ds.std(dim="fwe_iter"))
            return ranking_ds.to_dataset()

    def _calculate_family_wise_error_rates(self):
        random_values = np.random.randn(*self.x_count_data.shape)
        shadowed_values = np.hstack([self.x_count_data, random_values])

        y_values = self.y_annotation_data

        if isinstance(self.annotation_variables, list):
            y_values = y_values.to_dataframe().values

        model = self.model.fit(shadowed_values, y_values)
        ranks = model.feature_importances_
        real, shadow = ranks.reshape((2, self.x_count_data.shape[1]))

        null_rank_dist = _null_rank_distribution(real, shadow)

        attrs = {'Ranking Model': str(model),
                 "count_variable": self.count_variable,
                 "annotation_variables": self.annotation_variables}

        return xr.DataArray(data=null_rank_dist, coords=[self.get_gene_index()],
                            dims=[self.gem.gene_index_name], attrs=attrs,
                            name="family_wise_error_rate")


# class cv_score_model(OperationInterface):
#     """
#     Document me!
#     """
#     # TODO: Cite and update docstring.
#
#     model = param.Parameter()
#     cv = param.Parameter()
#
#     def process(self):
#         return cross_val_score(self.model, self.x_count_data, self.y_annotation_data, cv=self.cv)
