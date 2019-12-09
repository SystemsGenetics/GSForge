"""
`Analytics` are injective functions, *i.e.* those that returns only one set of results for each set of
arguments provided.

They usually operate by requesting X and Y data via the `GEMprospector.Interface` class in their `process()`
function, and then return an `Xarray.Dataset` object.
"""

import warnings

import param
import xarray as xr
import numpy as np

from sklearn.feature_selection import chi2, f_classif
from sklearn.model_selection import cross_val_score, train_test_split

from ..models import OperationInterface

__all__ = [
    "get_data",
    "rank_genes_by_model",
    "calculate_null_rank_distribution",
    "chi_squared_test",
    "f_classification_test",
    "basic_DESeq2_test"
]


class get_data(OperationInterface):
    """
    Gets the GEM matrix and an optional annotation column.
    """

    def process(self):
        return self.x_data, self.y_data


class gs_train_test_split(OperationInterface):
    """
    Performs an `sklearn.preprocessing.train_test_split()` call on the subset of data
    specified by the interface options (the same options passed to `get_data()`.
    """
    train_test_split_options = param.Parameter(default=dict())

    def process(self):
        # Get the subset of data selected by this operation.
        y_index = self.y_data[self.gem.sample_index_name].copy(deep=True)

        # Get the sample index and make the train and test indexes.
        train_idx, test_idx = train_test_split(y_index, **self.train_test_split_options)

        x_train = self.x_data.sel({self.gem.sample_index_name: train_idx})
        x_test = self.x_data.sel({self.gem.sample_index_name: test_idx})

        y_train = self.y_data.sel({self.gem.sample_index_name: train_idx})
        y_test = self.y_data.sel({self.gem.sample_index_name: test_idx})

        return x_train, x_test, y_train, y_test


class chi_squared_test(OperationInterface):
    """
    Compute chi-squared stats between each non-negative feature and class.
    See the `Scikit-learn documentation <https://scikit-learn.org/>`_
    """

    def process(self):
        x_data, y_data = self.x_data, self.y_data

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

    def process(self):
        x_data, y_data = self.x_data, self.y_data

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
    model = param.Parameter()
    n_iterations = param.Integer(default=1)

    def process(self):
        if self.n_iterations == 1:
            return self._rank_genes_by_model()
        else:
            rankings = [self._rank_genes_by_model() for _ in range(self.n_iterations)]
            ranking_ds = xr.concat(rankings, "feature_importance_iter")
            ranking_ds.name = "feature_importances"
            ranking_ds["feature_importance_iter"] = (["feature_importance_iter", ], np.arange(self.n_iterations))
            ranking_ds["feature_importance_mean"] = ([self.gem.gene_index_name],
                                                     ranking_ds.mean(dim="feature_importance_iter"))
            ranking_ds["feature_importance_std"] = ([self.gem.gene_index_name],
                                                    ranking_ds.std(dim="feature_importance_iter"))
            return ranking_ds.to_dataset()

    def _rank_genes_by_model(self):
        x_data = self.x_data
        y_data = self.y_data

        if isinstance(self.y_variables, list):
            y_data = y_data.to_dataframe().values

        model = self.model.fit(x_data, y_data)

        attrs = {'Ranking Model': str(model),
                 "x_variable": self.x_variable,
                 "y_variables": self.y_variables}

        data = xr.DataArray(data=model.feature_importances_,
                            dims=[self.gem.gene_index_name],
                            coords=[self.get_gene_index()],
                            attrs=attrs)
        return data


def _null_rank_distribution(a, b):
    return np.sum(a[:, None] <= b[None, :], axis=-1) / b.shape[0]


class calculate_null_rank_distribution(OperationInterface):
    """
    Probability for an irrelevant feature to be ranked above or at the same position as a given gene.
    """
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
        x_values = shadowed_subset[self.x_variable].values
        y_values = shadowed_subset[self.y_variables].values

        if isinstance(self.y_variables, list):
            y_values = y_values.to_dataframe().values

        model = self.model.fit(x_values, y_values)
        ranks = model.feature_importances_
        real, shadow = ranks.reshape((2, selected_genes.shape[0]))
        null_rank_dist = _null_rank_distribution(real, shadow)

        attrs = {'Ranking Model': str(model),
                 "x_variable": self.x_variable,
                 "y_variables": self.y_variables}

        return xr.DataArray(data=null_rank_dist, coords=[self.get_gene_index()],
                            dims=[self.gem.gene_index_name], attrs=attrs,
                            name="null_rank_distribution")


class calculate_family_wise_error_rates(OperationInterface):
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
        random_values = np.random.randn(*self.x_data.shape)
        shadowed_values = np.hstack([self.x_data, random_values])

        y_values = self.y_data

        if isinstance(self.y_variables, list):
            y_values = y_values.to_dataframe().values

        model = self.model.fit(shadowed_values, y_values)
        ranks = model.feature_importances_
        real, shadow = ranks.reshape((2, self.x_data.shape[1]))

        null_rank_dist = _null_rank_distribution(real, shadow)

        attrs = {'Ranking Model': str(model),
                 "x_variable": self.x_variable,
                 "y_variables": self.y_variables}

        return xr.DataArray(data=null_rank_dist, coords=[self.get_gene_index()],
                            dims=[self.gem.gene_index_name], attrs=attrs,
                            name="family_wise_error_rate")


class cv_score_model(OperationInterface):
    """
    Document me!
    """
    model = param.Parameter()
    cv = param.Parameter()

    def process(self):
        return cross_val_score(self.model, self.x_data, self.y_data, cv=self.cv)


# TODO: Cite original gist. Consider moving entirely to R for brevity.
#       Leave this implementation or a reference to it for scripting
#       purposes.
# TODO: Add support or a way to get more than one formula output.
class basic_DESeq2_test(OperationInterface):
    """
    Just a simple API to a standard DESeq2 run.

    Modified from this gist:
    https://gist.github.com/wckdouglas/3f8fb27a3d7a1eb24c598aa04f70fb25
    """

    formula = param.Parameter()

    def process(self):
        # Import R packages.
        try:
            import rpy2.robjects as robjects
            from rpy2.robjects import pandas2ri, Formula
            from rpy2.robjects.packages import importr

        except ImportError:
            # TODO: Expand upon packages required.
            warnings.warn("You must install 'rpy2' and other R packages for this"
                          "analytic to function.")

        # Setup parameters for DESeq2.
        self.set_param(count_mask="dropped")
        deseq = importr('DESeq2')
        pandas2ri.activate()

        # Ready the count matrix for R / DESeq2.
        x_data, y_data = self.x_data, self.y_data

        # Convert to a pandas dataframe and transpose.
        count_df = x_data.to_pandas().copy().transpose()

        # Ensure the counts are integers.
        count_df = count_df.astype(dtype='int32')
        r_count_matrix = pandas2ri.py2rpy(count_df)

        # Ready the design matrix for R / DESeq2.
        label_df = self.gem.data[self.y_variables].to_dataframe().copy()
        r_design_matrix = pandas2ri.py2rpy(label_df)

        r_design_formula = Formula(self.formula)

        dds = deseq.DESeqDataSetFromMatrix(countData=r_count_matrix,
                                           colData=r_design_matrix,
                                           design=r_design_formula)

        dds = deseq.DESeq(dds)
        deseq_result = deseq.results(dds)
        to_dataframe = robjects.r('function(x) data.frame(x)')
        deseq_result = to_dataframe(deseq_result)

        # Convert the output to an Xarray.Dataset object.
        deseq_result["Gene"] = x_data["Gene"].values
        deseq_result = deseq_result.set_index("Gene")

        data = deseq_result.to_xarray()

        attrs = {'model': "deseq2",
                 'formula': str(self.formula),
                 "x_variable": self.x_variable,
                 "y_variables": self.y_variables}

        data = data.assign_attrs(attrs)
        return data
