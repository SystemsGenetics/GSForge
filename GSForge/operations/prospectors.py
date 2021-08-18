"""
``Prospector`` operations return either boolean support arrays or arrays of selected genes.
Prospector operations differ from analytics, in that they are not required to return a 'result' for every gene,
or return the same result each call.
"""

import enum
import json
from textwrap import dedent

import numpy as np
import param
import xarray as xr
from boruta import boruta_py
from sklearn.feature_selection import chi2, f_classif

from ..models import Interface
from ..utils import kwargs_overlap

__all__ = [
    "parse_boruta_model",
    "BorutaProspector",
    "ChiSquaredTest",
    "FClassificationTest"
]


class ChiSquaredTest(Interface, param.ParameterizedFunction):
    """
    Compute chi-squared stats between each non-negative feature and class.
    See the `Scikit-learn documentation <https://scikit-learn.org/>`_
    """

    def __call__(self, *args, **params):
        super().__init__(*args, **params)
        x_data, y_data = self.x_count_data, self.y_annotation_data

        chi2_scores, chi2_pvals = chi2(np.nan_to_num(x_data), y_data)
        attrs = {"Method": "Chi-Squared",
                 "count_variable": self.count_variable,
                 "annotation_variables": self.annotation_variables}
        data = xr.Dataset({"chi2_scores": (["Gene"], chi2_scores),
                           "chi2_pvals": (["Gene"], chi2_pvals)},
                          coords={"Gene": x_data[self.gem.gene_index_name]},
                          attrs=attrs)
        return data


class FClassificationTest(Interface, param.ParameterizedFunction):
    """
    Compute the ANOVA F-value for the provided sample.
    See the `Scikit-learn documentation <https://scikit-learn.org/>`_
    """

    def __call__(self, *args, **params):
        super().__init__(*args, **params)
        x_data, y_data = self.x_count_data, self.y_annotation_data

        f_scores, f_pvals = f_classif(np.nan_to_num(x_data), y_data)
        attrs = {"Method": "ANOVA F-value",
                 "count_variable": self.count_variable,
                 "annotation_variables": self.annotation_variables}
        data = xr.Dataset({"f_scores": (["Gene"], f_scores),
                           "f_pvals": (["Gene"], f_pvals)},
                          coords={"Gene": x_data[self.gem.gene_index_name]},
                          attrs=attrs)
        return data


class _model_parsers(enum.Enum):
    boruta = {
        "support": lambda d, m: ([d], np.array(getattr(m, "support_"), dtype=bool)),
        "support_weak": lambda d, m: ([d], np.array(getattr(m, "support_weak_"), dtype=bool)),
        "ranking": lambda d, m: ([d], np.array(getattr(m, "ranking_"))),
    }


def _parse_model(model, model_type, dim) -> dict:
    """Parse a model into a dictionary appropriate for constructing an ``xarray.Dataset``.

    :param model: The model object to be parsed.

    :param model_type: The key to the Enum ``_model_parsers``.

    :param dim:

    :return:
    """
    parsing_dict = _model_parsers[model_type].value
    return {key: func(dim, model) for key, func in parsing_dict.items()}


def parse_boruta_model(boruta_model, gene_coords, attrs=None, dim="Gene") -> xr.Dataset:
    """Convert a boruta model into an ``xarray.Dataset`` object.

    :param boruta_model: A boruta_py model.

    :param attrs: A dictionary to be assigned to the output dataset attrs.

    :param gene_coords: An array (index) of the genes passed to the boruta_model.

    :param dim: The name of the coordinate dimension.

    :return: An ``xarray.Dataset`` object.
    """
    model_data = _parse_model(boruta_model, "boruta", dim=dim)
    return xr.Dataset(model_data, coords={dim: gene_coords}, attrs=attrs)


# TODO: Consider a direct to GeneSet object option.
class BorutaProspector(Interface, param.ParameterizedFunction):
    """Runs a single instance of BorutaPy feature selection.

    This is just a simple wrapper for a boruta model that produces an
    ``xarray.Dataset`` object suitable for use in the creation of a
    ``GSForge.GeneSet`` object."""

    estimator = param.Parameter(doc=dedent("""\
    A supervised learning estimator, with a 'fit' method that returns the
    ``feature_importances_`` attribute. Important features must correspond to
    high absolute values in the `feature_importances_`."""))

    n_estimators = param.Parameter(default=1000, doc=dedent("""\
    If int sets the number of estimators in the chosen ensemble method.
    If 'auto' this is determined automatically based on the size of the
    dataset. The other parameters of the used estimators need to be set
    with initialisation."""))

    perc = param.Integer(default=100, doc=dedent("""\
    Instead of the max we use the percentile defined by the user, to pick
    our threshold for comparison between shadow and real features. The max
    tend to be too stringent. This provides a finer control over this. The
    lower perc is the more false positives will be picked as relevant but
    also the less relevant features will be left out. The usual trade-off.
    The default is essentially the vanilla Boruta corresponding to the max.
    """))

    alpha = param.Number(default=0.05, doc=dedent("""\
    Level at which the corrected p-values will get rejected in both
    correction steps."""))

    two_step = param.Boolean(default=True, doc=dedent("""\
    If you want to use the original implementation of Boruta with Bonferroni
    correction only set this to False."""))

    max_iter = param.Integer(default=100, doc=dedent("""\
    The number of maximum iterations to perform."""))

    random_state = param.Parameter(default=None, doc=dedent("""\
    If int, random_state is the seed used by the random number generator;
    If RandomState instance, random_state is the random number generator;
    If None, the random number generator is the RandomState instance used
    by ``np.random``."""))

    verbose = param.Integer(default=0, doc=dedent("""\
    Controls verbosity of output:
    - 0: no output
    - 1: displays iteration number
    - 2: which features have been selected already"""))

    def __call__(self, *args, **params):
        super().__init__(*args, **params)

        if len(self.annotation_variables) > 1:
            raise ValueError(f"This operation only accepts a single entry for `annotation_variables`.")
        x_data, y_data = self.x_count_data, self.y_annotation_data
        boruta_kwargs = kwargs_overlap(self, boruta_py.BorutaPy)
        boruta_model = boruta_py.BorutaPy(**boruta_kwargs)
        boruta_model.fit(x_data, y_data)

        # Convert attribute dictionaries to .json format.
        # We can't store nested dictionaries in the attrs.
        boruta_json_attrs = json.dumps(boruta_kwargs, default=lambda o: '<not serializable>')
        ranking_model_json_attrs = json.dumps(self.estimator.get_params(), default=lambda o: '<not serializable>')

        attrs = {'boruta_model': boruta_json_attrs,
                 'ranking_model': ranking_model_json_attrs,
                 "selected_count_variable": self.count_variable,
                 "selected_annotation_variables": self.annotation_variables}

        return parse_boruta_model(boruta_model,
                                  attrs=attrs,
                                  gene_coords=self.get_gene_index(),
                                  dim=self.gem.gene_index_name)
