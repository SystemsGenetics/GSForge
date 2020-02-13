"""
Test configuration for GSForge.
"""
import numpy as np
import pandas as pd
import xarray as xr

from sklearn.datasets import make_multilabel_classification

import GSForge as gsf
import pytest

# Control the size of the generated sample data.
N_SAMPLES = 100
N_GENES = 1000
N_CLASSES = 4
N_LABELS = 2
GENE_SETS = 4


def _gene_names():
    """Placeholder gene names."""
    return [f"gene_{x}" for x in range(N_GENES)]


def _sample_names():
    """Placeholder sample names."""
    return [f"sample_{x}" for x in range(N_SAMPLES)]


@pytest.fixture(scope="session")
def random_annotated_gem():
    # TODO: Replace this with a call to the utility function to be created.
    #       It can generate its own gene and sample names.
    data, labels = make_multilabel_classification(
        n_samples=N_SAMPLES,
        n_features=N_GENES,
        n_classes=N_CLASSES,
        n_labels=N_LABELS
    )

    count_df = pd.DataFrame(data)
    count_df = count_df.reindex(index=_sample_names(), columns=_gene_names())

    label_df = pd.DataFrame(labels)
    label_df = label_df.reindex(index=_sample_names())
    return gsf.AnnotatedGEM.from_pandas(count_df.transpose(), label_df, name="Generated GEM")


# Create a factory for creating more than one random GeneSet object.
# See this discussion: https://github.com/pytest-dev/pytest/issues/2703
@pytest.fixture(name="make_gene_set")
def make_gene_set_(random_annotated_gem):

    def make_gene_set(name="Random Gene Set"):
        scores = np.random.random(size=N_GENES)
        xr_scores = xr.Dataset(
            {"scores": (["Gene"], scores),
             "support": (["Gene"], scores > 0.95)},
            coords={"Gene": _gene_names()})
        return gsf.GeneSet(name=name, data=xr_scores)

    yield make_gene_set


# This is the 'normal' singular fixture for a GeneSet.
@pytest.fixture
def random_gene_set(make_gene_set):
    return make_gene_set()


@pytest.fixture
def random_gene_set_collection(random_annotated_gem, make_gene_set):
    gene_sets = {f"set_{x}": make_gene_set(f"set_{x}") for x in range(GENE_SETS)}
    return gsf.GeneSetCollection(gem=random_annotated_gem, gene_sets=gene_sets, name="random collection")
