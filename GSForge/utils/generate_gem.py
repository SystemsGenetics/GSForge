from sklearn.datasets import make_multilabel_classification
import pandas as pd
import numpy as np
import param

from ..models import AnnotatedGEM


class RandomGEMSimulator(param.Parameterized):

    def __init__(self, **params):
        super().__init__(**params)
        # Prepare random number generator.
        self.rng = np.random.default_rng()

    def build_mu_array(self):
        basis_set = np.array([-3.0, 7.0, 0.0])
        self.rng.choice(basis_set, 5000)

        # log2_base_gene_means = np.hstack((
        #     self.rng.normal(loc=-3.0, scale=2.0, size=2000),
        #     self.rng.normal(loc=7.0, scale=2.0, size=2000),
        #     self.rng.normal(loc=0.0, scale=2.0, size=1000),
        # ))
        # return np.exp2(log2_base_gene_means)


    def build_count_matrix(self):
        pass
