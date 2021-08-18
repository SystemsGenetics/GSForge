from sklearn.datasets import make_multilabel_classification
import pandas as pd
import numpy as np
import param

from ..models import AnnotatedGEM


class RandomGEMSimulator(param.Parameterized):
    """
    
    Based on the publication from Freytag et al. BMC Bioinformatics (2015) 16:309.
    
    doi.org/10.1371/journal.pbio.3000481

    
    """
    n_samples = param.Integer()
    n_genes = param.Integer()
    
    n_negative_control_genes = param.Integer()
    n_strongly_expressed_genes = param.Integer()
    
    n_signal_dims = param.Integer()
    n_random_noise_dims = param.Integer()
    n_signal_noise_dims = param.Integer()
    
    random_noise_sigma = param.Parameter()


    def __init__(self, **params):
        super().__init__(**params)
        # Prepare random number generator.
        self.rng = np.random.default_rng()
        
    def beta_matrix(self):
        return np.matrix(np.hstack((
            rng.uniform(
                low= -2 / np.sqrt(self.n_signal_dims), 
                high= 2 / np.sqrt(self.n_signal_dims), 
                size=(self.n_signal_dims,self.n_genes - self.n_negative_control_genes)),
            np.zeros((self.n_signal_dims, self.n_negative_control_genes))
        )))
    
    def noise_matrix(self):
        return np.block([
            [np.eye(self.n_signal_noise_dims), 
             np.zeros((self.n_signal_noise_dims, self.n_random_noise_dims - self.n_signal_noise_dims))],
            [np.zeros((self.n_signal_dims - self.n_signal_noise_dims, self.n_random_noise_dims))],
        ])
    
    def noise_and_design_correlation_matrix(self):
        noise_matrix = self.noise_matrix()
        return np.block([
            [np.eye(self.n_signal_dims), noise_matrix],
            [noise_matrix.T, np.eye(self.n_random_noise_dims)],
        ])


    def build_count_matrix(self):
        
        gene_signal_matrix = np.matrix(
            rng.multivariate_normal(mean=np.zeros(self.n_signal_dims + self.n_random_noise_dims), 
                                    size=self.n_samples, 
                                    cov=self.noise_and_design_correlation_matrix()))
        X_signal = np.matrix(gene_signal_matrix[:, :signal_dims])
        noise_signal = np.matrix(gene_signal_matrix[:, signal_dims:])
        
        alpha = np.matrix(
            rng.uniform(low=-2 * self.random_noise_sigma / np.sqrt(self.n_random_noise_dims), 
                        high=2 * self.random_noise_sigma / np.sqrt(self.n_random_noise_dims), 
                        size=(self.n_random_noise_dims, self.n_genes)))
            
        XB = X_signal * self.beta_matrix()
        W_alpha = noise_signal * alpha
            
        Y = XB + W_alpha + noise
        
        return Y
