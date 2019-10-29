"""

Normalization Methods
=====================

**A replication of Table 1 from https://doi.org/10.1371/journal**

Summary of current normalization methods to correct the technical biases for RNA-Seq data.

Total count (TC)
Median
UQ Upper Quartile
TMM Trimmed Mean ofM-values
RLE Relative Log Expression
Quantile (Q)
PCA Principal Component Analysis
RUV Remove Unwanted Variation
SVA Surrogate Variable Analysis
TPM (transcripts per million)
FPKM / RPKM (reads/fragments per kilo-base per million mapped reads)

"""

import param
import numpy as np
from ..models import OperationInterface
from .analytics import get_data
from ..utils import kwargs_overlap

__all__ = [
    "ReadsPerKilobaseMillion"
]


class ReadsPerKilobaseMillion(OperationInterface):
    """

    ***Citation(s):***
    +[**]()

    """

    length_variable = param.String(default="lengths")

    @staticmethod
    def xr_reads_per_kilobase_million(counts, lengths, sample_dim="Sample"):
        scaling_factor = counts.sum(dim=sample_dim) / 1e6
        reads_per_million = counts / scaling_factor
        normalized_counts = reads_per_million / lengths
        return normalized_counts

    @staticmethod
    def np_reads_per_kilobase_million(counts, lengths):
        scaling_factor = counts.sum(axis=0) / 1e6
        reads_per_million = counts / scaling_factor
        normalized_counts = reads_per_million / lengths
        return normalized_counts

    def process(self):
        counts = self.x_data
        lengths = self.gem.data[self.length_variable].sel({self.gene_index_name: self.get_gene_index()}).copy()
        return self.xr_reads_per_kilobase_million(counts=counts, lengths=lengths, sample_dim=self.sample_index_name)


class QuantileNormalization(OperationInterface):

    quantile = param.Magnitude(default=0.75)

    @staticmethod
    def np_quantile_normalization():
        pass

    @staticmethod
    def xr_quantile_normalization():
        pass

    def process(self):
        self.set_param(count_mask="dropped")
        counts = self.x_data
        sample_quantile = counts.quantile(q=self.quantile, dim=self.sample_index_name)
        return counts / (sample_quantile / sample_quantile.mean())


class NormalizationMethod(OperationInterface):
    """
    Some normalization method.
    """
    k = param.Integer(default=100)

    @staticmethod
    def normalization_method(*args, **kwargs):
        raise NotImplemented

    def process(self):
        random_genes = np.random.choice(self.gene_index, self.k)
        return random_genes, {"random_size": self.k}
