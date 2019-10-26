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
from ..utils import kwargs_overlap


class ReadsPerKilobaseMillion:
    """

    ***Citation(s):***
    +[**]()

    """

    length_variable = param.String(default="lengths")

    @staticmethod
    def reads_per_kilobase_million(counts, lengths, sample_dim="Sample"):
        scaling_factor = counts.sum(dim=sample_dim) / 1e6
        reads_per_million = counts / scaling_factor
        normalized_counts = reads_per_million / lengths
        return normalized_counts


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

