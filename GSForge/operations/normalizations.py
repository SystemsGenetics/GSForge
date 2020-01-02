"""


"""
import param
import numpy as np
from ..models import OperationInterface

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

__all__ = [
    "ReadsPerKilobaseMillion",
    "UpperQuartile",
]


# class FullQuantile(OperationInterface):
#     """
#
#     """


class UpperQuartile(OperationInterface):
    """
    Under this normalization method, after removing genes having zero read counts for all samples,
    the remaining gene counts are divided by the upper quartile of counts different from zero in
    the computation of the normalization factors associated with their sample and multiplied by
    the mean upper quartile across all samples of the dataset.

    **Original R code [method_compare]_** ::

        uq<-function(X){

          #excluding zero counts in each sample
          UQ<-function(y){
            quantile(y, 0.75)
          }
          X<-X+0.1
          upperQ<-apply(X,2,UQ)
          f.uq<-upperQ/mean(upperQ)
          upq.res<-scale(X,center=FALSE,scale=f.uq)
          return(upq.res)
        }

    .. [method_compare] `A comparison of per sample global scaling and per gene normalization methods for differential
    expression analysis of RNA-seq data <https://doi.org/10.1371/journal.pone.0176185>`_
    """

    # For clarity of testing I write normalization methods as their own, self-contained functions.
    # I often create a numpy and xarray version.
    # @staticmethods decorated functions do not depend on the state of the class, and may be called
    # without instantiating the class.
    @staticmethod
    def np_upper_quartile(counts):
        """
        Perform the upper quartile normalization.

        :param counts: A numpy array containing the raw count values. The shape is assumed to be
            (samples * genes). Zero counts are expected to be present as zeros.

        :returns: The upper quartile normalized count matrix.
        """
        adjusted_counts = counts + 0.1
        per_sample_quantile = np.quantile(adjusted_counts, q=0.75, axis=1)
        per_sample_quantile_mean = np.mean(per_sample_quantile)
        return adjusted_counts / (per_sample_quantile[:, np.newaxis] / per_sample_quantile_mean)

    @staticmethod
    def xr_upper_quartile(counts):
        """
        Perform the upper quartile normalization.

        :param counts: An xarray.DataArray containing the raw count values. The shape is
            assumed to be (samples * genes). Zero counts are expected to be present as zeros.

        :returns: The upper quartile normalized count matrix.
        """
        adjusted_counts = counts + 0.1
        per_sample_quantile = adjusted_counts.quantile(0.75, dim="Gene")
        per_sample_quantile_mean = per_sample_quantile.mean()
        return (adjusted_counts / (per_sample_quantile / per_sample_quantile_mean)).drop("quantile")

    def process(self):
        """
        Perform the upper quartile normalization.

        :return: The upper quartile normalized count matrix.
        """
        counts = self.x_count_data
        return self.xr_upper_quartile(counts)


class ReadsPerKilobaseMillion(OperationInterface):
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
        counts = self.x_count_data
        lengths = self.gem.data[self.length_variable].sel({self.gene_index_name: self.get_gene_index()}).copy()
        return self.xr_reads_per_kilobase_million(counts=counts, lengths=lengths, sample_dim=self.sample_index_name)

# class QuantileNormalization(OperationInterface):
#
#     quantile = param.Magnitude(default=0.75)
#
#     @staticmethod
#     def np_quantile_normalization():
#         pass
#
#     @staticmethod
#     def xr_quantile_normalization():
#         pass
#
#     def process(self):
#         self.set_param(count_mask="dropped")
#         counts = self.x_data
#         sample_quantile = counts.quantile(q=self.quantile, dim=self.sample_index_name)
#         return counts / (sample_quantile / sample_quantile.mean())
#
#
# class NormalizationMethod(OperationInterface):
#     """
#     Some normalization method.
#     """
#     k = param.Integer(default=100)
#
#     @staticmethod
#     def normalization_method(*args, **kwargs):
#         raise NotImplemented
#
#     def process(self):
#         random_genes = np.random.choice(self.gene_index, self.k)
#         return random_genes, {"random_size": self.k}
