"""
Normalization functions inherit from the ``OperationInterface`` class.
This means that they can all be called upon an ``AnnotatedGEM`` or a ``GeneSetCollection``.

These (classes) functions have static methods that implement the transform on a ``numpy``
or ``xarray`` source.
"""

import numpy as np
import param

from ..models import Interface

__all__ = [
    "ReadsPerKilobaseMillion",
    "UpperQuartile",
]


class UpperQuartile(Interface, param.ParameterizedFunction):
    """
    Under this normalization method, after removing genes having zero read counts for all samples,
    the remaining gene counts are divided by the upper quartile of counts different from zero in
    the computation of the normalization factors associated with their sample and multiplied by
    the mean upper quartile across all samples of the dataset. [method_compare]_

    Original R code.  ::

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

    .. [method_compare] `A comparison of per sample global scaling and per gene normalization methods for differential expression analysis of RNA-seq data <https://doi.org/10.1371/journal.pone.0176185>`_

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
            (samples by genes). Zero counts are expected to be present as zeros.

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

        :param counts: An ``xarray.DataArray`` containing the raw count values. The shape is
            assumed to be (samples by genes). Zero counts are expected to be present as zeros.

        :returns: The upper quartile normalized count matrix.
        """
        adjusted_counts = counts + 0.1
        per_sample_quantile = adjusted_counts.quantile(0.75, dim="Gene")
        per_sample_quantile_mean = per_sample_quantile.mean()
        return (adjusted_counts / (per_sample_quantile / per_sample_quantile_mean)).drop("quantile")

    def __call__(self, *args, **params):
        """
        Perform the upper quartile normalization.

        :return: The upper quartile normalized count matrix.
        """
        super().__init__(*args, **params)
        counts = self.x_count_data
        return self.xr_upper_quartile(counts)


class ReadsPerKilobaseMillion(Interface, param.ParameterizedFunction):
    """
    RPKM or FPKM -- Reads or Fragments per per Kilobase Million.

    These methods attempt to compensate for sequencing depth and gene length.
    The utility of this method is disputed in the literature [cite me].
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

    def __call__(self, *args, **params):
        super().__init__(*args, **params)
        counts = self.x_count_data
        lengths = self.gem.data[self.length_variable].sel({self.gene_index_name: self.get_gene_index()}).copy()
        return self.xr_reads_per_kilobase_million(counts=counts, lengths=lengths, sample_dim=self.sample_index_name)
