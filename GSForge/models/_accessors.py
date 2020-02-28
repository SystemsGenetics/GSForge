# import xarray as xr
#
# __all__ = [
#     "CountMatrixAccessor"
# ]
#
#
# @xr.register_dataset_accessor('gem')
# class CountMatrixAccessor:
#     """
#     http://xarray.pydata.org/en/stable/internals.html
#     """
#
#     def __init__(self, xarray_obj):
#         self._obj = xarray_obj
#
#     @property
#     def complete_counts(self):
#         count_variable = self.gem.data.attrs.get("count_array_name", "counts")
#         return self._obj[count_variable].fillna(0)
#
#     @property
#     def masked_counts(self):
#         count_variable = self.gem.data.attrs.get("count_array_name", "counts")
#         counts = self._obj[count_variable]
#         return counts.where(counts > 0.0)
#
#     @property
#     def complete_genes(self):
#         gene_dim = self.gem.data.attrs.get("gene_index_name", "Gene")
#         return self.complete_counts[gene_dim].values.copy()
#
#     @property
#     def masked_genes(self):
#         count_variable = self.gem.data.attrs.get("count_array_name", "counts")
#         sample_dim = self.gem.data.attrs.get("sample_index_name", "Sample")
#         gene_dim = self.gem.data.attrs.get("gene_index_name", "Gene")
#         return self._obj.where(self._obj > 0)
#
#     @property
#     def dropped_genes(self):
#         return self.masked_genes().dropna(dim="Gene")
