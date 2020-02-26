import xarray as xr

__all__ = [
    "GSFAccessor"
]


@xr.register_dataarray_accessor('gem')
class GSFAccessor:

    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    @property
    def masked_genes(self):
        return self._obj.where(self._obj > 0)

    @property
    def dropped_genes(self):
        return self.masked_genes().dropna(dim="Gene")
