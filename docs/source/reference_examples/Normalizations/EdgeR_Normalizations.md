---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.2
kernelspec:
  display_name: gsfenv
  language: python
  name: gsfenv
---

# Using EdgeR to create Normalizations

This notebook provides an example of using data selected from `GSForge` in an R session, 
and returning the results back to Python to store and view the results with `GSForge`.

```{code-cell} ipython3
import xarray as xr
import GSForge as gsf
import holoviews as hv
hv.extension('bokeh')
```

***R integration setup***

```{code-cell} ipython3
import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
%load_ext rpy2.ipython
pandas2ri.activate()
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Supresses verbose R output.
```

```{code-cell} ipython3
%%R
library("edgeR")
```

***Declare paths used***

```{code-cell} ipython3
# OS-independent path management.
from os import fspath, environ
from pathlib import Path
OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage/oryza_sativa")).expanduser()
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hisat2_raw.nc")
```

***Load an AnnotatedGEM***

```{code-cell} ipython3
agem = gsf.AnnotatedGEM(GEM_PATH)
agem
```

***Select counts using `get_gem_data()`***

```{code-cell} ipython3
counts, _ = gsf.get_gem_data(agem)
```

***Prepare the counts for R***

Notice the counts are transposed after this step to the form more common in R. (features by samples).

```{code-cell} ipython3
ri_counts = gsf.utils.R_interface.Py_counts_to_R(counts)
ri_counts.shape
```

***Run the normalization within R***

```{code-cell} ipython3
%%R -i ri_counts -o tmm_counts

dge_list <- DGEList(counts=ri_counts)
dge_list <- calcNormFactors(dge_list, method="TMM")
tmm_counts <- cpm(dge_list, normalized.lib.sizes=TRUE, log=FALSE)
```

```{code-cell} ipython3
tmm_counts = xr.DataArray(tmm_counts.T, coords=counts.coords, name='tmm_counts')
tmm_counts
```

***Add the counts to the GEM .data attribute.***

```{code-cell} ipython3
agem.data['tmm_counts'] = tmm_counts
```
