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

# DGE Log-Fold Change vs Mean

A common way to visualize the results of a DGE analysis.

***Plotting Guide Setup***

A shared setup for all plotting guides.

```{code-cell}
# OS-independent path management.
from os import environ
from pathlib import Path
import numpy as np
import GSForge as gsf
import holoviews as hv
hv.extension('bokeh')

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage")).expanduser()
HYDRO_GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hydro_raw.nc")
HYDRO_NORMED_GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hydro_normed.nc")
BOR_COLL_PATH = OSF_PATH.joinpath("Collections", "nf_boruta")
```

```{code-cell}
agem = gsf.AnnotatedGEM(HYDRO_GEM_PATH)
agem
```

```{code-cell}
gsc = gsf.GeneSetCollection.from_folder(
    gem=agem, target_dir=BOR_COLL_PATH, name="Boruta Results")
gsc
```

***Load Differential Gene Expression Analysis Results into a `GeneSetCollection`***

```{code-cell}
deg_gsc = gsf.GeneSetCollection.from_folder(gem=agem, target_dir=LIT_DGE_GSC_PATH, name="DEG Results")
deg_gsc
```

***Select a particular result set of interest***

```{code-cell}
deg_gs = deg_gsc.gene_sets["DROUGHT_UP"]
deg_gs
```

***View the data stored within this GeneSet result***

```{code-cell}
deg_gs.data
```

## Plot gene means vs log-fold change

```{code-cell} ipython3
gsf.plots.results.MeanVsLFC(deg_gs)
```
