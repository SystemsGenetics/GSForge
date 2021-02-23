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

# Grouped Gene Covariance

Here we group samples by an annotation category, then plot the gene means against one another.

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


## View Annotation Group Covariance

```{code-cell}
gsf.plots.gem.GroupedGeneCovariance(
    agem, group_variable="treatment", 
    x_group_label="CONTROL", y_group_label="HEAT",
    count_transform=lambda counts: np.log10(counts + 0.25))
```
