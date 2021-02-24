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

# GeneSet Overlap Heatmap

Here we use a heatmap to visualize the overlap between gene sets within a single collection.

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
LIT_DGE_GSC_PATH = OSF_PATH.joinpath("Collections", "literature", "DGE")
```

```{code-cell}
agem = gsf.AnnotatedGEM(HYDRO_GEM_PATH)
agem
```

```{code-cell}
gsc = gsf.GeneSetCollection.from_folder(
    gem=agem, target_dir=LIT_DGE_GSC_PATH, name="DGE Results")
gsc
```

## Viewing GeneSet Overlaps within a GeneSetCollection

```{code-cell} ipython3
gsf.plots.collections.WithinCollectionOverlapHeatMap(gsc)
```
