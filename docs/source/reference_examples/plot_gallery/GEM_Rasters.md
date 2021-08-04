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

# GEM Rasters

Here we display gene intensity values as an image.
Using the [datashader](https://datashader.org/) library,
we can create a dynamic raster of the entire GEM.

Specific results from this plot are not the goal,
rather hope to *not* see strange patterns in expression or missing values that may indicate
a problem in a previous processing step.

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

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage/oryza_sativa")).expanduser()
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_raw.nc")
TOUR_BORUTA = OSF_PATH.joinpath("GeneSetCollections", "tour_boruta")
```

```{code-cell}
agem = gsf.AnnotatedGEM(GEM_PATH)
agem
```

```{code-cell}
gsc = gsf.GeneSetCollection.from_folder(
    gem=agem, target_dir=TOUR_BORUTA, name="Boruta Results")
gsc
```

## Creating a Count Matrix Raster

```{code-cell}
gsf.plots.gem.RasterGEM(agem)
```

Selecting a `GeneSet` from a collection provides a raster of just those supported genes.

```{code-cell}
gsf.plots.gem.RasterGEM(gsc, selected_gene_sets=["Boruta_treatment"], hue="genotype")
```
