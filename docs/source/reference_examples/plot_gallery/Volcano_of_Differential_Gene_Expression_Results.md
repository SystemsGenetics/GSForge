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

# Volcano Plot of Differential Gene Expression Results

A common way to visualize the results of a DGE analysis.

+++

***Plotting Guide Setup***

Shared setup for all plotting guides.

```{code-cell} ipython3
# OS-independent path management.
from os import environ
from pathlib import Path
import numpy as np
import GSForge as gsf
import holoviews as hv
hv.extension('bokeh')

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage/oryza_sativa")).expanduser()
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hisat2_raw.nc")
TOUR_DGE = OSF_PATH.joinpath("GeneSetCollections", "DEG_gene_sets")
```

```{code-cell} ipython3
agem = gsf.AnnotatedGEM(GEM_PATH)
agem
```

```{code-cell} ipython3
deg_gsc = gsf.GeneSetCollection.from_folder(
    gem=agem, target_dir=TOUR_DGE, name="DGE Results")
deg_gsc
```

***Select a particular result set of interest***

```{code-cell} ipython3
deg_gs = deg_gsc.gene_sets["'0 + treatment:genotype'__treatment[HEAT]"]
deg_gs
```

***View the data stored within this GeneSet result***

```{code-cell} ipython3
deg_gs.data
```

## Create a Volcano Plot

```{code-cell} ipython3
gsf.plots.results.Volcano(deg_gs)
```

```{code-cell} ipython3

```
