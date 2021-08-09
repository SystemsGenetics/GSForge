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

# Gene vs Counts

Here genes form the x-axis, with intensity values on the y-axis.
This method is intended to be used to compare a subset of selected genes,
and will give a warning if you try to pass too many genes.
This can be overridden by setting ``soft_max=<some integer>``.

***Plotting Guide Setup***

Shared setup for all plotting guides.

```{code-cell}
# OS-independent path management.
from os import environ
from pathlib import Path
import numpy as np
import GSForge as gsf
import holoviews as hv
hv.extension('bokeh')

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/")).expanduser().joinpath("osfstorage", "oryza_sativa")
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hisat2_raw.nc")
TOUR_BORUTA = OSF_PATH.joinpath("GeneSetCollections", "tour_boruta")
```

```{code-cell}
agem = gsf.AnnotatedGEM(GEM_PATH)
agem
```

```{code-cell}
gsc = gsf.GeneSetCollection.from_folder(
    gem=agem, target_dir=TOUR_BORUTA, name="Boruta Tour Results")
gsc
```

## Plot Genes vs Gene Counts

**Select Genes of Interest**

There are a few too many genes for us to look at all the genes in our collection.
Rather we will get the 10 largest (by mean) intensity genes.

```{code-cell}
counts, _ = gsf.get_gem_data(gsc, selected_gene_sets=["Boruta_treatment"])
selected_genes = counts.isel(Gene=np.argsort(counts.mean(dim="Sample").values)[-10:])["Gene"].values
selected_genes
```

```{code-cell}
plot = gsf.plots.gem.GeneVsCountsScatter(
    gsc, selected_genes=selected_genes, hue="treatment",
    count_transform=lambda counts: np.log2(counts.where(counts > 0))
)
plot
```
