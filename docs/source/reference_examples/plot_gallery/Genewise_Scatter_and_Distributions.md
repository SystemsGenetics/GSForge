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

# Gene-wise Aggregate Scatter and Distribution Plots

Displays the output of selected aggregations upon the count array on a scatter plot with optional
adjoined kernel density estimates. e.g. mean counts vs mean variance etc. By default, such outputs are
$log_2$ transformed as well.

The available axis aggregation functions are:
    + frequency
    + mean
    + variance
    + standard_dev
    + fano
    + mean_rank
    + cv_squared

The default axis selections are 'mean' vs 'variance'.

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
## Creating a Gene-wise scatter plot

```{code-cell}
gsf.plots.gem.GenewiseAggregateScatter(agem, datashade=True)
```

## Select an Aggregation method

```{code-cell}
gsf.plots.gem.GenewiseAggregateScatter(agem, y_axis_selector="cv_squared",  datashade=True)
```

## Select and Group Genes

We can also limit and select to sets of genes.
Note that the aggregation transform occurs on the `x_count_data`.
Selecting more than one gene set returns the union of those sets by default.

```{code-cell}
gsf.plots.gem.GenewiseAggregateScatter(
    gsc, datashade=False, selected_gene_sets=["treatment", "genotype"])
```
