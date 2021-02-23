---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Plotting Guide

This notebook describes the types of plotting functions included with `GSForge`,
as well as how to apply customizations to those plots.

Plotting functions are delineated by the (primary) data source:

+ **gem** the count array from an `AnnotatedGEM` or the same via a `GeneSetCollection`.
+ **collections** based on membership of `GeneSet` objects within a `GeneSetCollection`.
+ **results** Custom functions for specific analytical methods, e.g. volcano plots for DGE analysis.

`GSForge` uses the `Holoviews` package for creating plots.
`Holoviews` is a common API to create plots using two popular backends, `matplotlib` and `bokeh`.

+++

***Plotting Guide Setup***

Shared setup for all plotting guides.

```{code-cell} ipython3
# OS-independent path management.
from os import environ
from pathlib import Path
import GSForge as gsf
import holoviews as hv
hv.extension('bokeh', 'matplotlib')

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage")).expanduser()
HYDRO_GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hydro_raw.nc")
HYDRO_NORMED_GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hydro_normed.nc")
BOR_COLL_PATH = OSF_PATH.joinpath("Collections", "boruta")
```

```{code-cell} ipython3
agem = gsf.AnnotatedGEM(HYDRO_GEM_PATH)
agem
```

```{code-cell} ipython3
---
pycharm:
  name: '#%%

    '
---
gsc = gsf.GeneSetCollection.from_folder(
    gem=agem, target_dir=BOR_COLL_PATH, name="Boruta Results")
gsc
```

Select a backend, or simply only load the one you wish to use.

```{code-cell} ipython3
hv.output(backend="matplotlib")
gsf.plots.gem.GenewiseAggregateScatter(agem, datashade=False, apply_default_opts=False)
```

```{code-cell} ipython3
hv.output(backend="bokeh")
gsf.plots.gem.GenewiseAggregateScatter(agem, datashade=True)
```

### Saving Plots

```python
hv.save(my_plot.options(toolbar=None),  # This sometimes does not remove the toolbar.
        "my_filename.png")
```

+++

**Where to go for more help?**

See the excellent [Holoviews documentation](http://holoviews.org/) for more on creating and styling plots.

+++

**How to return a plot with no styling applied**

It is then up to you to modify the `holoviews` object, and to apply desired options from your selected backend.

All `GSForge` plotting functions have an `apply_default_opts`, which attempts to apply options based on
the currently loaded backend extension.
If you want to apply your own styling, plots can be returned with their default settings only by setting
``apply_default_opts=False``.
