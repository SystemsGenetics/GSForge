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

# Sample-wise Count Distributions

View the distribution of count intensities on a per-sample basis.

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
BOR_COLL_PATH = OSF_PATH.joinpath("Collections", "boruta")
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

### Creating a Sample-wise Distribution plot

```{code-cell}
gsf.plots.gem.SamplewiseDistributions(
    agem,
    selected_samples=agem.data.Sample[:25],
    hue_key="treatment",
    count_transform=lambda counts: np.log2(counts.where(counts > 0))
)
```
