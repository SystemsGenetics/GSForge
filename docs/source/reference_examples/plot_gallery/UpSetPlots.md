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

# UpSetPlots

```{code-cell} ipython3
from os import  environ
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import GSForge as gsf
import upsetplot
import matplotlib.pyplot as plt
import holoviews as hv
%matplotlib inline
hv.extension('bokeh')

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage/oryza_sativa")).expanduser()
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hisat2_raw.nc")
LIT_DGE_GSC_PATH = OSF_PATH.joinpath("GeneSetCollections", "literature", "DGE")
LIT_TF_PATH = OSF_PATH.joinpath("GeneSetCollections", "literature", "TF")
BORUTA_GSC_PATH = OSF_PATH.joinpath("GeneSetCollections", "boruta")
```

```{code-cell} ipython3
agem = gsf.AnnotatedGEM(GEM_PATH)
agem
```

```{code-cell} ipython3
%%time
lit_dge_coll = gsf.GeneSetCollection.from_folder(gem=agem, target_dir=LIT_DGE_GSC_PATH, name="Literature DGE")
lit_tf_coll = gsf.GeneSetCollection.from_folder(gem=agem, target_dir=LIT_TF_PATH, name="Literature TF")
boruta_gsc = gsf.GeneSetCollection.from_folder(gem=agem, target_dir=BORUTA_GSC_PATH, name="Boruta Results")
tf_geneset = gsf.GeneSet.from_GeneSets(*list(lit_tf_coll.gene_sets.values()), name='transcription factors')
```

```{code-cell} ipython3
combined_gsc = gsf.GeneSetCollection(gem=agem, 
                                     gene_sets={**lit_dge_coll.gene_sets, 
                                                **boruta_gsc.gene_sets, 
                                                'transcription factors': tf_geneset})
```

```{code-cell} ipython3
gsf.plots.collections.WithinCollectionOverlapHeatMap(combined_gsc)
```

```{code-cell} ipython3
gsf.plots.collections.UpsetPlotInterface(combined_gsc, min_overlap_size=2, 
                                         upset_kwargs=dict(orientation='vertical'))
```

```{code-cell} ipython3
data = gsf.plots.collections.UpsetPlotInterface.build_membership_series(combined_gsc.as_dict(), min_size=5)
# Create the figure at the right size and resolution.
# fig_inches = 3.5
fig_dpi = 300  # 300 is a common DPI requirement.

# Construct the figure.
fig, ax = plt.subplots(dpi=fig_dpi)
ax.axis('off')
upsetplot.plot(data, fig=fig, orientation='vertical')
```

```{code-cell} ipython3
data = gsf.plots.collections.UpsetPlotInterface.build_membership_series(lit_dge_coll.as_dict(), min_size=5)
data
# upsetplot.UpSet(data=data, orientation='vertical')
```

```{code-cell} ipython3
upsetplot.UpSet(data=data)
```
