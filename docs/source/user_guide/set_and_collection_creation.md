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

***Notebook setup***

```{code-cell} ipython3
from os import environ
from pathlib import Path
import numpy as np
import pandas as pd
import GSForge as gsf
rng = np.random.default_rng(0)
```

# Creating Feature Sets and Collections

In this example we have an annotated GEM already constructed:

```{code-cell} ipython3
OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/")).expanduser().joinpath("osfstorage", "oryza_sativa")
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hisat2_raw.nc")
```

```{code-cell} ipython3
agem = gsf.AnnotatedGEM(GEM_PATH)
agem
```

## Creating GeneSets

See the API reference of [`GSForge.GeneSet`](../API/GSForge.models) for all availble creation functions, which include:
+ `from_pandas`
+ `from_GeneSets`
+ `from_bool_array`
+ `from_gene_array`
+ `from_xarray_dataset`
+ `from_netcdf`


### GeneSets from Lists or Arrays

A minimal GeneSet is just that: a set of genes.

Here we draw random features to demonstrate this:

```{code-cell} ipython3
random_features = rng.choice(agem.data.Gene.values, 10, replace=False)
random_features
```

Provide this to the `from_gene_array()` constructor to create a simple GeneSet.

```{code-cell} ipython3
set_example = gsf.GeneSet.from_gene_array(random_features, name='Random Set')
set_example
```

### From `pandas.DataFrames` or `xarray.DataSet` objects

Commonly there is some information associated with a set of features.
Differential gene expression results often contain information about many (or all) of the genes, but only identify a few as 'differentially expressed'.
We can store all such information in the `GeneSet` object, and indicate which genes are selected by setting a boolean array named `support`.

Here we simulate an example `DataFrame`.

```{code-cell} ipython3
n = 100
random_features = rng.choice(agem.data.Gene.values, n, replace=False)

df = pd.DataFrame(
    {
        'sim_LFC': rng.normal(size=n),
        'sim_pvalue': np.abs(rng.normal(size=n)),
    },
    index=random_features
)

df['support'] = (df['sim_pvalue'] < 0.05) | (df['sim_LFC'] > 1.0)

df.head()
```

Create the `GeneSet` using the `from_pandas()` constructor.

```{code-cell} ipython3
dge_gs = gsf.GeneSet.from_pandas(df, name='Sim DGE')
dge_gs
```

```{code-cell} ipython3
dge_gs.data
```

## Creating and Saving GeneSetCollections

We only need to provided an AnnotatedGEM and a name to create a `GeneSetCollection`.
Then add `GeneSet` objects like you would entries to a dictionary:

```{code-cell} ipython3
sample_coll = gsf.GeneSetCollection(gem=agem, name='Literature DGE')
sample_coll['set example'] = dge_gs
sample_coll['simulated DGE example'] = set_example
sample_coll
```

GeneSetCollections are saved as a directory, each set saved as a separate netcdf file.

```{code-cell} ipython3
save = False
if save == True:
    sample_coll.save('my_collection')
```
