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

# Creating an Annotated Expression Matrix

Here we demonstrate the most common ways to load expression data into `GSForge`.

***Notebook setup***

```{code-cell} ipython3
from os import environ
from pathlib import Path
import numpy as np
import pandas as pd
import GSForge as gsf

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage/oryza_sativa")).expanduser()
COUNT_PATH = OSF_PATH.joinpath("GEMmakerGEMs", "rice_heat_drought.GEM.raw.txt")
LABEL_PATH = OSF_PATH.joinpath("GEMmakerGEMs", "raw_annotation_data", "PRJNA301554.hydroponic.annotations.txt")

# Output path.
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hisat2_raw.nc")
```

## From Text Files

If your count and annotation files have matching sample indices you can create an AnnotatedGEM in a single step:

```{code-cell} ipython3
agem = gsf.AnnotatedGEM.from_files(
    count_path=COUNT_PATH,
    label_path=LABEL_PATH,
    # These are the default arguments passed to from_files,
    # to the individual calls to `pandas.read_csv`.
    count_kwargs=dict(index_col=0, sep="\t"),
    label_kwargs=dict(index_col=1, sep="\t"),
)
```

It is not uncommon to have to wrangle sample or gene names to some degree.
Once complete you may supply the data as a pair of `pandas.Dataframe` or a single `xarray.DataSet` object.

```{code-cell} ipython3
count_df = pd.read_csv(COUNT_PATH, sep="\t", index_col=0)
# Wrangle data here...

label_df = pd.read_csv(LABEL_PATH, index_col=1, sep="\t")
label_df['genotype'] = label_df['genotype'].str.split(" ", expand=True).iloc[:, 0]
label_df['time'] = label_df['time'].str.split(' ', expand=True).iloc[:, 0].astype(int)
# Perhaps even more wrangling...

# Then provide them to GSForge:
del agem
agem = gsf.AnnotatedGEM.from_pandas(count_df=count_df, label_df=label_df, name="Oryza sativa")

if not GEM_PATH.exists():
    agem.save(GEM_PATH)
    
agem
```

