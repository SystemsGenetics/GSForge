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

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage")).expanduser()
COUNT_PATH = OSF_PATH.joinpath("GEMmaker_GEMs", "Osativa_heat_drought_PRJNA301554.GEM.raw.txt")
LABEL_PATH = OSF_PATH.joinpath("raw_annotation_data", "GEMmaker_GEMs", "PRJNA301554.hydroponic.annotations.txt")

# Output path.
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_raw.nc")
```

## From Text Files

### My counts are in a single file

Probably the easiest scenario to deal with.

```{code-cell} ipython3
count_df = pd.read_csv(COUNT_PATH, sep="\t", index_col=0)
```

### My counts are in multiple files

```{code-cell} ipython3

```
