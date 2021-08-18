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

# Interface Guide

Set operations, intersections, unions and differences and combinations thereof answer many of the most basic questions.
Moreover, we often desires to examine subsets of the original GEM based a selection determined by some set operation.
The ``Interface`` provided by ``GSForge`` provides a uniform API access to both the ``AnnotatedGEM`` and the
``GeneSetCollection`` objects for retrieving count values and sample annotations.


***Notebook setup***

```{code-cell}
import numpy as np
import pandas as pd
import xarray as xr
import GSForge as gsf
import holoviews as hv
hv.extension('bokeh')

# OS-independent path management.
from os import fspath, environ
from pathlib import Path

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/")).expanduser().joinpath("osfstorage", "oryza_sativa")
GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hisat2_raw.nc")
TOUR_DGE = OSF_PATH.joinpath("GeneSetCollections", "DEG_gene_sets")
```

## Load Data

```{code-cell}
agem = gsf.AnnotatedGEM(GEM_PATH)
agem 
```

```{code-cell}
agem.data
```

## Selecting Data using the Interface

Select data through the interface via the `get_gem_data` function.

The simplest possible call to this function returns the zero filled count matrix within a two-item tuple.

```{code-cell}
counts, empty_labels = gsf.get_gem_data(agem)
```

```{code-cell}
counts
```

When `annotation_variables` are not provided, the second item is the `None` singleton.

```{code-cell}
empty_labels == None
```

This is for the sake of consistency in handling the number of expected output objects when `annotation_variables` 
are provided:

```{code-cell}
counts, label_ds = gsf.get_gem_data(agem, annotation_variables=['treatment'])
```

`counts` remain the same, but this time a corresponding label dataset is also returned.

```{code-cell}
label_ds
```

### Masking or Dropping Samples or Genes

Count values can be returned in one of three 'forms' with respect to NaN or zero values:
+ zero counts 'masked' as `NaN`.
+ zero counts 'dropped' gene-wise.
+ zero counts 'complete', or included within the matrix, this is the default selection.

This is done using the `count_mask` argument.

```{code-cell}
counts, _ = gsf.get_gem_data(agem, count_mask='complete')
counts.shape
```

```{code-cell}
counts, _ = gsf.get_gem_data(agem, count_mask='dropped')
counts.shape
```

For samples there are only two options:
+ Samples with a missing annotation are 'dropped'.
+ Use the 'complete' set of samples.

This only has an effect when used with `annotation_variables`.

```{code-cell}
counts, label_ds = gsf.get_gem_data(agem, annotation_mask='complete', annotation_variables=['treatment'])
```

### Transforming the Count Matrix

A transform can be applied when selecting data, such as a log transform. ``GSFoge`` allows users to supply a function 
to transform the subset of counts returned, this function only operates on the count values returned.

```{code-cell}
counts, label_ds = gsf.get_gem_data(agem, annotation_mask='complete', annotation_variables=['treatment'],
                                    count_transform=lambda c: np.log(c + 1.0))
```
