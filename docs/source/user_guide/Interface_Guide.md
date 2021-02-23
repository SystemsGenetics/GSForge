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

# Interface Guide

+++

This notebook describes interactions with the `Interface` object of GSForge, the primary way in which users access count and annotation data.

The `Interface` class is a common application interface (API) for interacting with `AnnotatedGEM` and `GeneSetCollection` objects.

+++

## Notebook Setup

### Python Imports

```{code-cell} ipython3
import numpy as np
import pandas as pd
import xarray as xr

import GSForge as gsf

import holoviews as hv
hv.extension('bokeh')
```

+++ {"jupyter": {"outputs_hidden": false}, "pycharm": {"is_executing": false, "name": "#%%\n"}}

### Data Path Management

```{code-cell} ipython3
# OS-independent path management.
from os import fspath, environ
from pathlib import Path

OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data")).expanduser()
AGEM_PATH = OSF_PATH.joinpath("osfstorage", "rice.nc")
DEG_COLL_PATH = OSF_PATH.joinpath("osfstorage", "DEG_gene_sets")
assert AGEM_PATH.exists()
```

## Logging Setup

Logging output can be enabled so as to clarify the steps being preformed by `gsforge`.
This can set for an entire notebook or set on a per-function basis (where supported).

```{code-cell} ipython3
# This cell enables logging for gsforge within this notebook.
import logging
import sys

logger = logging.getLogger("GSForge")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter('%(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)
```

In other cases, you can try passing a log level the `verbose` argument.
The available levels are: `['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG']`.

**New users are encouraged to set the log-level to** `'INFO'` **notebook-wide as above.**

+++

## Load Data

```{code-cell} ipython3
agem = gsf.AnnotatedGEM(AGEM_PATH)
agem 
```

```{code-cell} ipython3
agem.data
```

## Selecting Data using the Interface

Selecting data through the interface is done via the `get_gem_data` function.

The simplest possible call to this function returns the zero filled count matrix within a two-item tuple.

```{code-cell} ipython3
counts, empty_labels = gsf.get_gem_data(agem)
```

```{code-cell} ipython3
counts
```

When `annotation_variables` are not provided, the second item is the `None` singleton.

```{code-cell} ipython3
empty_labels == None
```

This is for the sake of consistency in handling the number of expected output objects when `annotation_variables` are provided:

```{code-cell} ipython3
counts, label_ds = gsf.get_gem_data(agem, annotation_variables=['Treatment'])
```

`counts` remain the same, but this time a corresponding label dataset is also returned.

```{code-cell} ipython3
label_ds
```

### Masking or Dropping Samples or Genes

It is often desired that the returned count values be returned in one of three 'forms' with respect to NaN or zero values:
+ zero counts 'masked' as `NaN`.
+ zero counts 'dropped' gene-wise.
+ zero counts 'complete', or included within the matrix, this is the default selection.

This is done using the `count_mask` argument.

```{code-cell} ipython3
counts, _ = gsf.get_gem_data(agem, count_mask='complete')
counts.shape
```

```{code-cell} ipython3
counts, _ = gsf.get_gem_data(agem, count_mask='dropped')
counts.shape
```

For samples there are only two options:
+ Samples with a missing annotation are 'dropped',
+ or the 'complete' set of samples is used.

This only has an effect when used with `annotation_variables`.

```{code-cell} ipython3
counts, label_ds = gsf.get_gem_data(agem, annotation_mask='complete', annotation_variables=['Treatment'])
```

### Transforming the Count Matrix

Often times a simple transform is required when using a count matrix, e.g. log transforms.

```{code-cell} ipython3

```
