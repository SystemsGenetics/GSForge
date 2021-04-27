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

# Creating DESeq2 GeneSets

This notebook covers how to run and load a basic `DESeq2` DEG result as a `GSForge.GeneSet`.

***Import Python packages***

```{code-cell} ipython3
from os import fspath, environ
from pathlib import Path
import GSForge as gsf
```

***R integration setup***

```{code-cell} ipython3
import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
%load_ext rpy2.ipython
pandas2ri.activate()
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
```

***Import R Packages***

```{code-cell} ipython3
%%R
library("DESeq2")
```

```{code-cell} ipython3
OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage")).expanduser()
HYDRO_GEM_PATH = OSF_PATH.joinpath("AnnotatedGEMs", "oryza_sativa_hydro_raw.nc")
DEG_COLL_PATH = OSF_PATH.joinpath("Collections", "DEG_gene_sets")
```

***Loading an AnnotatedGEM***

```{code-cell} ipython3
agem = gsf.AnnotatedGEM(HYDRO_GEM_PATH)
agem
```

```{code-cell} ipython3
ldf = agem.data[['treatment' ,'time']].to_dataframe()
ldf['time'] = ldf['time'].str.split(' ', expand=True).iloc[:, 0].astype(int)
ldf.head()
```

```{code-cell} ipython3

```

## Prepare input data for DESeq2

This requires us to drop genes that have counts of zero.

```{code-cell} ipython3
dropped_counts, labels = gsf.get_gem_data(agem, count_mask="dropped", annotation_variables=["treatment"])
```

These counts were made with Kallisto, so we must round them for use in `DEseq2`.

```{code-cell} ipython3
dropped_counts.shape
```

*Round counts to intergers*

Since our example GEM was aligned with Kallisto, which can give non-integer counts.

```{code-cell} ipython3
ri_dropped_counts = gsf.utils.R_interface.Py_counts_to_R(dropped_counts)
ri_dropped_counts = ri_dropped_counts.round()
ri_labels = labels.to_dataframe()
```

```{code-cell} ipython3
ri_dropped_counts.head(2)
```

```{code-cell} ipython3
ri_labels.head(2)
```

## Run `DESeq2`

```{code-cell} ipython3
%%R -i ri_dropped_counts -i ri_labels -o deseq_df

dds <- DESeqDataSetFromMatrix(countData = ri_dropped_counts,
                              colData = ri_labels,
                              design= ~ treatment)
dds <- DESeq(dds)
deseq_results <- results(dds)
deseq_df = data.frame(deseq_results)
```

We can examine the $\beta$ matrix produced by `DESeq2` directly:

```{code-cell} ipython3
%%R

resultsNames(dds)
```

```{code-cell} ipython3
deseq_df.head()
```

And save the entire dataset as a `GeneSet`.

```{code-cell} ipython3
deseq2_treatment = gsf.GeneSet(deseq_df, 
                               name="deseq2_treatment", 
                               attrs={"DESeq2_formula": "~ Treatment"})
deseq2_treatment
```

```{code-cell} ipython3
deseq2_treatment.data
```

This lets us keep all the data produced by the analysis, while setting the `support` to whatever we deem this method selected.

```{code-cell} ipython3
def top_n_abs(dataframe, n=100, padj_cuttoff=0.05):
    """Returns the top n most (absolutely) differentially expressed genes from a deseq2 result.
    This also filters by p-values."""
    filtered_df = dataframe[dataframe["padj"] < padj_cuttoff]
    filtered_df = filtered_df.reindex(filtered_df["log2FoldChange"].abs().sort_values().index)
    return filtered_df.tail(n).index
```

```{code-cell} ipython3
top_n_abs(deseq_df, n=5)
```

```{code-cell} ipython3
deseq2_treatment = deseq2_treatment.set_support_by_genes(top_n_abs(deseq_df))
deseq2_treatment
```

```{code-cell} ipython3
deseq2_treatment.save_as_netcdf(DEG_COLL_PATH)
```
