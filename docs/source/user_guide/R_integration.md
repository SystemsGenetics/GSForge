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
```{code-cell}
from pathlib import Path
import numpy as np
import pandas as pd
import GSForge as gsf
```


```{code-cell}
from rpy2.robjects import pandas2ri
%load_ext rpy2.ipython
pandas2ri.activate()
```

# R Integration

```{code-cell}
#%%R
#suppressMessages({
#  library("sleuth")
#})
```

```{code-cell}
#%%R
#counts, labels = gsf.get_gem_data(agem, annotation_variables=['treatment', 'time', 'genotype'])
#labels = labels.to_dataframe()
#ri_counts = gsf.utils.R_interface.Py_counts_to_R(counts)
```


```{code-cell}
#%%R -i ri_counts -i dmatrix_list -i contrast_list -o results -o keep
```
