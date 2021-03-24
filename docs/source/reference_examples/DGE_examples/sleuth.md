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

# Creating `sleuth` GeneSets

This notebook covers how to run and load a basic `sleuth` DEG result as a `GSForge.GeneSet`.

[`sleuth` documentation](https://github.com/pachterlab/sleuth).

```{code-cell} ipython3
from pathlib import Path
import numpy as np
import pandas as pd
import GSForge as gsf
```

***R integration setup***

```{code-cell} ipython3
from rpy2.robjects import pandas2ri
%load_ext rpy2.ipython
pandas2ri.activate()
```

***Import R Packages***

```{code-cell} ipython3
%%R
suppressMessages({
  library("sleuth")
})
```

```{code-cell} ipython3
label_df = pd.read_csv('Maturity_Annotations.csv', sep='\t')
label_df.head()
```

```{code-cell} ipython3
label_df.shape
```

```{code-cell} ipython3
count_df = pd.read_csv('/scidas/John/08-Apple_Maturity_Project/02-GEMmaker/08-output/GEMs/Kamiak.GEM.raw.txt', sep='\t')
count_df
```

```{code-cell} ipython3
count_df.shape
```

```{code-cell} ipython3
count_df.columns.values[-5:]
```

```{code-cell} ipython3
vals = pd.Series(count_df.columns).str.split('_S', expand=True).iloc[:, 0].values
vals[:5]
```

```{code-cell} ipython3
match = np.intersect1d(vals, label_df['RNASample.ID'].str.replace("-", "_").values)
match.shape
```

```{code-cell} ipython3
diff = np.setdiff1d(vals, label_df['RNASample.ID'].str.replace("-", "_").values)
diff.shape
```

```{code-cell} ipython3
diff
```

```{code-cell} ipython3
np.setdiff1d(label_df['RNASample.ID'].str.replace("-", "_").values, vals)
```

```{code-cell} ipython3
label_df.shape
```

```{code-cell} ipython3
label_df.columns
```

Get the pcolumnss to kallisto output:

```{code-cell} ipython3
# kallisto_files = list(Path("/scidas/John/08-Apple_Maturity_Project/02-GEMmaker/08-output").glob("*"))
# kallisto_files = [x for x in kallisto_files
#                   if x.name not in ["reports", "GEMs"]]
# kallisto_files[:5]
```

```{code-cell} ipython3

```

```{code-cell} ipython3

```
