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

# Selecting and Comparing Genes with GSForge


Our goal in examining RNA-seq data sets often reduces to "feature selection" -- to borrow a term from machine learning.
Examining the count data should give us an idea of which genes correlate (by a given measure) with a phenotype of intrest.
We can then examine that 'selected' set more closely and form biological or chemical hypothesis that explain the expression-phenotype link. Unfortunately there is not a gold standard measure to preform this selection, and practicing researchers must juggle a number of different methods and results.


---

**Enter `GSForge`**: 
A tool that helps collate and compare gene 'selection' resuls from a given method.
Here we present a brief analysis of a *oryza sativa* cultivar set[1], with the explicit purpose of demonstrating `GSForge`.

---


## Contents

1. Create an Annotated Gene Expression Matrix
    + The `AnnotatedGEM` object stores gene expression data alongside sample annotations.
    + [optional] save normalizations or transforms.
2. Select genes / features / tags
    + Using an R script: `edgeR`.
    + Using the Boruta algorithm with a random forest model from `sklearn`.
    + From the literature.
3. Compare results
    + Set operations on selection indices.
    + Compare and rank between sets.
    + Compare and rank within sets.
4. Next Steps
    + Create a comparative specifiaction.
    + Visualize results.

---

```{code-cell} ipython3
from os import  environ
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import patsy
import GSForge as gsf
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import quantile_transform

import matplotlib.pyplot as plt
# import seaborn as sns
import holoviews as hv
hv.extension('matplotlib')


# Declare paths.
OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage")).expanduser()
RAW_COUNT_PATH = OSF_PATH.joinpath("GEMmaker_GEMs", "Osativa_heat_drought_PRJNA301554.GEM.raw.txt")
HYDRO_LABEL_PATH = OSF_PATH.joinpath("raw_annotation_data", "PRJNA301554.hydroponic.annotations.txt")
SI_FILE_1_PATH = OSF_PATH.joinpath('GEMmaker_GEMs', 'raw_annotation_data', 'TPC2016-00158-LSBR2_Supplemental_File_1.csv')
```

## 1. Create an Annotated Gene Expression Matrix

The `AnnotatedGEM` contains our expression matrix and any annotations that can be sample or gene indexed.
This object can be created directly from `GEMmaker` output and (properly formatted) annotation matrixes.
See [this notebook](../walkthroughs/oryza_sativa/01-AnnotatedGEM_from_pandas) for detail.

Our [`GEMmaker` workflow](https://gemmaker.readthedocs.io/en/latest/) for the alignment and quantification presented here.

```{code-cell} ipython3
agem = gsf.AnnotatedGEM.from_files(
    count_path=RAW_COUNT_PATH,
    label_path=HYDRO_LABEL_PATH,
    # These are the default arguments passed to from_files,
    # to the individual calls to `pandas.read_csv`.
    count_kwargs=dict(index_col=0, sep="\t"),
    label_kwargs=dict(index_col=1, sep="\t"),
)
agem
```

Under the hood this is a light-weight wrapper for the `xarray.DataSet` object, which does all of the heavy lifting.

```{code-cell} ipython3
agem.data
```

### Select Counts and Annotations using `get_gem_data()`

The `AnnotatedGEM` object (and the `GeneSetCollection`, introduced further down) can have data subsets pulled from them easily using the `get_gem_data()` interface.

```{code-cell} ipython3
counts, labels = gsf.get_gem_data(agem, annotation_variables=['treatment', 'time', 'genotype'])
labels = labels.to_dataframe()
# Get time as an inetger.
labels['time'] = labels['time'].str.split(' ', expand=True).iloc[:, 0].astype(int)
# Infer missing replicate column.
rep_dfs = []
for g, df in labels.groupby(['treatment', 'time', 'genotype']):
    df['replicate'] = np.arange(1, df.shape[0] + 1)
    rep_dfs.append(df)

labels = pd.concat(rep_dfs)
labels.head()
```

```{code-cell} ipython3
counts.shape
```

```{code-cell} ipython3
import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
%load_ext rpy2.ipython
pandas2ri.activate()
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Supresses verbose R output.
```

```{code-cell} ipython3
%%R
library("edgeR")
```

### Add a Normalized or Transformed Count Matrix

The `AnnotatedGEM` object can hold more than one count matrix, so long as they share the same gene and sample coordinates.
Here we demonstrate adding a TPM normalized matrix as produced by `edgeR`.
This is more usefull for transforms that are computationally expensive, or that require data not easily stored in the `AnnotatedGEM` object.

We can then access a given count matrix by passing `count_variable='NAME'` to `get_gem_data()`.

```{code-cell} ipython3
agem.data['qt_counts'] = xr.DataArray(
    quantile_transform(counts.values, axis=1), 
    coords=counts.coords, 
    name='qt_counts')
agem.data['qt_counts']
```

## 2. Select genes / features / tags

While we can use the formula interface in R, it is often desirable to re-use the same design matrix and contrasts accross feature selection methods.
We use the [`patsy`](https://patsy.readthedocs.io/en/latest/) package to accomplish this.


---

**Note**:  
I end up creating multiple `GeneSetCollection` objects here.
I could have just as easily used one. 
How many collection objects used depends on your preferences and the complexity of your analysis.

---

### Create Design Matrixes and Constrasts

Here I make some simple contrasts that select for treatment from each of our designs. 
Newer users may find it easier to use the [`makeContrasts` function](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/makeContrasts) provided by the popular `limma` package.

```{code-cell} ipython3
ri_counts = gsf.utils.R_interface.Py_counts_to_R(counts)

forumla_designs = [
    "0 + treatment",
    "0 + treatment:genotype",
    # These models take too long to run for a demo-notebook.
    # "0 + treatment:genotype:scale(time, center=False)",
    # "0 + treatment:genotype:C(time):C(replicate)",
]

# I use nested lists instead of dictionaries to facilitate transfering to and from R.
design_list = [[f, patsy.dmatrix(f, labels)] for f in forumla_designs]
dmatrix_list = [[f, np.asarray(v)] for f, v in design_list]


def make_contrast(design, factors):
    return np.sum([(np.core.defchararray.find(design.design_info.column_names, f)!=-1).astype(int)
                   for f in factors], axis=0)


treatment_factors = design_list[0][1].design_info.column_names[1:]
control_factor = design_list[0][1].design_info.column_names[0]

contrast_list = list()

for i, (formula, design) in enumerate(design_list):
    design_contrast = list()
    for j, col_name in enumerate(treatment_factors):
        contrast = make_contrast(design, [col_name]) - make_contrast(design, [control_factor])
        design_contrast.append([col_name, contrast])
    contrast_list.append(design_contrast)        
```

Let's view the simplest model:

+++

### Select via `edgeR`

Here we present a 'standard' feature selection with `edgeR`. 

```{code-cell} ipython3
:tags: []

%%R -i ri_counts -i dmatrix_list -i contrast_list -o results -o keep

results <- list()

di <- 1

print('Applying filter.')
# Apply default gene filter.
keep <- filterByExpr(ri_counts, design=dmatrix_list[[1]][[2]])
Y <- ri_counts[keep,]

for (design in dmatrix_list)
{
    name <- design[[1]]
    design_matrix <- design[[2]]
    print(sprintf('Examing design %s', name))
    
    print('Construct the edgeR DGEList object.')
    # Construct the edge list object.
    dge_list <- DGEList(counts=Y)
    dge_list <- calcNormFactors(dge_list)
    dge_list <- estimateDisp(dge_list, design=design_matrix)

    # Plot variation.
    dir.create('disp_plots', showWarnings = FALSE)
    png(sprintf("disp_plots/%s_counts_disp.png", name))
    plotBCV(dge_list)
    dev.off()

    # Fit the negative binomial model.
    fit <- glmQLFit(dge_list, design=design_matrix)
    
    contrast_results <- list()
    ci <- 1
    
    for (contrast in contrast_list[[di]])
    {
        cname <- contrast[[1]]
        contrast_array <- contrast[[2]]
        print(sprintf('Examing contrast %s', cname))
        contrast_result  = glmQLFTest(fit, contrast=contrast_array)
        p_values <- contrast_result$table$PValue
        FDR <- p.adjust(p_values, method="BH")
        contrast_result$table$FDR <- FDR
        contrast_result$table$support_dir <- c(decideTests(contrast_result, p.value=0.001, lfc=2.0))
        contrast_result$table$support <- abs(c(decideTests(contrast_result, p.value=0.001, lfc=2.0)))
        contrast_results[[ ci ]] = contrast_result$table
        ci <- ci + 1
    }
    
    results[[di]] <- contrast_results
    di <- di + 1
}
```

```{code-cell} ipython3
dge_collection = gsf.GeneSetCollection(gem=agem, name='DGE Results')

for result_dfs, formula, contrasts in  zip(results, forumla_designs, contrast_list):
    for contrast_df, (contrast_name, _) in zip(list(result_dfs), contrasts):
        key = f"'{formula}'__{contrast_name}"
        dge_collection[key] = gsf.GeneSet(contrast_df, name=key)

dge_collection
```

```{code-cell} ipython3
# dge_collection["'0 + treatment'__treatment[DROUGHT]"].get_n_top_genes(threshold=2.0, mode="above_absolute_threshold")
```

### Random Forest via Boruta Feature Selection

Boruta is
Random forests..
sklearn...
...

```{code-cell} ipython3
boruta_gsc = gsf.GeneSetCollection(gem=agem, name='Boruta Results')

selecting_model = RandomForestClassifier(
    class_weight='balanced',
    max_depth=4, 
    n_estimators=1000, 
    n_jobs=-1)

for target in ["treatment", "genotype"]:
    boruta_treatment_ds = gsf.operations.BorutaProspector(
        agem,
        estimator=selecting_model,
        annotation_variables=target,
        max_iter=25)
    
    boruta_gsc[f"Boruta_{target}"] = gsf.GeneSet(boruta_treatment_ds, name=f"Boruta_{target}")
    
boruta_gsc
```

```{code-cell} ipython3
boruta_gsc
```

### Gene Sets from Literature

In practice there is often a gene set that is externally defined that we wish to compare.
Fortunately the only requirement for a GeneSet to be useful is that has an index that matches that in our GEM.
See [this notebook](../walkthroughs/oryza_sativa/04-GeneSets_from_Literature) for more details.

Here we read in a supporting information file from the EGRIN study by Wilkins *et. al.*, and after a bit of parsing it is converted to a series of `GeneSet` objects.

```{code-cell} ipython3
si1_df = pd.read_csv(SI_FILE_1_PATH, skiprows=3, index_col=0)
si1_df.head()
```

```{code-cell} ipython3
mappings = {'ChrSy.fgenesh.gene.37': 'ChrSy.fgenesh.mRNA.37'}


def parse_gene_splices(gene, gene_index: pd.Series, replacement_mappings=mappings):
    """Convert to an existing splice site, if possible."""
    if gene in gene_index:
        return gene
    
    for splice_site in range(1, 4):
        gene_splice = f'{gene}.{splice_site}'
        if gene_splice in gene_index:
            return gene_splice
        
    if replacement_mappings is not None:
        if replacement_mappings.get(gene):
            return replacement_mappings.get(gene)

        
lit_dge_coll = gsf.GeneSetCollection(gem=agem, name='Literature DGE')


for col in si1_df.columns:
    genes = si1_df[si1_df[col] > 0].index.values
    genes = np.asarray([parse_gene_splices(gene, agem.data.Gene.to_series())
                        for gene in genes])
    
    diff = np.setdiff1d(genes, agem.data.Gene.values)
    if diff.shape[0] > 0:
        print(diff)
    
    lit_dge_coll[col] = gsf.GeneSet.from_gene_array(genes, name=col)
    
lit_dge_coll.summarize_gene_sets()
```

## 3. Compare Results

While this step is ultimately experiment-specific, it usually entails some combination of:

1. Comparing selected to unselected features by one or more measure.
    + p-values and log-fold-change.
2. Ranking or comparing within a selected set.
    + Random forest feature importance.
    + p-values and log-fold-change.
3. Comparing set memberships, these usually take the form of set operations, such as:
    + union
    + intersection
    + difference
    + unique
4. Comparing between selection sets.
    + Model prediction scores.
    
    
Entries 1 and 2 are typically routine within the same selection set, as the tool used to create that set should provide the needed measures.
More complications emerge with steps 3 and 4, especially when comparing selection sets that derive from different methods.
Without additional wet lab experimentation we are limited in declaring our success in feature identification.
Instead we can score our selected subset as it preforms in prediction models.

For this demonstration we will combine each collection into its own set by taking the union of their selections.
Then we will examine the features selected for all treatments by each method.

```{code-cell} ipython3
lit_geneset = gsf.GeneSet.from_GeneSets(*lit_dge_coll.gene_sets.values(), name='literature_union')
lit_geneset
```

```{code-cell} ipython3
union_coll = gsf.GeneSetCollection(gem=agem, name='Combnied Collection')
union_coll.gene_sets.update(boruta_gsc.gene_sets)
union_coll.gene_sets.update({'literature_union': lit_geneset})

# Here I get the keys for each of the two DGE models we ran above.
dge_keys = pd.Series(dge_collection.gene_sets.keys())
key_sets = [dge_keys[dge_keys.str.startswith(f"'{f}'")].values for f in forumla_designs]

for f, keys in zip(forumla_designs, key_sets):
    name = f"combined '{f}'"
    union_coll[name] = gsf.GeneSet.from_GeneSets(
        *[dge_collection[k] for k in keys],
        name=name)
    
union_coll
```

### Visualize Set Overlap

Instead of a Venn diagram we use an 'Upset plot'.
This allows us to view overlaps of sets larger than three.

```{code-cell} ipython3
gsf.plots.collections.UpsetPlotInterface(union_coll)
```

## 4. Next Steps

```{code-cell} ipython3
umap_interface = gsf.panels.UMAP_Interface(union_coll, random_state=42)
layout = umap_interface.view(hue='treatment') + umap_interface.view(hue='genotype')
layout.opts(fig_size=200)
```

```{code-cell} ipython3
umap_interface = gsf.panels.UMAP_Interface(union_coll, random_state=42, selected_gene_sets=['all'])
layout = umap_interface.view(hue='treatment') + umap_interface.view(hue='genotype')
layout.opts(fig_size=200)
```

## References

1. Wilkins, O. et al. EGRINs (Environmental gene regulatory influence networks) in rice that function in the response to water deficit, high temperature, and agricultural environments. Plant Cell 28, 2365–2384 (2016).

2. LIMMA

3. UMAP

```{code-cell} ipython3

```
