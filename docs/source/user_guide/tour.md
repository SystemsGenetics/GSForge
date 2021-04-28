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
    + Save normalizations or transforms.
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
from sklearn import model_selection
from sklearn import linear_model
import umap
import umap.plot

import matplotlib.pyplot as plt
import colorcet as cc
import seaborn as sns
import holoviews as hv
hv.extension('matplotlib')


# Declare paths.
OSF_PATH = Path(environ.get("GSFORGE_DEMO_DATA", default="~/GSForge_demo_data/osfstorage")).expanduser()
RAW_COUNT_PATH = OSF_PATH.joinpath("GEMmaker_GEMs", "Osativa_heat_drought_PRJNA301554.GEM.raw.txt")
HYDRO_LABEL_PATH = OSF_PATH.joinpath("raw_annotation_data", "PRJNA301554.hydroponic.annotations.txt")
SI_FILE_1_PATH = OSF_PATH.joinpath('GEMmaker_GEMs', 'raw_annotation_data', 'TPC2016-00158-LSBR2_Supplemental_File_1.csv')
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

### Add a Normalized or Transformed Count Matrix

The `AnnotatedGEM` object can hold more than one count matrix, so long as they share the same gene and sample coordinates.
Here we demonstrate adding a TPM normalized matrix as produced by `edgeR`.
This is more usefull for transforms that are computationally expensive, or that require data not easily stored in the `AnnotatedGEM` object.

We can then access a given count matrix by passing `count_variable='NAME'` to `get_gem_data()`.

```{code-cell} ipython3
counts, _ = gsf.get_gem_data(agem)
agem.data['qt_counts'] = xr.DataArray(
    quantile_transform(counts.values, axis=1), 
    coords=counts.coords, 
    name='qt_counts')
# agem.data['qt_counts']
```

### Select Counts and Annotations using `get_gem_data()`

The `AnnotatedGEM` object (and the `GeneSetCollection`, introduced further down) can have data subsets pulled from them easily using the `get_gem_data()` interface.

+++

UMAP DESCRIPTION.

```{code-cell} ipython3
counts, labels = gsf.get_gem_data(agem, annotation_variables=['treatment', 'genotype'])
mapper = umap.UMAP(densmap=True).fit(counts.values)
fig, axes = plt.subplots(1, 2, figsize=(20, 10))
umap.plot.points(mapper, labels=labels['treatment'], background='black', ax=axes[0], color_key_cmap='Set1');
umap.plot.points(mapper, labels=labels['genotype'], background='black', ax=axes[1], color_key_cmap='Set2');
```

```{code-cell} ipython3
counts, labels = gsf.get_gem_data(agem, annotation_variables=['treatment', 'genotype'], count_variable='qt_counts')
mapper = umap.UMAP(densmap=True).fit(counts.values)
fig, axes = plt.subplots(1, 2, figsize=(20, 10))
umap.plot.points(mapper, labels=labels['treatment'], background='black', ax=axes[0], color_key_cmap='Set1');
umap.plot.points(mapper, labels=labels['genotype'], background='black', ax=axes[1], color_key_cmap='Set2');
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

Here we present a 'standard' feature selection with `edgeR`, see the [edgeR Bioconductor](https://bioconductor.org/packages/release/bioc/html/edgeR.html) page for more information.

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

### Random Forest via Boruta Feature Selection

See this [FAQ](https://notabug.org/mbq/Boruta/wiki/FAQ) for details on Boruta.

In brief, boruta lets us use a 'minimal optimal' feature selection model as an 'all relevant' method.
Random forest models typically work by finding a minimum optimal set of features.
In our use, we are interested in all features that may be relevant, not the minimum required to train a model.
`GSForge` provides a warpper which returns `xarray.DataSet` objects suitable for immediate `GeneSet` construction.

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
        max_iter=100)
    
    boruta_gsc[f"Boruta_{target}"] = gsf.GeneSet(boruta_treatment_ds, name=f"Boruta_{target}")
    
boruta_gsc
```

We can check if boruta still has features to resolve into important or not by examining the 'support_weak' variable.
In practice one should increase the number of iterations until there are no remaining values to be resolved.

```{code-cell} ipython3
boruta_gsc['Boruta_genotype'].data.support_weak.values.sum()
```

### Gene Sets from Literature

In practice there is often a gene set that is externally defined that we wish to compare.
Fortunately the only requirement for a GeneSet to be useful is that has an index that matches that in our GEM.
See [this notebook](../walkthroughs/oryza_sativa/04-GeneSets_from_Literature) for more details.

Here we read in a supporting information file from the EGRIN study by Wilkins *et. al.*, and after a bit of parsing it is converted to a series of `GeneSet` objects.
This SI table appears to be counts of significant differential expression for time contrasts within each treatment.
The full model is not (exactly) specified in the paper, but we can infer that it was probably:

$$
~ 0 + genotype:treatment:C(time):C(replicate)
$$

Where $C(\text{label})$ indicates treatment as a categorical, rather than scalar variable.
The union for each set of treatment contrasts was probably then used.

```{code-cell} ipython3
si1_df = pd.read_csv(SI_FILE_1_PATH, skiprows=3, index_col=0)
# Manual index fix.
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
union_coll = gsf.GeneSetCollection(gem=agem, name='Combnied Collection')
union_coll.gene_sets.update(boruta_gsc.gene_sets)
lit_geneset = gsf.GeneSet.from_GeneSets(*lit_dge_coll.gene_sets.values(), name='literature_union')
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

See the [upsetplot documentation](https://upsetplot.readthedocs.io/en/stable/) for more details.

Instead of a Venn diagram we use an 'Upset plot'.
This allows us to view overlaps of sets larger than three.

```{code-cell} ipython3
gsf.plots.collections.UpsetPlotInterface(union_coll)
```

### Comparing Selection Sets

We can estimate how well a given subset of genes 'describes' a sample (phenotype) label by comparing how well they perform using a given machine learning model.

```{code-cell} ipython3
:tags: []

results = dict()

for key in list(union_coll.gene_sets.keys()) + ['all']:
    counts, treatment = gsf.get_gem_data(union_coll, selected_gene_sets=[key], annotation_variables=["treatment"])
    x_train, x_test, y_train, y_test = model_selection.train_test_split(counts, treatment)
    # model = RandomForestClassifier(class_weight='balanced', n_estimators=1000, n_jobs=-1, max_depth=6)
    model = linear_model.Perceptron()
    model.fit(x_train, y_train) 
    results[key] = model.score(x_test, y_test)
    

hv.Bars(results, kdims=["Gene Selection Group"]).opts(
    xrotation=90, invert_axes=True, ylim=(0, 1.1), ylabel='Score', fig_size=150, 
    aspect=2, title='Scores vs Treatment Labels')
```

### Ranking Features within Sets

GSForge provides helper functions to extract genes by a score values.
Note we get the data from the original dge collection, as those logFC values are intact.

```{code-cell} ipython3
dge_ds = dge_collection["'0 + treatment:genotype'__treatment[HEAT]"]
# dge_ds.data
```

```{code-cell} ipython3
dge_ds.get_top_n_genes("logFC", 10)
```

```{code-cell} ipython3
dge_ds.get_genes_by_threshold(3.0, "logFC")
```

### Rank Genes with a Random Forest

Random forests and feature ranks.
Robust enough to function in our case.
Some values filtered prior to dge analysis...

```{code-cell} ipython3
union_coll["'0 + treatment:genotype'__treatment[HEAT]"] = dge_collection["'0 + treatment:genotype'__treatment[HEAT]"]

gene_rank_mdl = RandomForestClassifier(class_weight='balanced', n_estimators=1000, n_jobs=-2)

treatment_nFDR = gsf.operations.nFDR(
    union_coll,
    selected_gene_sets=["Boruta_treatment", "'0 + treatment:genotype'__treatment[HEAT]"],
    gene_set_mode="union",
    annotation_variables=["treatment"],
    model=gene_rank_mdl,
    n_iterations=5
)

treatment_feature_importance = gsf.operations.RankGenesByModel(
    union_coll,
    selected_gene_sets=["Boruta_treatment", "'0 + treatment:genotype'__treatment[HEAT]"],
    gene_set_mode="union",
    annotation_variables=["treatment"],
    model=gene_rank_mdl,
    n_iterations=5
)
```

```{code-cell} ipython3
gene_union = union_coll.union(["Boruta_treatment", "'0 + treatment:genotype'__treatment[HEAT]"])

boruta_support = np.isin(gene_union, union_coll["Boruta_treatment"].gene_support())
dge_support = np.isin(gene_union, 
                      union_coll["'0 + treatment:genotype'__treatment[HEAT]"].gene_support())

support = np.zeros_like(gene_union)
support[boruta_support] = "Boruta"
support[dge_support] = "DGE"
support[(boruta_support * dge_support) == True] = "Both"

df = pd.DataFrame({
    "logFC": dge_ds.data['logFC'].reindex(Gene=gene_union).values,
    "Gene": gene_union,
    "mean feature importance": treatment_feature_importance.feature_importance_mean.values,
    "mean nFDR": treatment_nFDR.nFDR_mean.values,
    "support source": support,
}).set_index("Gene")


# sns.pairplot(df, hue="support source", markers='.', diag_kind="kde", 
#              plot_kws=dict(edgecolor=None, alpha=0.25));

fig, axes = plt.subplots(1, 2, figsize=(20, 10))
sns.scatterplot(data=df, x='mean feature importance', y='mean nFDR', hue='support source',
                edgecolor=None, alpha=0.5, marker='.', ax=axes[0]);
sns.scatterplot(data=df, x='mean feature importance', y='logFC', hue='support source',
                edgecolor=None, alpha=0.5, marker='.', ax=axes[1]);
```

### UMAP Embeddings of Selections

```{code-cell} ipython3
counts, labels = gsf.get_gem_data(union_coll, selected_gene_sets=['Boruta_treatment'], count_variable='qt_counts', 
                                  annotation_variables=['treatment', 'genotype'])
mapper = umap.UMAP(densmap=True).fit(counts.values)
fig, axes = plt.subplots(1, 2, figsize=(20, 10))
umap.plot.points(mapper, labels=labels['treatment'], background='black', ax=axes[0], color_key_cmap='Set1');
umap.plot.points(mapper, labels=labels['genotype'], background='black', ax=axes[1], color_key_cmap='Set2');
```

```{code-cell} ipython3
counts, labels = gsf.get_gem_data(union_coll, selected_gene_sets=['literature_union'], count_variable='qt_counts', 
                                  annotation_variables=['treatment', 'genotype'])
mapper = umap.UMAP(densmap=True).fit(counts.values)
fig, axes = plt.subplots(1, 2, figsize=(20, 10))
umap.plot.points(mapper, labels=labels['treatment'], background='black', ax=axes[0], color_key_cmap='Set1');
umap.plot.points(mapper, labels=labels['genotype'], background='black', ax=axes[1], color_key_cmap='Set2');
```

### Clustermap Selection

```{code-cell} ipython3
def series_to_colors(series, cmap, categorical=True):
    keys = series.unique()
    colors = hv.plotting.util.process_cmap(cmap, len(keys), categorical=categorical)
    mapper = {k: c for k, c in zip(keys, colors)}
    return series.map(mapper)

# az_samples = agem.data.sel(Sample=agem.data.genotype == 'Azuenca (AZ; IRGC#328, Japonica)')['Sample'].values
df, labels = gsf.get_gem_data(union_coll, 
#                               sample_subset=az_samples, 
                              annotation_variables=['treatment', 'genotype'],
                              selected_gene_sets=["Boruta_treatment", 
                                                  "literature_union"],
                              gene_set_mode='intersection',
                              count_transform=lambda counts: np.log2(counts.where(counts > 0)),
                              output_type="pandas")


color_df = pd.DataFrame({
    "treatment": series_to_colors(labels['treatment'], "Set1"),
    "genotype": series_to_colors(labels['genotype'], "Set2")
})

sns.clustermap(df.fillna(0), 
               metric='cityblock', 
               row_colors=color_df,
               row_cluster=False, 
               dendrogram_ratio=0.1,
               cmap='jet',
               figsize=(10, 10));
```

## 4. Conclusion & Next Steps

Links to other notebooks and resources.

+++

## References

1. Wilkins, O. et al. EGRINs (Environmental gene regulatory influence networks) in rice that function in the response to water deficit, high temperature, and agricultural environments. Plant Cell 28, 2365–2384 (2016).

2. Ritchie, M. E. et al. Limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 43, e47 (2015).

3. McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. (2018).

4. Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R. & Pfister, H. UpSet: Visualization of Intersecting Sets Europe PMC Funders Group. IEEE Trans Vis Comput Graph 20, 1983–1992 (2014).
