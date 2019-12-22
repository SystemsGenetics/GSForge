**********
User Guide
**********

*Under Construction...*


Getting Started
---------------

There are three core models that hold data in GSForge.

The AnnotatedGEM holds the gene expression matrix (GEM) as well as any gene or sample annotations.

* AnnotatedGEM
  A data class for a gene expression matrix and any associated sample or gene annotations.
* GeneSet
  A data class for a the result of a gene selection or analysis.
* GeneSetCollection
  A data class that holds an AnnotatedGEM and a dictionary of associated GeneSet objects.

So in using GSForge to help analyze your GEM:

+ Load the count matrix and any annotations into an AnnotatedGEM.
+ Perform one or more analyses, or gene-selections, and store each as a GeneSet.
+ Load the created AnnotatedGEM and GeneSet(s) into a GeneSetCollection.
+ Compare the selected GeneSets.


Creating and working with an AnnotatedGEM
-----------------------------------------

* `AnnotatedGEM from pandas <AnnotatedGEM_from_pandas.ipynb>`_
  How to load your count and label .csv or .txt files into an `xarray.Dataset` using `GSForge`, and save them as a single netcdf (.nc) file for future use.

* `GEM Normalization <GEM_Normalization.ipynb>`_
  How to run and save different normalizations to the same `AnnotatedGEM` object.


Selecting Genes & Creating a GeneSet
------------------------------------

* `Selecting Genes with Boruta <Selecting_Genes_with_Boruta.ipynb>`_
  How to select genes based on a sample label using random forests via the Boruta algorithm.


Creating a GeneSetCollection and Analyzing Results
--------------------------------------------------

* `GeneSet Analysis <GeneSet_Analysis.ipynb>`_
  How to run a basic gene (feature) selection using random forests and the Boruta all-relevant feature selection algorithm.


Other Guides
------------

* `Plotting Guide <plotting_guide/index.rst>`_
  Plotting examples by plot type.

* `Workflow Guide <workflow_guide/index.rst>`_
  Using nextflow to run many different boruta analyses.

* `R-Integration Guide <R_integration/index.rst>`_
  Integrating results from R packages into GeneSet objects.



.. toctree::
    :titlesonly:
    :maxdepth: 2

    AnnotatedGEM from pandas <AnnotatedGEM_from_pandas>
    GEM Normalization <GEM_Normalization>
    Selecting Genes with Boruta <Selecting_Genes_with_Boruta>
    GeneSet Analysis <GeneSet_Analysis>
    Plotting Guide <plotting_guide/index>
    Workflow Guide <workflow_guide/index>
    R-Integration Guide <R_integration/index>
