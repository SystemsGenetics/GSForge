===================================
Welcome to GSForge's documentation!
===================================

.. note::
    This is a pre-release version of ``GSForge``, and some features and function names may not yet be stable.

``GSForge`` provides the following utilities to aid in collating and analyzing feature selection results from
expression count matrices:

1. Data structures that index expression and annotation in a single file.
2. A common interface for membership set selection operators.
3. Plotting functions for comparing expression as well as set intersections.

It aims to provide tools that are species-agnostic, that allow a practicing researcher to more easily collate and
compare selection results across a variety of methods (from both R and Python). ``GSForge`` is not an analysis tool
itself, but the documentation provides example implementations.

.. image:: ./static/rna_seq_generation_highlights.svg


Installation
============


``GSForge`` works with Python 3.7+.
It is currently only tested on Linux.


Install via pip
---------------

::

    pip install GSForge

Install via Github source
-------------------------

Clone the repository::

    git clone git@github.com:SystemsGenetics/GSForge.git

Then install locally via pip::

    pip install ./GSForge


Overview
========

You should consider using ``GSForge`` when you have:

* More than one feature selection method.
* Many parameter sets to compare using one (or more) feature selection methods.
* More than one normalization to explore.
* You have a feature selection method that is non-deterministic.

You should seriously consider using ``GSForge`` if you find yourself in more than one of the above categories.
The general workflow is then:

1. Import Data and create an :ref:`AnnotatedGEM<user_guide/gem_creation>`.
2. Find or select genes of interest.
3. Create a GeneSetCollection.
4. Compare GeneSets.


1. Data Structures
------------------

The core structural components ``GSForge`` provides are:

**AnnotatedGEM**
A data object that contains:
gene expression matrix with named coordinates for genes and samples.
Sample annotations / descriptors.
Can store more than one transform of the GEM (e.g. raw counts and TPM normalized counts).

**GeneSet**
A list of gene names, as well as any desired metadata or numerical analysis.
Examples of data or analysis that fit within this model include:
differential gene expression,
machine learning feature selections,
and literature sets.

**GeneSetCollection**
An interface that provides utility access to an AnnotatedGEM based on provided GeneSets.
e.g. GEM subset selection by GeneSet membership, or by some set operation (union, intersection, difference)
from multiple sets.


2. Set Selection Interface
--------------------------

In comparing feature selection results the most basic needs are often set operations.
e.g., the intersection, union and unique members between one or more sets.
Often time subsets of the original data are examined based on a set-operation of one or more feature sets.

``GSForge`` provides this functionality through the ``GSForge.Interface`` object,
and the ``GSForge.get_gem_data`` function.


3. Plotting Functions
---------------------

``GSForge`` provides a number of plots common to gene expression analysis.
The ``Interface`` object powers these plotting functions, which allows users to easily plot only subsets of data based
on set membership.


.. toctree::
   :titlesonly:
   :maxdepth: 1
   :caption: User Guide / Tutorial:

    Installation <user_guide/installation>
    Tour <user_guide/tour>
    Annotated Expression Matrix Creation <user_guide/gem_creation>
    Set and Collection Creation <user_guide/set_and_collection_creation>
    Interface Guide <user_guide/Interface_Guide>
    Plotting Guide <user_guide/Plotting_Guide>
    R Integration Guide <user_guide/R_integration>
    Analysis Resources <user_guide/analysis_resources>


.. toctree::
   :maxdepth: 2
   :caption: Reference Examples:

    Reference Plots <reference_examples/plot_gallery/ref_plots>
    Normalization Examples <reference_examples/Normalizations/norm_index>


.. toctree::
   :maxdepth: 2
   :caption: Walkthroughs:

    Oryza sativa <walkthroughs/oryza_sativa/oryza_sativa_index>


.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: API Reference:

    API <API/modules>
    Development <development>
