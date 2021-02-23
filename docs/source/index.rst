===================================
Welcome to GSForge's documentation!
===================================

.. note::
    This is a pre-release version of ``GSForge``, and some features and function names may not yet be stable.

``GSForge`` provides the following utilities to aid in collating and analyzing feature selection results
from expression count matrices:

1. Data structures
2. Membership set selection operators
3. Plotting functions

It aims to provide tools that are species-agnostic, that allow a practicing researcher to more easily collate and
compare selection results across a variety of methods (from both R and Python).
``GSForge`` is not an analysis tool itself, but the documentation provides example implementations.


Overview
========

You should consider using ``GSForge`` when you have:

+ More than one feature selection method.
+ Many parameter sets to compare using one (or more) feature selection methods.
+ More than one normalization to explore.
+ You have a feature selection method that is non-deterministic.

You should seriously consider using ``GSForge`` if you find yourself in more than one of the above categories.


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

2. Set Selection Operators
--------------------------

todo...

3. Plotting Functions
---------------------

todo...

Examples
========

Examples are available...

Installation
============

.. include:: user_guide/installation


.. toctree::
    :titlesonly:
    :maxdepth: 1

    Welcome <self>
    User Guide <user_guide/user_guide>
    Reference Examples <reference_examples/ref_examples_index>
    Walk-through Examples <walkthroughs/walkthroughs_index>
    API <API/GSForge>
    About <about>
