===================================
Welcome to GSForge's documentation!
===================================

.. note::
    This is a pre-release version of ``GSForge``, and some function names may not be stable.

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

``GSForge`` works with Python 3.7+. It is currently only tested on Linux.


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
2. Find or select genes of interest and store the results in a :ref:`Gene Set<user_guide/gem_creation>`.
3. Collate both the AnnotatedGEM and the sets of interest in a :ref:`GeneSetCollection<user_guide/set_and_collection_creation>`.
4. Use the :ref:`Interface<user_guide/Interface_Guide>` to access and compare sets.


.. toctree::
   :titlesonly:
   :maxdepth: 1
   :caption: User Guide / Tutorial:

    Introduction <self>
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
