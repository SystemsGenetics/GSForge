===========
Development
===========

Notes and instructions on updating and contributing to GSForge.


Travis CI Integration
=====================

``GSForge`` is maintained via travis ci integration.
Upon pushing a tagged commit (that matches the regex string in travis.yml file) to the master branch,
travis will build and update the documents, then update the pip, docker and conda repositories.

The same deployment flow can be triggered under development, and occurs when either an
tag that matches the development regex, or when ``travis_dev`` is in a commit message.
This development run deploys documentation to https://systemsgenetics.github.io/GSForgeDev/index.html.


Notebooks
=========

The example jupyter notebooks are prepared using the ``MyST`` markdown format, which allows the notebooks
to be saved as plain text, separate from their output.

``https://myst-parser.readthedocs.io/en/latest/index.html``

They can be linked to normal ``.ipynb`` files for use via jupyter notebooks or jupyter lab.

https://jupytext.readthedocs.io/en/latest/

Building Documentation
======================

From the top directory of the repository::

    # Construct the API documentation.
    sphinx-apidoc --separate --force -o docs/source/API/ ./GSForge
    # Build the documentation.
    python -m sphinx docs/source/ ../gsforge_docs/


Convert an existing jupyter notebook to an all-text ``MyST`` format::

    jupytext my_file.ipynb --to myst
