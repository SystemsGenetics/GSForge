===========
Development
===========

.. note::
    Under construction.

Notes and instructions on updating and contributing to GSForge.


Documentation
=============

A number of tools are required to build the package components and documentation:

* `sphinx <https://domain.invalid/>`
* `jupyter <https://domain.invalid/>`
* `nbconvert <https://domain.invalid/>`
* `paramdoc <https://domain.invalid/>`
* `MyST <https://domain.invalid/>`
* jupyter book?
* jupyter-cache
* `Github pages <https://domain.invalid/>`
* `Docker <https://domain.invalid/>`
* `DockerHub <https://domain.invalid/>`
* `Travis Continuous Integration <https://domain.invalid/>`

Most of these components are managed using travis continuous integration. The exception being those demonstrations of
tools external to ``GSForge``, such as the tour notebook, workflow executions, and the walkthroughs, which must have
their results committed to the git repository.


Travis CI Integration
=====================

``GSForge`` is maintained via travis ci integration. A build will initialize upon:
travis will build and update the documents, then update the pip, docker and conda repositories.

* Upon pushing a tagged commit (that matches the regex string in travis.yml file) to the master branch
* a commit tag that matches the development regex, or when ``travis_dev`` is in a commit message.

This development run deploys documentation to  https://systemsgenetics.github.io/GSForgeDev/index.html.


Notebooks
=========

The example jupyter notebooks are prepared using the ``MyST`` markdown format, which allows the notebooks
to be saved as plain text, separate from their output.

``https://myst-parser.readthedocs.io/en/latest/index.html``

They can be linked to normal ``.ipynb`` files for use via jupyter notebooks or jupyter lab.

https://jupytext.readthedocs.io/en/latest/

Building Documentation
======================

There are two stages to building the documentation.

    1. Build the API with ``sphinx-apidoc``.
    2. Build the website and run the notebooks with ``sphinx``.

Most reference examples are executed during the travis continuous integration build process, while the walkthroughs
and some of the notebooks take too long to run, and must be executed locally. This is enforced by saving files that
should be run every time as markdown files ``.md`` to be processed as notebooks automatically by MyST. All jupyter
notebooks should be run prior to updating the main branch.

https://myst-parser.readthedocs.io/en/latest/index.html

..
    jupyter nbconvert --to notebook --execute --inplace docs/source/user_guide/tour.ipynb
    jupyter nbconvert --to notebook --execute --inplace docs/source/walkthroughs/oryza_sativa/*.ipynb

    jupytext notebook.ipynb --to myst


We can run the travis process locally, see here: https://stackoverflow.com/a/49019950

..
    BUILDID="build-$RANDOM"
    INSTANCE="travisci/ci-sardonyx"

    docker run --name "build-local" -dit "travisci/ci-sardonyx" /sbin/init
    docker exec -it $BUILDID bash -l


From the top directory of the repository::

    # Construct the API documentation.
    sphinx-apidoc --separate --force -o docs/source/API/ ./GSForge
    # Build the documentation.
    python -m sphinx docs/source/ ../gsforge_docs/


Convert an existing jupyter notebook to an all-text ``MyST`` format::

    jupytext my_file.ipynb --to myst

If a notebook should run everytime travis integration is called, ensure the checkpoint files are not commited.
For notebooks that should not be run (many take too long for travis), ensure the checkpoint file is commited.


OSF Data Repository
===================


