=============================
Documentation and Development
=============================


Building the Docs
=================

The recommended way to build the docs is through a docker image, which is managed with ``docker-compose``::

    # optional, if you built previously:
    docker-compose rm -f
    # Create and run the documentation service:
    docker-compose -f docker-compose.yml up --build documentation
    # build the documentation service
    docker-compose -f docker-compose.yml build documentation
    # Create the documentation.
    docker-compose run --rm documentation

This runs the ``CMD`` of the ``systemsgenetics/gsforge_documentation`` Docker file which:

    1. Builds the API with ``sphinx-apidoc``.
    2. Builds the website and run the notebooks with ``sphinx``.

Most reference examples are executed during the travis continuous integration build process, while the walkthroughs
and some of the notebooks take too long to run, and must be executed locally. This is enforced by saving files that
should be run every time as markdown files ``.md`` to be processed as notebooks automatically by MyST. All jupyter
notebooks should be run prior to updating the main branch.

Notebooks can be executed from the command line via::

    # Execute a single notebook.
    jupyter nbconvert --to notebook --execute --inplace docs/source/user_guide/tour.ipynb
    # Execute multiple notebooks.
    jupyter nbconvert --to notebook --execute --inplace docs/source/walkthroughs/oryza_sativa/*.ipynb


Resources
---------

A number of tools are required to build the package components and documentation:

* `sphinx <https://domain.invalid/>`
* `jupyter <https://domain.invalid/>`
* `nbconvert <https://domain.invalid/>`
* `paramdoc <https://domain.invalid/>`
* `MyST <https://domain.invalid/>`
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


OSF Data Repository
===================


