***********
Development
***********

Overview
--------

*GSForge* is maintained via travis ci integration.
Upon pushing a tagged commit (that matches the regex string in travis.yml file) to the master branch,
travis will build and update the documents, then update the pip, docker and conda repositories.

The same deployment flow can be triggered under development, and occurs when either an
tag that matches the development regex, or when 'travis_dev' is in a commit message.
This development run deploys documentation to https://systemsgenetics.github.io/GSForgeDev/index.html.


PyPi Repository
---------------




-------
