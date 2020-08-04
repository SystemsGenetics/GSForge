***********
Development
***********


Travis CI Overview
------------------

*GSForge* is maintained via travis ci integration.
Upon pushing a tagged commit (that matches the regex string in travis.yml file) to the master branch,
travis will build and update the documents, then update the pip, docker and conda repositories.

The same deployment flow can be triggered under development, and occurs when either an
tag that matches the development regex, or when 'travis_dev' is in a commit message.
This development run deploys documentation to https://systemsgenetics.github.io/GSForgeDev/index.html.

*Developer 'release' tag regex and examples:*

- ``^v(\d+|\.)+[a-z]+\d*$``
- ``v1.2rc4``

*Release tag examples:*

^v(\d+|\.)+[^a-z]\d+$

git tag -a v1.4 -m "my version 1.4"

Tags are not pushed by default. See documentation on `git tags
<https://git-scm.com/book/en/v2/Git-Basics-Tagging>`_. for more information.


*Commit message activation of travis development pipeline:*

git commit -m 'some message, travis_dev'


Travis integration handles updates to all of the other repositories below.


PyPi Repository
---------------

Under construction...


Conda Repository
----------------

Under construction...


DockerHub Repository
--------------------

Under construction...


Documentation Website
---------------------

Under construction...


-------
