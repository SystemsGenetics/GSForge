# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# Add gsforge directory to the path.
sys.path.insert(0, os.path.abspath('../../GSForge/'))

# Add param docstring tools to the path.
sys.path.insert(0, os.path.abspath('.'))
# Import param docstring tools. You can ignore linting warnings re this import.
from paramdoc import param_formatter, param_skip


# -- param docstring setup ---------------------------------------------------
# Docstring formatting for classes that inherit from param.Parameterized.
#
def setup(app):
    app.connect("autodoc-process-docstring", param_formatter)
    app.connect("autodoc-process-docstring", param_skip)


# -- Project information -----------------------------------------------------

project = 'gsforge'
copyright = '2021, Tyler Biggs, Stephen Ficklin'
author = 'Tyler Biggs, Stephen Ficklin'

# The full version, including alpha/beta/rc tags
release = 'v0.5beta2'
version = '0.52'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.napoleon",
    'sphinx.ext.autodoc',
    # 'recommonmark',
    "sphinx_rtd_theme",
    # 'nbsphinx',
    # 'sphinx_gallery.load_style',
    'sphinx.ext.mathjax',
    "myst_nb",
    # 'sphinx_copybutton',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = [
# ]

# -- nbsphinx configuration ---------------------------------------------------

# nbsphinx_allow_errors = True
# nbsphinx_kernel_name = 'gsfenv'
# nbsphinx_execute = 'always'
# jupyter_execute_notebooks = "force"
# execution_excludepatterns = ['R_integration_guides', 'walkthroughs', 'panel_gallery', 'workflow_guide']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# _NAV = (
#     ('Welcome', 'index'),
#     ('User Guide', 'user_guide/user_guide'),
#     ('Reference Examples', 'reference_examples/ref_examples_index'),
#     ('Walk-throughs', 'walkthroughs/walkthroughs_index'),
#     ('API', 'API/GSForge'),
#     # ('FAQ', 'FAQ'),
#     ('About', 'about')
# )
