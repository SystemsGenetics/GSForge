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
release = 'v0.6beta'
version = '0.6'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.napoleon",
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    "sphinx_rtd_theme",
    # 'sphinx_gallery.load_style',
    'sphinx.ext.mathjax',
    "myst_nb",
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = [
#     'in_dev/*',
# ]


# -- sphinx autodoc config ------------------------------------------------
autodoc_member_order = 'bysource'
autosummary_generate = True

# -- MyST configuration ---------------------------------------------------
myst_enable_extensions = [
    "amsmath",
    "deflist",
    "dollarmath",
]

jupyter_execute_notebooks = "off"
execution_timeout = 1500


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['static']

html_css_files = [
    'custom.css',
]
