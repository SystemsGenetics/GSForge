# -*- coding: utf-8 -*-
# flake8: noqa (hacky way of sharing config, etc...)
# See here: https://github.com/pyviz-dev/nbsite/blob/master/nbsite/shared_conf.py
# For values imported from shared_conf
import os
import sys

from nbsite.shared_conf import *

sys.path.insert(0, os.path.abspath('..'))

extensions += [
    "nbsphinx",
    "sphinx.ext.napoleon",  # numpy-style docstrings.
    'nbsite.gallery'
]

nbsite_gallery_conf = {
    'galleries': {
        # This key must be the same as a folder name within /examples.
        'plot_gallery': {
            'title': 'Plot Gallery',
            'intro': 'Demonstrations of plotting functions provided by GSForge.',
            'orphans': ['Plotting_Guide.ipynb'],
            # The orphans key allows the user to pass a list of files that will be rendered to html without
            # being thumbnailed and linked from the gallery page. The main use-case for this is when a section
            # has an index which provides an overview of the section and directs users through
            # the notebooks in a particular order.
        },
        "panel_gallery": {
            'title': 'Panel Gallery',
            'intro': 'Demonstrations of "panel-ized" applications provided by GSForge.',
        },
        "how_to_galleries": {
            'title': 'How-to Guides',
            'intro': 'How-to guides provided by GSForge.',
            'sections': [
                {'path': 'core_guides',
                 'title': 'Core How-to Guides'},
                {
                    'path': 'R_integration_guides',
                    'title': 'R Integration Guides',
                    # 'skip': [
                    #     'DESeq2_GeneSets.ipynb',
                    #     'EdgeR_GeneSets.ipynb',
                    #     'R_GEM_normalizations.ipynb'
                    # ]
                },
                {'path': 'workflow_guide',
                 'title': 'Workflow Integration How-to Guides'},

            ],
        },
    },
    'github_org': 'SystemsGenetics',
    'github_project': 'GSForge',
    'thumbnail_url': 'https://github.com/SystemsGenetics/GSForgeDev/tree/gh-pages/assets/thumbnails',
}

project = u'GSForge'
authors = u'Tyler Biggs'
copyright = u'2019 - 2020 ' + authors
description = 'GSForge is a Python software package that assists researchers in the selection of ' \
              'gene sets with potential association to an experimental condition or phenotypic trait, ' \
              'which offers new potential hypotheses for gene-trait causality.'

version = '0.5'
release = 'alpha'

html_static_path += ['_static']
html_theme = 'sphinx_ioam_theme'
# logo file etc should be in html_static_path, e.g. assets
html_theme_options = {
    #    'custom_css':'bettercolors.css',
    #    'logo':'amazinglogo.png',
    #    'favicon':'amazingfavicon.ico'
}

_NAV = (
    ('Welcome', 'Welcome'),
    ('User Guide', 'user_guide/index'),
    ('API', 'Reference_Manual/GSForge'),
    ('Developer Guide', 'Development'),
    ('About', 'About')
)

html_context.update({
    'PROJECT': project,
    'DESCRIPTION': description,
    'AUTHOR': authors,
    'VERSION': version,
    'NAV': _NAV,
    # by default, footer links are same as those in header
    'LINKS': _NAV,
    'SOCIAL': (
        ('Github', 'https://github.com/SystemsGenetics/GSForge'),
    )
})
