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
        'plot_gallery': {  # This key must be the same as a folder name within /examples
            'title': 'Plot Gallery',
            'intro': 'Demonstrations of plotting functions provided by GSForge.',
        },
        "panel_gallery": {
            'title': 'Panel Gallery',
            'intro': 'Demonstrations of "panel-ized" applications provided by GSForge.',
        },
        "user_guide": {  # This key must be the same as a folder name within /examples
            'title': 'User Guide',
            'intro': 'How-to guides provided by GSForge.',
            'orphans': ['overview.ipynb'],
            # The orphans key allows the user to pass a list of files that will be rendered to html without
            # being thumbnailed and linked from the gallery page. The main usecase for this is when a section
            # has an index which provides an overview of the section and directs users through
            # the notebooks in a particular order.
            'sections': [
                {
                    'path': 'R_integration_guide',
                    'title': 'R Integration Guide',
                    'skip': True,
                }
            ],
        },
    },
    'github_org': 'SystemsGenetics',
    'github_project': 'GSForge',
    # 'thumbnail_url': 'https://assets.holoviews.org/thumbnails',
}

project = u'GSForge'
authors = u'Tyler Biggs'
copyright = u'2019 - 2020' + authors
description = 'Short description for html meta description.'

version = '0.4'
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
    # ('Gallery', 'gallery/index'),
    ('API', 'Reference_Manual/index'),
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
