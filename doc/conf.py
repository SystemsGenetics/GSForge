# -*- coding: utf-8 -*-
# flake8: noqa (hacky way of sharing config, etc...)

from nbsite.shared_conf import *
import os
import sys

sys.path.insert(0, os.path.abspath('../../gemprospector'))

###################################################
# edit things below as appropriate for your project

# extensions += ['nbsite.gallery']
#
# nbsite_gallery_conf = {
#     'backends': ['bokeh', 'matplotlib'],
#     'default_extensions': ['*.ipynb', '*.py'],
#     'enable_download': True,
#     'examples_dir': os.path.join('..', 'examples'),
#     'galleries': {
#         'user_guide': {'title': 'User Guide'}
#     },
#     'github_org': 'SystemsGenetics',
#     'github_project': 'GemProspector',
#     # 'thumbnail_url': 'https://assets.holoviews.org/thumbnails',
#     'within_subsection_order': lambda key: key
# }

project = u'GEMprospector'
authors = u'Tyler Biggs'
copyright = u'2019 ' + authors
description = 'Short description for html meta description.'

version = '0.0.1'
release = '0.0.1'

html_static_path += ['_static']
html_theme = 'sphinx_ioam_theme'
# logo file etc should be in html_static_path, e.g. _static
html_theme_options = {
    #    'custom_css':'bettercolors.css',
    #    'logo':'amazinglogo.png',
    #    'favicon':'amazingfavicon.ico'
}

_NAV = (
    # ('Getting Started', 'getting_started/index'),
    ('Feature Tour', 'Feature_Tour'),
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
    # will work without this - for canonical (so can ignore when building locally or test deploying)
    # 'WEBSITE_SERVER': 'https://ceball.github.io',
    'VERSION': version,
    'NAV': _NAV,
    # by default, footer links are same as those in header
    'LINKS': _NAV,
    # 'SOCIAL': (
    #     ('Gitter', '//gitter.im/ioam/holoviews'),
    #     ('Twitter', '//twitter.com/holoviews'),
    #     ('Github', '//github.com/ioam/holoviews'),
    # )
})
