# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os, sys
sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------
project = 'FP-Tools'
copyright = '2026, Los Alamos National Lab'
author = 'Leah Hartman, Alfred Farris, Jerome Daligault'
version = u'0.1.0'
release = u'0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'breathe',
    'sphinx_copybutton'
]
breathe_projects = {"fp-tools": "../xml"}
breathe_default_project = "FP-Tools"

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'breeze'
html_theme_options = {
    "header_tabs": False,
    "collapse_navigation": False,
    "sticky_navigation": True,      # Sidebar stays visible when scrolling
    "navigation_depth": -1,         # Show all levels (no limit)
    "includehidden": True,          # Include hidden TOCs
    "titles_only": False,           # Show full TOC, not just titles
}
html_static_path = ['_static']

# Tell Sphinx what the primary language being documented is
primary_domain = 'cpp'

# Tell Sphinx what the pygments highlight language should be
highlight_language = 'cpp'


