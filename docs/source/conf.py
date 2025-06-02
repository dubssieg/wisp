# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

print("SYS PATH", sys.path)
project = 'MicroTaxo'
copyright = '2025, PNRIA'
author = 'PNRIA'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'nbsphinx',  # pour les notebooks
    'myst_parser', #pour le readme.md
]

nbsphinx_execute = 'never'  # ou 'auto' si tu veux exécuter les notebooks

myst_enable_extensions = ['colon_fence', 'deflist'] # Exemples d'extensions myst que vous pouvez activer

templates_path = ['_templates']
exclude_patterns = []

language = 'fr'
# source_suffix = ['.rst', '.md']
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static', '../../wisp_light/notebooks']
html_css_files = [
    'custom.css',
]