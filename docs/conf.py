# Configuration file for the Sphinx documentation builder.
#
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = "sptxinsight"
copyright = "2026, Chao-Hui Huang"
author = "Chao-Hui Huang"

try:
    import sptxinsight

    release = getattr(sptxinsight, "__version__", "0.1.0")
except Exception:  # docs must build without the package installed
    release = "0.1.0"
version = release

# -- General configuration ---------------------------------------------------

# Only stdlib Sphinx extensions are used so the docs build in any environment
# without optional theme / autoapi packages.
extensions = [
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

language = "en"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "alabaster"
html_static_path = ["_static"]
