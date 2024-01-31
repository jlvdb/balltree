# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------

import os
import sys

try:
    try:  # user has installed the package
        import balltree
    except ImportError:  # try local package location

        sys.path.insert(0, os.path.abspath("../.."))
        import balltree
except ImportError as e:
    if "core._math" in e.args[0]:
        raise RuntimeError("balltree must be compiled") from e

from numpy.typing import ArrayLike, NDArray

sys.path.insert(0, os.path.abspath("."))
from parse_c_doc import c_doc_to_rst

# -- Project information -----------------------------------------------------

project = "balltree"
copyright = "2024, Jan Luca van den Busch"
author = "Jan Luca van den Busch"
release = balltree.__version__
version = ".".join(release.split(".")[:2])


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

master_doc = "index"
extensions = [
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = []

autodoc_inherit_docstrings = True
autodoc_type_aliases = {
    ArrayLike: "ArrayArray",
    NDArray: "NDArray"
}
autodoc_member_order = "bysource"
autosummary_generate = True
autoclass_content = "both"

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
copybutton_only_copy_prompt_lines = True
copybutton_line_continuation_character = "\\"

# -- Options for HTML output -------------------------------------------------

pypi = "https://pypi.org/project/balltree"
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_favicon = "_static/icon.ico"
html_css_files = ["css/custom.css"]
# html_favicon = "_static/icon.ico"
html_theme_options = {
    "github_url": "https://github.com/jlvdb/balltree.git",
    "collapse_navigation": True,
    "navigation_depth": 3,
    "show_nav_level": 3,
    "show_toc_level": 3,
    "navbar_align": "content",
    "secondary_sidebar_items": ["page-toc"],
    "logo": {
        "image_light": "_static/logo.svg",
        "image_dark": "_static/logo.svg",
    },
    "pygment_light_style": "xcode",
    "pygment_dark_style": "github-dark",
}
html_sidebars = {
    "**": ["search-field.html", "sidebar-nav-bs.html", "sidebar-ethical-ads.html"]
}
html_context = {
    "default_mode": "auto",
}

# -- Build custom files ------------------------------------------------------

# generate the index page
index_text = """
Documentation of balltree
=========================

%README%

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
"""
toc = """
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   examples
   api
"""
with open("../../README.rst") as f:
    content = f.read()
    readme_text = "".join(content.split("|", 2)[2])  # drop header
with open("index.rst", "w") as f:
    f.write(index_text.replace("%README%", readme_text).replace(".. toc", toc))

# generate the BallTree documentation
with open("api/balltree.BallTree.rst", "w") as f:
    f.write("balltree.BallTree\n=================\n\n\n")
    f.write(c_doc_to_rst("../../balltree/balltree.c", "balltree"))
