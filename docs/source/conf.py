# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------

import os

try:
    try:  # user has installed the package
        import balltree
    except ImportError:  # try local package location
        import sys

        sys.path.insert(0, os.path.abspath("../.."))
        import balltree
except ImportError as e:
    if "core._math" in e.args[0]:
        raise RuntimeError("balltree must be compiled") from e

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
    # "logo": {
    #     "image_light": "_static/logo-light.svg",
    #     "image_dark": "_static/logo-dark.svg",
    # },
    "pygment_light_style": "xcode",
    "pygment_dark_style": "github-dark",
}
html_sidebars = {
    "**": ["search-field.html", "sidebar-nav-bs.html", "sidebar-ethical-ads.html"]
}
html_context = {
    "default_mode": "auto",
}
