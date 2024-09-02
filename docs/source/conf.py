# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "sigmazerosearch"
copyright = "2024, Niam Patel"
author = "Niam Patel"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autosummary",
    "sphinx.ext.githubpages",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "autodoc2",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "smartquotes",
    "strikethrough",
]

myst_url_schemes = {
    "http": None,
    "https": None,
    "doi": "https://doi.org/{{path}}",
    "github": {
        "url": "https://github.com/{{path}}",
        "title": "{{path}}",
        "classes": ["github"],
    },
    "docdb": {
        "url": "https://microboone-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid={{path}}",
        "title": "DocDB-{{path}}",
    },
    "pypi": "https://pypi.org/project/{{path}}",
}

autodoc2_packages = ["../../sigmazerosearch"]
autodoc2_output_dir = "api"
autodoc2_render_plugin = "myst"
autodoc2_docstring_parser_regexes = [
    # this will render all docstrings as Markdown
    (r".*", "myst"),
]
autodoc2_class_docstring = "both"

templates_path = ["_templates"]
exclude_patterns = []

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "uproot": ("https://uproot.readthedocs.io/en/latest", None),
    "awkward-array": ("https://awkward-array.org/doc/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "hist": ("https://hist.readthedocs.io/en/latest", None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = "furo"
html_static_path = ["_static"]
html_title = "sigmazerosearch"

html_theme_options = {
    "dark_logo": "logo-dark.png",
    "light_logo": "logo.png",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/playonverbs/sigmazerosearch",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}