# This code includes contributions from the sphinx-book-theme
# Original Repository: https://github.com/executablebooks/sphinx-book-theme
# Original License: BSD 3-clause
# https://github.com/executablebooks/sphinx-book-theme/blob/master/LICENSE

"""A lightweight book theme based on the pydata sphinx theme."""

# pylint: disable=unused-argument

import hashlib
import os
from functools import lru_cache
from pathlib import Path

from docutils import nodes as docutil_nodes
from docutils.parsers.rst.directives.body import Sidebar
from sphinx.application import Sphinx
from sphinx.locale import get_translation
from sphinx.util import logging

from .header_buttons import (
    add_header_buttons,
    prep_header_buttons,
    update_context_with_repository_info,
    update_sourcename,
)
from .header_buttons.launch import add_launch_buttons
from .header_buttons.utils import get_theme_options_dict

SPHINX_LOGGER = logging.getLogger(__name__)


def add_metadata_to_page(app, pagename, templatename, context, doctree):
    """Adds some metadata about the page that we re-use later."""
    # Add the site title to our context so it can be inserted into the navbar
    if not context.get("root_doc"):
        # TODO: Sphinx renamed master to root in 4.x, deprecate when we drop 3.x
        context["root_doc"] = context.get("master_doc")
    context["root_title"] = app.env.titles[context["root_doc"]].astext()

    # Update the page title because HTML makes it into the page title occasionally
    if pagename in app.env.titles:
        title = app.env.titles[pagename]
        context["pagetitle"] = title.astext()

    # Add a shortened page text to the context using the sections text
    if doctree:
        description = ""
        for section in doctree.traverse(docutil_nodes.section):
            description += section.astext().replace("\n", " ")
        description = description[:160]
        context["page_description"] = description

    # Add the author if it exists
    if app.config.author != "unknown":
        context["author"] = app.config.author


@lru_cache(maxsize=None)
def _gen_hash(path: str) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


def update_mode_thebe_config(app):
    """Update thebe configuration with SBT-specific values"""
    theme_options = get_theme_options_dict(app)
    if theme_options.get("launch_buttons", {}).get("thebe") is True:
        # In case somebody specifies they want thebe in a launch button
        # but has not activated the sphinx_thebe extension.
        if not hasattr(app.env.config, "thebe_config"):
            SPHINX_LOGGER.warning(
                "Thebe is activated but not added to extensions list. "
                "Add `sphinx_thebe` to your site's extensions list."
            )
            return
        # Will be empty if it doesn't exist
        thebe_config = app.env.config.thebe_config
    else:
        return

    if not theme_options.get("launch_buttons", {}).get("thebe"):
        return

    # Update the repository branch and URL
    # Assume that if there's already a thebe_config, then we don't want to over-ride
    if "repository_url" not in thebe_config:
        thebe_config["repository_url"] = theme_options.get("repository_url")
    if "repository_branch" not in thebe_config:
        branch = theme_options.get("repository_branch")
        if not branch:
            # Explicitly check in case branch is ""
            branch = "master"
        thebe_config["repository_branch"] = branch

    app.env.config.thebe_config = thebe_config


class Margin(Sidebar):
    """Goes in the margin to the right of the page."""

    optional_arguments = 1
    required_arguments = 0

    def run(self):
        """Run the directive."""
        if not self.arguments:
            self.arguments = [""]
        nodes = super().run()
        nodes[0].attributes["classes"].append("margin")

        # Remove the "title" node if it is empty
        if not self.arguments:
            nodes[0].children.pop(0)
        return nodes


def update_templates(app, pagename, templatename, context, doctree):
    """Update template names and assets for page build.

    This is a copy of what the pydata theme does here to include a new section
    - https://github.com/pydata/pydata-sphinx-theme/blob/0a4894fab49befc59eb497811949a1d0ede626eb/src/pydata_sphinx_theme/__init__.py#L173 # noqa: E501
    """
    # Allow for more flexibility in template names
    template_sections = ["theme_footer_content_items"]
    for section in template_sections:
        if context.get(section):
            # Break apart `,` separated strings so we can use , in the defaults
            if isinstance(context.get(section), str):
                context[section] = [
                    ii.strip() for ii in context.get(section).split(",")
                ]

            # Add `.html` to templates with no suffix
            for ii, template in enumerate(context.get(section)):
                if not os.path.splitext(template)[1]:
                    context[section][ii] = template + ".html"


def setup(app: Sphinx):
    # Events
    app.connect("builder-inited", update_mode_thebe_config)
    app.connect("builder-inited", update_sourcename)
    app.connect("builder-inited", update_context_with_repository_info)
    app.connect("html-page-context", add_metadata_to_page)
    app.connect("html-page-context", update_templates)

    # Header buttons
    app.connect("html-page-context", prep_header_buttons)
    # Bump priority so that it runs after the pydata theme sets up the edit URL func.
    app.connect("html-page-context", add_launch_buttons, priority=501)
    app.connect("html-page-context", add_header_buttons, priority=501)

    # Directives
    app.add_directive("margin", Margin)

    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
