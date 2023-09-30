# This code includes contributions from the pydata-sphinx-theme
# Original Repository: https://github.com/pydata/pydata-sphinx-theme
# Original License: BSD 3-clause
# https://github.com/pydata/pydata-sphinx-theme/blob/main/LICENSE

"""General helpers for the management of config parameters."""

from typing import Any, Dict

from collections.abc import Iterator

from docutils.nodes import Node
from sphinx.application import Sphinx


def get_theme_options_dict(app: Sphinx) -> dict[str, Any]:
    """Return theme options for the application w/ a fallback if they don't exist.

    The "top-level" mapping (the one we should usually check first, and modify
    if desired) is ``app.builder.theme_options``. It is created by Sphinx as a
    copy of ``app.config.html_theme_options`` (containing user-configs from
    their ``conf.py``); sometimes that copy never occurs though which is why we
    check both.
    """
    if hasattr(app.builder, "theme_options"):
        return app.builder.theme_options
    elif hasattr(app.config, "html_theme_options"):
        return app.config.html_theme_options
    else:
        return {}


def config_provided_by_user(app: Sphinx, key: str) -> bool:
    """Check if the user has manually provided the config."""
    # pylint: disable-next=protected-access
    return any(key in ii for ii in [app.config.overrides, app.config._raw_config])


def traverse_or_findall(node: Node, condition: str, **kwargs) -> Iterator[Node]:
    """Triage node.traverse (docutils <0.18.1) vs node.findall.

    TODO: This check can be removed when the minimum supported docutils version
    for numpydoc is docutils>=0.18.1.
    """
    return (
        node.findall(condition, **kwargs)
        if hasattr(node, "findall")
        else node.traverse(condition, **kwargs)
    )
