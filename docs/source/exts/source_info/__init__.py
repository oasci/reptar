# This code includes contributions from the pydata-sphinx-theme
# Original Repository: https://github.com/pydata/pydata-sphinx-theme
# Original License: BSD 3-clause
# https://github.com/pydata/pydata-sphinx-theme/blob/main/LICENSE

"""Bootstrap-based sphinx theme from the PyData community."""

from typing import Dict

import json
import os
from pathlib import Path
from urllib.parse import urlparse

import requests

# pylint: disable-next=redefined-builtin
from requests.exceptions import ConnectionError, HTTPError, RetryError
from sphinx.application import Sphinx
from sphinx.errors import ExtensionError
from sphinx.util import logging

# pylint: disable-next=import-self
from . import edit_this_page, utils

logger = logging.getLogger(__name__)


def update_config(app):
    """Update config with new default values and handle deprecated keys."""
    # By the time `builder-inited` happens, `app.builder.theme_options` already exists.
    # At this point, modifying app.config.html_theme_options will NOT update the
    # page's HTML context (e.g. in jinja, `theme_keyword`).
    # To do this, you must manually modify `app.builder.theme_options`.
    theme_options = utils.get_theme_options_dict(app)

    # Validate icon links
    if not isinstance(theme_options.get("icon_links", []), list):
        raise ExtensionError(
            "`icon_links` must be a list of dictionaries, you provided "
            f"type {type(theme_options.get('icon_links'))}."
        )

    # Set the anchor link default to be # if the user hasn't provided their own
    if not utils.config_provided_by_user(app, "html_permalinks_icon"):
        app.config.html_permalinks_icon = "#"

    # check the validity of the theme switcher file
    is_dict = isinstance(theme_options.get("switcher"), dict)
    should_test = theme_options.get("check_switcher", True)
    if is_dict and should_test:
        theme_switcher = theme_options.get("switcher")

        # raise an error if one of these compulsory keys is missing
        json_url = theme_switcher["json_url"]
        theme_switcher["version_match"]  # pylint: disable=pointless-statement

        # try to read the json file. If it's a url we use request,
        # else we simply read the local file from the source directory
        # display a log warning if the file cannot be reached
        reading_error = None
        if urlparse(json_url).scheme in ["http", "https"]:
            try:
                request = requests.get(json_url, timeout=30)
                request.raise_for_status()
                content = request.text
            except (ConnectionError, HTTPError, RetryError) as e:
                reading_error = repr(e)
        else:
            try:
                # pylint: disable-next=unspecified-encoding
                content = Path(app.srcdir, json_url).read_text()
            except FileNotFoundError as e:
                reading_error = repr(e)

        if reading_error is not None:
            logger.warning(
                f'The version switcher "{json_url}" file cannot be read due to the following error:\n'
                f"{reading_error}"
            )
        else:
            # check that the json file is not ill formed,
            # throw a warning if the file is ill formed and an error if it's not json
            switcher_content = json.loads(content)
            missing_url = any(["url" not in e for e in switcher_content])
            missing_version = any(["version" not in e for e in switcher_content])
            if missing_url or missing_version:
                logger.warning(
                    f'The version switcher "{json_url}" file is malformed'
                    ' at least one of the items is missing the "url" or "version" key'
                )

    # Add an analytics ID to the site if provided
    analytics = theme_options.get("analytics", {})
    if analytics:
        # Plausible analytics
        plausible_domain = analytics.get("plausible_analytics_domain")
        plausible_url = analytics.get("plausible_analytics_url")

        # Ref: https://plausible.io/docs/plausible-script
        if plausible_domain and plausible_url:
            kwargs = {
                "loading_method": "defer",
                "data-domain": plausible_domain,
                "filename": plausible_url,
            }
            app.add_js_file(**kwargs)

        # Google Analytics
        gid = analytics.get("google_analytics_id")
        if gid:
            gid_js_path = f"https://www.googletagmanager.com/gtag/js?id={gid}"
            gid_script = f"""
                window.dataLayer = window.dataLayer || [];
                function gtag(){{ dataLayer.push(arguments); }}
                gtag('js', new Date());
                gtag('config', '{gid}');
            """

            # Link the JS files
            app.add_js_file(gid_js_path, loading_method="async")
            app.add_js_file(None, body=gid_script)

    # Update ABlog configuration default if present
    fa_provided = utils.config_provided_by_user(app, "fontawesome_included")
    if "ablog" in app.config.extensions and not fa_provided:
        app.config.fontawesome_included = True

    # Handle icon link shortcuts
    shortcuts = [
        ("twitter_url", "fa-brands fa-square-twitter", "Twitter"),
        ("bitbucket_url", "fa-brands fa-bitbucket", "Bitbucket"),
        ("gitlab_url", "fa-brands fa-square-gitlab", "GitLab"),
        ("github_url", "fa-brands fa-square-github", "GitHub"),
    ]
    # Add extra icon links entries if there were shortcuts present
    # TODO: Deprecate this at some point in the future?
    icon_links = theme_options.get("icon_links", [])
    for url, icon, name in shortcuts:
        if theme_options.get(url):
            # This defaults to an empty list so we can always insert
            icon_links.insert(
                0,
                {
                    "url": theme_options.get(url),
                    "icon": icon,
                    "name": name,
                    "type": "fontawesome",
                },
            )
    theme_options["icon_links"] = icon_links

    # Prepare the logo config dictionary
    theme_logo = theme_options.get("logo")
    if not theme_logo:
        # In case theme_logo is an empty string
        theme_logo = {}
    if not isinstance(theme_logo, dict):
        raise ValueError(f"Incorrect logo config type: {type(theme_logo)}")
    theme_options["logo"] = theme_logo


def setup(app: Sphinx) -> dict[str, str]:
    """Setup the Sphinx application."""

    # app.connect("builder-inited", update_config)
    app.connect("html-page-context", edit_this_page.setup_edit_url)

    return {"parallel_read_safe": True, "parallel_write_safe": True}
