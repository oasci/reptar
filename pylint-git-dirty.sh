#!/usr/bin/env bash
# Runs pylint only on dirty git files.
pylint --rcfile=.pylintrc `git diff --name-only --diff-filter=d | grep -E '\.py$' | tr '\n' ' '`
