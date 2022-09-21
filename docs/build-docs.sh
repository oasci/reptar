#!/usr/bin/env bash
# Build sphinx documentation with multiversion.
# Will not build current branch documentation; only main and tagged commits.
cd "${0%/*}"
rm -rf ./html/
sphinx-multiversion -nT ./source/ ./html/
touch ./html/.nojekyll