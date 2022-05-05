#!/usr/bin/env bash
cd "${0%/*}"
rm -rf ./source/doc/
rm -rf ./html/
sphinx-apidoc --force -o ./source/doc/ ../reptar 
# sphinx-build -nT ./source/ ./html/
sphinx-multiversion -nT ./source ./html
touch ./html/.nojekyll