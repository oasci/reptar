#!/usr/bin/env python3

# # MIT License
#
# Copyright (c) 2022-2023, OASCI
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Convert YAML definitions into rst files"""

import os
import shutil

import yaml

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))


# Collect paths of all YAML files
def get_files(path, expression, recursive=True):
    r"""Returns paths to all files in a given directory that matches a provided
    expression in the file name.

    Commonly used to find all files of a certain type, e.g., output or xyz
    files.

    Parameters
    ----------
    path : :obj:`str`
        Specifies the directory to search.
    expression : :obj:`str`
        Expression to be tested against all file names in ``path``.
    recursive : :obj:`bool`, default: ``True``
        Recursively find all files in all subdirectories.

    Returns
    -------
    :obj:`list` [:obj:`str`]
        All absolute paths to files matching the provided expression.
    """
    # pylint: disable=redefined-outer-name
    if path[-1] != "/":
        path += "/"
    if recursive:
        all_files = []
        for dirpath, _, filenames in os.walk(path):
            index = 0
            while index < len(filenames):
                if dirpath[-1] != "/":
                    dirpath += "/"
                filenames[index] = dirpath + filenames[index]
                index += 1
            all_files.extend(filenames)
        files = []
        for f in all_files:
            if expression in f:
                files.append(f)
    else:
        files = []
        for f in os.listdir(path):
            filename = os.path.basename(f)
            if expression in filename:
                files.append(path + f)
    return files


# Retrieve file paths.
print("Finding YAML files")
yaml_file_paths = get_files("../reptar/definitions", ".yaml")
yaml_file_paths = [os.path.abspath(path) for path in yaml_file_paths]

print("Loading YAML files")
# Load all YAML files.
yaml_files = {}
for file_path in yaml_file_paths:
    file_name = os.path.splitext(os.path.basename(file_path))[0]

    with open(file_path, encoding="utf-8") as f:
        yaml_files[file_name] = yaml.safe_load(f)


WRITE_DIR = os.path.abspath("./source/definitions")

# Remove directory if it exists
if os.path.exists(WRITE_DIR):
    shutil.rmtree(WRITE_DIR)
os.makedirs(WRITE_DIR)

# Write definitions main rst file
print("Writing definitions.rst")
print(f"Directory: {WRITE_DIR}")
defs_main_path = os.path.join(WRITE_DIR, "definitions.rst")

DEFS_MAIN_HEADING = r"""===========
Definitions
===========

One challenging aspect of repositories is communicating the meaning of each label we give our data.
Reptar tries to suggest common data types encountered in the computational chemistry community.
Labels and their respective descriptions are in YAML files that could be stored directly in the data.
We consider three types of data:

- **identification:** unique labels that describe the within.
- **system_info:** specifies the atomistic system.
- **runtime_info:** results and information from calculations.

We further separate definitions into the following categories.

.. toctree::
    :maxdepth: 1

"""

with open(defs_main_path, mode="w", encoding="utf-8") as f:
    f.write(DEFS_MAIN_HEADING)
    labels = list(yaml_files.keys())
    labels.sort()
    for label in labels:
        f.write(f"    {label}\n")


# Write individual files
print("Writing definitions files")
for file_name, defs in yaml_files.items():
    file_path = os.path.join(WRITE_DIR, file_name + ".rst")
    header_len = int(len(file_name))
    header_string = "=" * header_len

    with open(file_path, mode="w", encoding="utf-8") as f:
        f.write(header_string + "\n")
        f.write(file_name + "\n")
        f.write(header_string + "\n")

        for defs_cat_name, defs_dicts in defs.items():
            cat_len = int(len(defs_cat_name))
            cat_string = "-" * cat_len

            f.write("\n" + defs_cat_name + "\n")
            f.write(cat_string + "\n\n")

            f.write(".. glossary::\n")
            for def_name, def_dict in defs_dicts.items():
                f.write("\n" + " " * 4 + def_name + "\n")
                f.write(" " * 8 + def_dict["description"] + "\n")
                f.write(" " * 8 + "**Type:** :obj:`" + def_dict["type"] + "`" + "\n")
