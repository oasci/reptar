# MIT License
# 
# Copyright (c) 2022, Alex M. Maldonado
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

import os
import numpy as np
from .parsers import parserORCA, parserXTB
from .utils import get_md5

def identify_parser(out_path):
    """Identifies the correct parser depending on some trigger.
    Only supported packaged should be included in `triggers`.

    Parameters
    ----------
    out_path : :obj:`str`
        Path to output file.
    
    Returns
    -------
    :obj:`reptar.creator`
        One of the supported creator classes.
    """
    with open(out_path, 'r') as f:
        for line in f:
            for parser, phrases, do_break in triggers:
                if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
                    filetype = parser
                    if do_break:
                        return filetype

###   Runtime Functions   ###

# Triggers to identify output files.
triggers = [
    (parserORCA, ["O   R   C   A"], True),
    (parserXTB, ["x T B"], True)
]

class creator:
    """Handles creating reptar data files.

    Parameters
    ----------
    data : ``reptar.data``, optional
        An initialized data object.
    """

    def __init__(self, data=None):
        if data is not None:
            self.data = data

    def parse_output(
        self, out_path, geom_path=None, traj_path=None, extractors=None
    ):
        """Parse output file using cclib and custom parser.

        Parameters
        ----------
        out_path : :obj:`str`
            Path to computational chemistry output file.
        """
        self.out_path = os.path.abspath(out_path)
        if geom_path is not None:
            geom_path = os.path.abspath(geom_path)
        if traj_path is not None:
            traj_path = os.path.abspath(traj_path)
        
        self.geom_path = geom_path
        self.traj_path = traj_path

        packageParser = identify_parser(self.out_path)
        self.parser = packageParser(
            self.out_path, geom_path=geom_path, traj_path=traj_path,
            extractors=extractors
        )
        self.parsed_info = self.parser.parse()
    
    def create_group(self, group_key):
        """Create group using parsed information.
        
        Parameters
        ----------
        group_key : :obj:`str`
            Key to the desired new group (including parent).
        """
        assert hasattr(self, 'data')
        parsed_info = self.parsed_info

        self.data.init_group(group_key)
        # Loop through each category of data.
        for cat_key in parsed_info.keys():
            for data_key in parsed_info[cat_key].keys():
                data = parsed_info[cat_key][data_key]
                self.data.add(f'{group_key}/{data_key}', data)
        
        md5 = get_md5(self.data.get(group_key))
        self.data.add(f'{group_key}/md5', md5)
        
        try:
            md5_group = get_md5(self.data.get(group_key), only_dataset=True)
            self.data.add(f'{group_key}/md5_data', md5_group)
        except Exception:
            pass

        return self.data.get(group_key)

