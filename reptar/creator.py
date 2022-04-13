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

class creator():
    """"""

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
    
    def create_group(self, name, parent):
        """Create exdir group using parsed information.
        
        Parameters
        ----------
        name : :obj:`str`
            New group name.
        parent : ``obj``
            An exdir File or Group this new Group will be nested under.
        """
        parsed_info = self.parsed_info

        grp = parent.create_group(name)

        # Loop through each category of data.
        for cat_key in parsed_info.keys():
            cat_info = parsed_info[cat_key]
            for data_key in cat_info.keys():
                data = cat_info[data_key]
                
                # Custom handling of parsed data to ensure logical behavior.
                if data_key == 'dipole_moment':
                    if isinstance(data, np.ndarray):
                        data = data.tolist()
                    grp.attrs[data_key] = data
                    continue

                # If data is not array/list or array/list of one dimension and
                # length we make it an exdir attribute.
                # We make anything else an exdir dataset.
                if isinstance(data, np.ndarray) or isinstance(data, list) or isinstance(data, tuple):
                    # If data is a list we check the types of data it contains.
                    # If all of them are strings, we do not convert to array.
                    if not isinstance(data, np.ndarray):
                        # If there is only one item in the list we can likely
                        # just store it as an attribute. However, we have
                        # to be careful when the single item is an array.
                        if len(data) == 1:
                            data = data[0]
                            if isinstance(data, np.ndarray):
                                grp.create_dataset(data_key, data=data)
                            else:
                                grp.attrs[data_key] = data
                        # If all the items are not all strings then we make
                        # the dataset.
                        elif not all(isinstance(i, str) for i in data):
                            data = np.array(data)
                            grp.create_dataset(data_key, data=data)
                        # At this point only data that contains all strings
                        # should be left. We put these as attributes (i.e., 
                        # comp_ids)
                        else:
                            grp.attrs[data_key] = data
                    # Handling arrays.
                    else:
                        # If only one item we store it as an attribute.
                        # Otherwise we create the dataset.
                        if data.shape == (1,):
                            grp.attrs[data_key] = data[0].item()
                        else:
                            grp.create_dataset(data_key, data=data)
                # Handles all non iterable data. 
                else:
                    grp.attrs[data_key] = data
        
        grp.attrs['md5'] = get_md5(grp)
        try:
            grp.attrs['md5_data'] = get_md5(grp, only_dataset=True)
        except Exception:
            pass

        return grp

