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

from .parser import parser
from .orca_extractor import extractorORCA
import cclib

class parserORCA(parser):
    """Custom parser for ORCA output files.
    """

    def __init__(self, out_path, geom_path=None, traj_path=None, extractors=None):
        if extractors is None:
            extractors = []
        extractors.insert(0, extractorORCA())
        super().__init__(out_path, extractors)
        # TODO: Handle geom and traj paths
    
    def parse(self):
        """Parses trajectory file and extracts information.
        """
        # cclib parsed information.
        try:
            self.cclib_data = cclib.io.ccread(self.out_path)
            self.map_cclib_data()
        except Exception as e:
            raise e

        # Extract information.
        with open(self.out_path, mode='r') as f:
            for line in f:
                for extractor in self.extractors:
                    for i in range(len(extractor.triggers)):
                        if extractor.triggers[i][0](line):
                            getattr(extractor, extractor.triggers[i][1])(f, line)
                            break
                    else:
                        continue
                    break
        self.combine_extracted()
        self.after_parse()
        # TODO: Figure out how to get driver.

        return self.parsed_info

    def after_parse(self):
        """Checks to perform after parsing output file.
        """

        """
        if 'scf_grid_level_final' not in self.parsed_info['runtime_info']['keywords'].keys():
            self.parsed_info['runtime_info']['scf_grid_level_final'] = self.parsed_info['runtime_info']['scf_grid_level']
        
        if 'conv_val_geo_energy' in self.parsed_info['runtime_info']['keywords'].keys():
            self.parsed_info['runtime_info']['keywords']['conv_val_geo_energy'].insert(0, 0.0)
        """
        pass
    
