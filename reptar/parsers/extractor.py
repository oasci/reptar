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

from abc import ABC, abstractmethod

class extractor(ABC):
    """Base class for extracting information from files.

    Attributes
    ----------
    triggers : :obj:`tuple`
        A collection of triggers that activate the corresponding extractor.
        The trigger is a lambda function that returns True or False depending
        on the criteria and the name of the extractor method.
    parsed_info : :obj:`dict`
        Information parsed from files. Contains the following keys.

        ``system_info``
            Information specifying the system prior to any computation. Such
            as the initial cartesian coordinates, total system charge and
            multiplicity, etc.
        
        ``runtime_info``
            Contains information about setting up the job/calculation or running
            the job. Defining convergence criteria, parameters, etc.
        
        ``outputs``
            Results, requested or not, of the job. For example, SCF
            cycle values, optimized coordinates, trajectory, number of
            electrons, generated structures, etc.
    """
    def __init__(self):
        self.parsed_info = {
            'system_info': {},
            'runtime_info': {},
            'outputs': {}
        }

    @property
    @abstractmethod
    def triggers(self):
        pass
    