# MIT License
#
# Copyright (c) 2023, OASCI
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

import logging
import random
import time

import numpy as np

from . import _version

__version__ = _version.get_versions()["version"]


def set_log_level(level):
    """Dynamically control the log level of reptar.

    Parameters
    ----------
    level : :obj:`int`
        The desired logging level.
    """
    for logger_name in logging.root.manager.loggerDict:  # pylint: disable=no-member
        if "reptar" in logger_name:
            logger = logging.getLogger(logger_name)
            logger.setLevel(level)
            for handler in logger.handlers:
                handler.setLevel(level)


class ReptarLogger:
    """Logger class for all reptar modules"""

    def __init__(self, name, level=logging.INFO):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)

        # '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        self.formatter = logging.Formatter("%(message)s")
        self.handler = logging.StreamHandler()
        self.handler.setLevel(level)
        self.handler.setFormatter(self.formatter)
        self.logger.addHandler(self.handler)

        self.t_hashes = {}

    def log(self, level, msg, *args, **kwargs):
        self.logger.log(level, msg, *args, **kwargs)

    def debug(self, msg, *args, **kwargs):
        self.logger.log(logging.DEBUG, msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self.logger.log(logging.INFO, msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self.logger.log(logging.WARNING, msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self.logger.log(logging.ERROR, msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self.logger.log(logging.CRITICAL, msg, *args, **kwargs)

    def log_array(self, array, level=logging.INFO):
        if isinstance(array, (list, tuple)):
            array = np.array(array)
        arr_str = np.array2string(array, separator=", ")
        self.logger.log(level, arr_str)

    def t_start(self):
        r"""Record the time and return a random hash."""
        t_hash = random.getrandbits(128)
        self.t_hashes[t_hash] = time.time()
        return t_hash

    def t_stop(self, t_hash, message="Took {time} s", precision=5, level=20):
        r"""Determine timing from a hash.

        Parameters
        ----------
        t_hash : :obj:`str`
            Timing hash generated from ``TimeTracker.start()``.
        message : :obj:`str`, default: ``'Took {time} s'``
            Timing message to be written to log. ``'{time}'`` will be replaced
            with for the elapsed time.
        precision : :obj:`int`, default: ``5``
            Number of decimal points to print.
        level : :obj:`int`, default: ``20``
            Log level.
        """
        t_stop = time.time()
        t_elapsed = t_stop - self.t_hashes[t_hash]
        del self.t_hashes[t_hash]
        # pylint: disable-next=no-member
        self.log(level, message.replace("{time}", f"%.{precision}f" % t_elapsed))
        return t_elapsed
