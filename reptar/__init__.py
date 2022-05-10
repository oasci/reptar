__author__ = """Alex M. Maldonado"""

from . import _version
__version__ = _version.get_versions()['version']

from .creator import creator
from .reptar_file import File

__all__ = ['creator', 'File']
