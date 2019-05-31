import os

if 'GMXLIB' not in os.environ:
    here = os.path.abspath(os.path.dirname(__file__))
    os.environ['GMXLIB'] = os.path.join(here, 'top')
#print(os.environ['GMXLIB'])
    
    
from ._version import __version__, gmx_version
from . import main
from .holeOutputProcessor import HoleOutputProcessor
