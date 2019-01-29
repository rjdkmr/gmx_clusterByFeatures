import os

if 'GMXLIB' not in os.environ:
    os.environ['GMXLIB'] = os.path.abspath('top')
    
from ._version import __version__, gmx_version
from . import main
