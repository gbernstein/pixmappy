
from ._version import __version__, __version_info__

import os

root_dir = os.path.dirname(__file__)
data_dir = os.path.join(root_dir, 'data')

from .PixelMapCollection import *
