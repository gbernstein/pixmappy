from importlib.metadata import version 
__version__ = version('pixmappy')
__version_info__ = tuple(map(int, __version__.split('.')))
