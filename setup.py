from __future__ import print_function
import sys,os,glob,re
import select

try:
    from setuptools import setup
    import setuptools
    print("Using setuptools version",setuptools.__version__)
except ImportError:
    from distutils.core import setup
    import distutils
    print("Using distutils version",distutils.__version__)

print('Python version = ',sys.version)
py_version = "%d.%d"%sys.version_info[0:2]  # we check things based on the major.minor version.

dependencies = ['numpy', 'future', 'astropy', 'scipy', 'pyyaml', 'LSSTDESC.Coord']

with open('README.md') as file:
    long_description = file.read()

# Read in the version from pixmappy/_version.py
# cf. http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
version_file=os.path.join('pixmappy','_version.py')
verstrline = open(version_file, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    pixmappy_version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (version_file,))
print('PixMapPy version is %s'%(pixmappy_version))

data = glob.glob(os.path.join('data','*'))

dist = setup(
        name="PixMapPy",
        version=pixmappy_version,
        author="Gary Bernstein",
        author_email="garyb@PHYSICS.UPENN.EDU",
        description="Python module for arbitrary mappings from pixels to sky coordinate",
        long_description=long_description,
        license = "BSD License",
        url="https://github.com/gbernstein/pixmappy",
        download_url="https://github.com/gbernstein/pixmappy/releases/tag/v%s.zip"%pixmappy_version,
        packages=['pixmappy'],
        package_data={'pixmappy' : data },
        install_requires=dependencies,
    )

