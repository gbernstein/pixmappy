# Test that PixelMapCollection works properly if galsim is not available.

from __future__ import print_function

# Mock up the import statement to fail for 'import galsim'
# cf. http://stackoverflow.com/questions/2481511/mocking-importerror-in-python
# Although I had to tweak their solution slightly.
try:
    import __builtin__ as builtins
except ImportError:
    import builtins
realimport = builtins.__import__
def myimport(name, globals={}, locals={}, fromlist=[], level=0):
    if name == 'galsim':
        raise ImportError("galsim is not available on purpose")
    return realimport(name, globals, locals, fromlist, level)
builtins.__import__ = myimport

def test_tpv():
    from test_reference import test_tpv
    test_tpv()

def test_complex():
    from test_reference import test_complex
    test_complex()

if __name__ == '__main__':
    test_tpv()
    test_complex()
