from __future__ import print_function
import sys
import traceback

def test_import():

    import proteus
    successful_import=True
    for m in proteus.__all__:
        try:
            module = __import__('proteus.'+m,fromlist=['proteus'])
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(repr(traceback.extract_tb(exc_traceback)))
            print("Failed to import proteus.",m)
            successful_import=False
    assert(successful_import)

if __name__ == '__main__':
    test_import()
