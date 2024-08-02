import sys
import traceback
import importlib

def test_import():
    import proteus
    successful_import=True
    for m in proteus.__all__:
        try:
            print(m)
            module = importlib.import_module("proteus.{0}".format(m))
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(repr(traceback.extract_tb(exc_traceback)))
            print("Failed to import proteus.{0}".format(m))
            successful_import=False
    assert(successful_import)

if __name__ == '__main__':
    test_import()
