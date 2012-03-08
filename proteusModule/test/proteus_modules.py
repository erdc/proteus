import sys

def test_import():
    import proteus
    successful_import=True
    for m in proteus.__all__:
        try:
            print m
            module = __import__('proteus.'+m,fromlist=['proteus'])
            del module
        except:
            print "Failed to import ",m
            successful_import=False
    del proteus
    assert(successful_import)

if __name__ == '__main__':
    test_import()
