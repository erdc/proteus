def test_import():
    import proteus
    successful_import=True
    for m in proteus.__all__:
        print m
        try:
            module = __import__('proteus.'+m,fromlist=['proteus'])
        except:
            print "Failed to import ",m
            successful_import=False
    assert(successful_import)

if __name__ == '__main__':
    test_import()
