#!/usr/bin/env python

def printFinalSolutionToFile(resFileName,outputFileName,
                             key='u_dof',component=0,meshLevel=0,verbose=0):
    """
    write solution component at last (or only time step) on a given mesh level
    as simple text file

    """
    import os
    if not os.path.exists(resFileName):
        print("""resFileName= %s not found! """ % resFileName)
        return True
    import shelve
    results = shelve.open(resFileName)
    try:
        value = results['solutionData'][component][meshLevel][key]
        if hasattr(value,"tofile"):
            value.tofile(outputFileName,sep="\n",format="%12.5e")
        else:
            output = open(outputFileName,'w')
            for val in value:
                output.write("%12.5e \n" % val)
            output.close()
        return False
    except KeyError:
        print("""results['solutionData'][%s][%s][%s] not found """ % (component,
                                                                   meshLevel,
                                                                   key))
        if verbose > 0:
            print("""results.keys() = %s """ % list(results.keys()))
            print("""results['solutionData'].keys() = %s """ % list(results['solutionData'].keys()))
        return True
    return True

if __name__ == '__main__':
    import optparse
    usage = "usage: %prog [options] resultsFileName outputFileName"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--component",
                      help="solution component",
                      action="store",
                      type="int",
                      default=0)
    parser.add_option("-k", "--keyname",
                      help="key name for solution",
                      action="store",
                      type="string",
                      default="u_dof")
    parser.add_option("--level",
                      help="mesh level number for solution",
                      action="store",
                      type="int",
                      default=0)
    parser.add_option("-v", "--verbose",
                      help="Print execution information to standard out",
                      action="store",
                      type="int",
                      default=0)


    (opts,args) = parser.parse_args()
    if len(args) < 2:
        raise RuntimeError("command line requires resultsFileName outputFileName")
    resultsFileName = args[0]
    outputFileName  = args[1]

    if opts.verbose > 0:
        print("""calling with resultsFileName=%s outputFileName=%s
        key=%s component=%s level=%s """ % (resultsFileName,
                                            outputFileName,
                                            opts.keyname,
                                            opts.component,
                                            opts.level))
    failed = printFinalSolutionToFile(resultsFileName,outputFileName,
                                      key=opts.keyname,
                                      component=opts.component,
                                      meshLevel=opts.level,
                                      verbose=opts.verbose)
