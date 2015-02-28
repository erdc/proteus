"""A global user context object"""
from collections import  namedtuple
from .Profiling import logEvent

"""The global context"""
context = None

contextOptionsString = None

def get():
    global context
    return context

def set(contextIn):
    """Set the global context as an object"""
    global context
    context = contextIn

def setFromModule(moduleIn):
    """Construct the context object from a module"""
    global context
    fields = {}
    for key,value in moduleIn.__dict__.iteritems():
        if  key[:2] != "__":
            fields[key]=value
    Context = namedtuple(moduleIn.__name__,fields.keys(),verbose=False)
    context = Context._make(fields.values())

def declareAndGetInputOptions(optionsList=None):
    """Return a named tuple of the context options"""
    global contextOptionsString,globalContextOptions
    contextOptionsDict = {}
    help="Context input options:\n"
    for optTuple  in optionsList:
        print optTuple
        help += """{0}[{1}] {2}\n\n""".format(optTuple[0],optTuple[1],optTuple[2])
        contextOptionsDict[optTuple[0]] = optTuple[1]
    logEvent(help)
    if contextOptionsString=="?":
        print(help)
        contextOptionsString=None
    if contextOptionsString != None:
        option_overides=contextOptionsString.split(" ")
        for option in option_overides:
            lvalue,rvalue = option.split("=")
            if contextOptionsDict.has_key(lvalue):
                logEvent("Processing context input options from commandline")
                try:
                    contextOptionsDict[lvalue] = eval(rvalue)
                    logEvent(lvalue+" = "+rvalue)
                except:
                    print "Failed setting context options from command line string."
                    raise 
            else:
                logEvent("IGNORING CONTEXT OPTION; DECLARE "+lvalue+" IF YOU WANT TO SET IT")
    #now set named tuple merging optionsList and opts_cli_dihelpct
    ContextOptions = namedtuple("ContextOptions",contextOptionsDict.keys(),verbose=False)
    return ContextOptions._make(contextOptionsDict.values())
