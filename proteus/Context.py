"""Support for a global user context

Sub-models in multi-physics simulations sometimes need to be consistent
with other sub-models through certain "global" information. This model
provides methods for setting and getting such a global user context.
It also alows declaring input options to the user context that can be
set from the command line.

Example (set the global context from an object)::

  #create a simple context class
  globalSettings = {"nnx":11, "T":10.0, "g"=9.8}
  MyContext = namedtuple("MyContext",globalSettings.keys())
  #set the global  context  object
  proteus.Context.set(MyContext._make(globalSettings.values())

Example (set the global context from a module)::

  import globalSettingsModule
  proteus.Context.setFromModule(globalSettingsModule)

Example (use the global context)::

  ct = proteus.Context.get()
  nnx = ct.nnx

"""
from collections import  namedtuple
from .Profiling import logEvent

"""The global context"""
context = None

"""A string containing Python setting input options"""
contextOptionsString = None

def get():
    """Get the global context object"""
    global context
    return context

def set(contextIn):
    """Set the global context as an object"""
    global context
    context = contextIn

def setFromModule(moduleIn):
    """Construct the global context object from a module"""
    global context
    fields = {}
    for key, value in moduleIn.__dict__.iteritems():
        if  key[:2] != "__":
            fields[key]=value
    Context = namedtuple(moduleIn.__name__.split('.')[-1], fields.keys(), verbose=False)
    context = Context._make(fields.values())

def Options(optionsList=None):
    """Construct an o
    from proteus.LinearAlgebraToptions object (named tuple)
    
    :param optionsList: An iteratble of options tuples. Each option is
                        a 3-tuple, with the first entry the option
                        name string, the second entry the default
                        value, and the third entry a help string

    Example::

      opts=Context.Options([("nnx", 11, "number of mesh nodes in x-direction"),
                            ("nny", 21, "number of mesh nodes in y-direction"])
      nnx = opts.nnx
      nny = opts.nny

    """
    import ast, sys
    global contextOptionsString
    contextOptionsDict = {}
    help="Context input options:\n"
    for optTuple  in optionsList:
        help += """{0}[{1}] {2}\n\n""".format(optTuple[0], optTuple[1], optTuple[2])
        contextOptionsDict[optTuple[0]] = optTuple[1]
    logEvent(help)
    if contextOptionsString=="?":
        print(help)
        contextOptionsString=None
    if contextOptionsString is not None:
        option_overides=contextOptionsString.split(" ")
        for option in option_overides:
            lvalue, rvalue = option.split("=")
            if contextOptionsDict.has_key(lvalue):
                logEvent("Processing context input options from commandline")
                try:
                    contextOptionsDict[lvalue] = ast.literal_eval(rvalue)
                    logEvent(lvalue+" = "+rvalue)
                except:
                    sys.stderr.write("Failed setting context options from command line string.")
                    raise 
            else:
                logEvent("IGNORING CONTEXT OPTION; DECLARE "+lvalue+" IF YOU WANT TO SET IT")
    #now set named tuple merging optionsList and opts_cli_dihelpct
    ContextOptions = namedtuple("ContextOptions", contextOptionsDict.keys(), verbose=False)
    return ContextOptions._make(contextOptionsDict.values())
