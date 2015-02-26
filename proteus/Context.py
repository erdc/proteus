"""A global user context object"""
from collections import  namedtuple

"""The global context"""
context = None

def set(contextIn):
    """Set the global context as an object"""
    global context
    context = contextIn

def get():
    global context
    return context

def setFromModule(moduleIn):
    """Construct the context object from a module"""
    global context
    fields = {}
    for key,value in moduleIn.__dict__.iteritems():
        if  key[:2] != "__":
            fields[key]=value
    Context = namedtuple(moduleIn.__name__,fields.keys(),verbose=False)
    context = Context._make(fields.values())
