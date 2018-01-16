import os, sys
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq

def ContextObject():
    from collections import namedtuple
    globalSettings = {"nnx":11, "T":10.0, "g":9.8}
    MyContext = namedtuple("MyContext",globalSettings.keys())
    return MyContext._make(globalSettings.values())

def check_eq(context):
    eq(context.nnx,11)
    eq(context.T,10.0)
    eq(context.g,9.8)

def test_set():
    from proteus import Context
    Context.set(ContextObject())
    check_eq(Context.context)

def test_setFromModule():
    import os
    from proteus import Context
    with open("context_module.py","w") as f:
        f.write("nnx=11; T=10.0; g=9.8\n")
    sys.path.append(os.getcwd())
    import context_module
    os.remove("context_module.py")
    Context.setFromModule(context_module)
    check_eq(Context.context)
    
def test_setMutableFromModule():
    import os
    from proteus import Context
    with open("context_module.py","w") as f:
        f.write("nnx=11; T=10.0; g=9.8\n")
    sys.path.append(os.getcwd())
    import context_module
    os.remove("context_module.py")
    Context.setFromModule(context_module, mutable=True)
    check_eq(Context.context)
    ct = Context.get()
    ct.T=11.0
    eq(ct.T,11.0)

def test_get():
    from proteus import Context
    Context.set(ContextObject())
    ct = Context.get()
    check_eq(ct)
    try:
        ct.T=11.0
    except Exception as e:
        assert(type(e) is AttributeError)

def test_Options():
    import os
    from proteus import Context
    Context.contextOptionsString="nnx=11"
    with open("context_module.py","w") as f:
        f.write('from proteus import Context; opts=Context.Options([("nnx",12,"number of nodes")]); nnx=opts.nnx; T=10.0; g=9.8\n')
    sys.path.append(os.getcwd())
    import context_module
    os.remove("context_module.py")
    Context.setFromModule(context_module)
    check_eq(Context.context)

if __name__ == '__main__':
    import nose
    nose.main()
