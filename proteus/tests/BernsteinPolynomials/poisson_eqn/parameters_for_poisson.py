from proteus import Context

##############
# PARAMETERS #
##############
ct=Context.Options([
    ("nd",2,"Number of dimensions"),
    ("pDegree",2,"Order of the polynomial approximation"),
    ("refinement",0,"Mesh refinement"),
    ("useHex",True,"Use quads?"),
    ("useBernstein",True,"Use Bernstein polynomials"),
    ("unstructured",False,"Use unstructured triangular mesh"),
    ("genMesh", False, "Generate a new mesh?")
],mutable=True)
