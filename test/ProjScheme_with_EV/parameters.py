from proteus import Context

ct=Context.Options([
    ("USE_SUPG_NS",0,"Use SUPG stabilization"),
    ("ARTIFICIAL_VISCOSITY_NS",2,"0: no art visc., 1: shock capturing, 2: entropy-viscosity"),
    ("INT_BY_PARTS_PRESSURE",0,"int or not by parts pressure in momentum equation")
],mutable=True)