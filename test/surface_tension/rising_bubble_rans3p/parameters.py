from proteus import Context

ct=Context.Options([
    ("USE_SUPG_NS",0,"Use SUPG stabilization"),
    ("ARTIFICIAL_VISCOSITY_NS",2,"0: no art visc., 1: shock capturing, 2: entropy-viscosity"),
    ("nd",2,"Physical dimensioln")
],mutable=True)