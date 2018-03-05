from proteus import Context

ct=Context.Options([
    ("NS_STABILIZATION_TYPE",0,"Stab on NS. 0: SUPG, 1, 2: EV via weak and strong residual"),
    ("nd",2,"Physical dimensioln")
],mutable=True)
