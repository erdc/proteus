class InitialCondition:
    """
    Initial condition class

    Attributes
    ----------
    uOfXT: func or None
        initial condition function of x (array_like) and t (float) or None for
        no initial condition
    """
    __frozen = False
    
    def __init__(self):
        self.uOfXT = None
