import cython
"""
Module for creating boundary conditions. Imported in SpatialTools.py

.. inheritance-diagram:: proteus.BoundaryConditions
   :parts: 1
"""

class BC_Base():
    """
    Generic class regrouping boundary conditions
    """
    def __init__(self, shape=None, name=None, b_or=None, b_i=0):
        self.Shape = shape
        self.name = name
        self.BC_type = 'None'
        if b_or is not None:
            print(b_or)
            self._b_or = b_or[b_i]  # array of orientation of all boundaries of shape

    # @staticmethod
    # def newGlobalBC(name, default_value):
    #     """
    #     Makes a new BoundaryCondition class instance attached to BC_Base.
    #     This creates a new class attributed that will be added to all BC_Base
    #     class instances.

    #     Parameters
    #     ----------
    #     name: string
    #     """
    #     setattr(BC_Base, name, default_value)

    def getContext(self, context=None):
        """
        Gets context from proteus.Context or

        Parameters
        ----------
        context: class, optional
             if set to None, the context will be created from proteus.Context
        """
        if context:
            from proteus import Context
            self.ct = Context.get()
        else:
            self.ct = context


class BoundaryCondition():
    """
    Boundary condition class

    Attributes
    ----------
    uOfXT: func or None
        boundary condition function of x (array_like) and t (float) or None for
        no boundary condition
    """
    def __init__(self):
        self.uOfXT = None

    def init_cython(self):
        return self.uOfXT

    def resetBC(self):
        self.uOfXT = None

    def setConstantBC(self, value):
        """
        function returning constant BC

        Parameters
        ----------
        value :

        """
        self.uOfXT = lambda x, t: value


    def setLinearBC(self, a0, a1, i):
        self.uOfXT = lambda x, t: a0+a1*x[i]
