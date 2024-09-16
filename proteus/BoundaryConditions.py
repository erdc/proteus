#!python
# distutils: language = c++
# cython: profile=True, binding=True, embedsignature=True
# cython: wraparound=False
# cython: boundscheck=False
# cython: initializedcheck=False
# cython: language_level=3
import cython
import numpy as np
"""
Module for creating boundary conditions. Imported in SpatialTools.py

.. inheritance-diagram:: proteus.BoundaryConditions
   :parts: 1
"""

__all__ = ['BC_Base',
           'BoundaryCondition']

class BC_Base:
    """
    Generic class regrouping boundary conditions
    """
    def __init__(self, shape=None, name=None, b_or=None, b_i=0, nd=None):
        self.Shape = shape
        self.name = name
        self.BC_type = 'None'
        if shape is not None:
          self.nd = self.Shape.Domain.nd
        elif nd is not None:
            self.nd = nd
        else:
            assert nd is not None, 'Shape or nd must be passed to BC'
        if b_or is not None:
            self._b_or = b_or[b_i]  # array of orientation of all boundaries of shape
        else:
            self._b_or = None

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


class BoundaryCondition:
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
        value : float
                Constant value

        """
        self.uOfXT = lambda x, t, n=np.zeros(3,): value


    def setLinearBC(self, a0, a):
        """
        function returning value=a0+ax*x+ay*y+az*z

        Parameters
        ----------
        a0 : float
             constant
        a: numpy.ndarray
            ax,ay,az

        """
        
        self.uOfXT = lambda x, t, n=np.zeros(3,): a0+sum(a[:]*x[:])

    def setLinearRamp(self,t1,value):
        """
        function setting a linear ramp from t=0 to t=t1

        Parameters
        -----------
        t1: float
            Ramp end time
        value: float
            Variable value
        """
        self.uOfXT = lambda x, t, n=np.zeros(3,): min( t/t1, 1)*value