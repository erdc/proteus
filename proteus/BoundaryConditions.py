"""
Module for creating boundary conditions. Imported in SpatialTools.py

.. inheritance-diagram:: proteus.BoundaryConditions
   :parts: 1
"""


class BC_Base(object):
    """
    class defining boundary conditions for a facet or a segment
    """
    def __init__(self, shape=None, name=None, b_or=None, b_i=None):
        self.Shape = shape
        self.name = name
        self.BC_type = 'None'
        self._b_or = b_or  # array of orientation of all boundaries of shape
        self._b_i = b_i  # indice for this boundary in list of boundaries

    @classmethod
    def newGlobalBC(cls, name, default_value=None):
        """
        Makes a new boundary condition with a default value. This creates a new
        class instance and adds it to all BoundaryConditions class instances.
        :param name: name of the new class attribute
        :param default: default value of BC for attr (usually a funtion of
                        (x, t) or None)
        """
        setattr(cls, name, default_value)


def constantBC(value):
    """
    function returning constant BC
    """
    return lambda x, t: value

def linearBC(a0, a1, i):
    """
    function returning linear BC
    :param a0:
    :param a1:
    :param i:
    """
    return lambda x, t: a0 + a1*x[i]
