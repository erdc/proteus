from collections import namedtuple

Level = namedtuple('Level', ['variableNames', 'nodeArray', 'u', 't', 'mesh', 'coefficients'])


class MockModel:
    """ Simple mock model for testing purposes.

    Currently, this object only provides the functionality needed for testing Gauges.

    A model contains a mesh and a series of coefficients as well as their corresponding finite element spaces.
    The elements create a non-overlapping partition, while ghost dofs are distributed and updated where needed.

    See Transport.py for more details.
"""

    def __init__(self, variableNames, nodeArray, u, t, mesh, coefficients):
        l = Level(variableNames, nodeArray, u, t, mesh, coefficients)
        self.levelModelList = [l]