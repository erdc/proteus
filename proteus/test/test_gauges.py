import numpy as np
from .mock_model import MockModel

variableNames = ['u0']

nodeArray = np.array([[0, 0, 0],
                      [1, 0, 0],
                      [0, 1, 0],
                      [1, 1, 0]])


# The majority of PointGauges.buildQuantityRow should be moved to either
# femTools or meshTools.  This will simplify test writing
# as that function can then be safely mocked.


# u = TODO
t = 0
# mesh = TODO
# coefficients = TODO

# model = MockModel(variableNames, nodeArray, u, t, mesh, coefficients)