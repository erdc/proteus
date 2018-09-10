"""
FE Spaces for 2PhaseFlow
"""
from proteus import FemTools as ft

class FESpace:
    def __init__(self,nd=2,spaceOrder=1):
        assert nd in [2,3], 'number of dimensions must be 2 or 3'
        assert spaceOrder in [1,2], 'spaceOrder must be 1 or 2'
        self.nd=nd
        self.spaceOrder=spaceOrder
    def getFESpace(self):
        # p1 space
        hFactor = 1.0
        basis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 3)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 3)
        if self.spaceOrder == 2: # p2 space
            hFactor = 0.5
            basis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
            elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 4)
            elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 4)        
        #
        return {'hFactor': hFactor,
                'basis': basis,
                'elementQuadrature': elementQuadrature,
                'elementBoundaryQuadrature': elementBoundaryQuadrature}
