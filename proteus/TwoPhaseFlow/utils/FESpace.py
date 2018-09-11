"""
Create FE Spaces. 
"""
from proteus import FemTools as ft
    
class FESpace:
    def __init__(self,ns_model,nd):
        assert ns_model == 'rans2p' or ns_model == 'rans3p', 'ns_model must be rans2p or rans3p'
        assert nd in [2,3], 'number of dimensions must be 2 or 3'
        self.ns_model=ns_model
        self.nd=nd
        # For now we just support rans2p with: p1-p1 and rans3p with: p2-p1
        if ns_model == 'rans2p':
            self.velSpaceOrder=1
            self.pSpaceOrder=1
        else:
            self.velSpaceOrder=2
            self.pSpaceOrder=1
            
    def getFESpace(self):
        ##################
        # VELOCITY SPACE #
        ##################
        if self.velSpaceOrder == 1: # p1 space
            hFactor = 1.0
            velBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        else: # p2 space
            hFactor = 0.5
            velBasis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        ##################
        # PRESSURE SPACE #
        ##################
        if self.pSpaceOrder == 1: # p1 space
            pBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        else: # p2 space
            pBasis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        ###################
        # LEVEL SET SPACE #
        ###################
        lsBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis # p1 space
        ###################
        # QUADRATURE RULE #
        ###################
        if max(self.velSpaceOrder,self.pSpaceOrder)==1:
            elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 3)
            elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 3)
        else:
            elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 5)
            elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 5)
        #
        return {'lsBasis': lsBasis,
                'hFactor': hFactor,
                'velBasis': velBasis,
                'pBasis': pBasis,
                'elementQuadrature': elementQuadrature,
                'elementBoundaryQuadrature': elementBoundaryQuadrature}
