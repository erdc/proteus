#!/usr/bin/env/python

from EGeometry import *
from MeshTools import *
from FemTools import *
from LinearAlgebraTools import *
from LinearSolvers import *
from ScalarTransport import *
from NormTools import *
import PoissonTestProblems
import DiagUtils
## \defgroup P1ncMixed P1ncMixed
#
# @{

"""
Implement P^1 nonconforming FEM approximations using Crouzeix-Raviart elements

Include postprocessing techniques to get equivalence with RT_0 MHFEM for piecewise
constant data

"""



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#higher level post-processing routines like getting RT0 fluxes and 
#potential
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def postProcessRT0flux(uh,A,f,nd,fluxdofs,verbose=0):
    """
    follow example in Chen '94 paper to
    take scalar finite element function, assuming it's from NC_Affine ...
    space and compute RT0 flux,  q_h \approx -A\del u
    A is the tensor conductivity and f is the scalar source term from
    the mass conservation equation. These have to be evaluated as the
    'averages' aka L_2 projections on element-wise constants

    For RT0 on a triangle have
    
      \vec q_h = \vec a_T + b_t\vec x

    which turns out to be on triangle T

      \vec q_h = -\mat{A}_h\grad u_h + \bar{f}_T/2(\vec x - \vec x_{T})

    so
      b_T = \bar{f}_T/d
    and
        \vec a_T = -\mat{A}_{h,T}\grad u_h - \bar{f}_{T}/d \vec x_{T}
    where
        \bar{f}_{T} is the average source over element T
        \mat{A}_{h} is the constant, projected conductivity over T
        \vec x_{T}  is the element barycenter
         d is the spatial dimension
    assumes that solution holds all degrees of freedom and boundary
    values.

    returns array fluxdofs that's nelements by nd+1, that stores
       fluxdofs[eN,:] = [\vec a_T,b_T]
    """
    ebary = Numeric.array([1./3.,1./3.,0.])
    if nd == 3:
        ebary = Numeric.array([1./4.,1./4.,1./4.])
    #end if
    ne    = uh.femSpace.elementMaps.mesh.nElements_global
    ndofs = uh.femSpace.referenceFiniteElement.localFunctionSpace.dim
    n_xi  = 1
    xiArray = Numeric.zeros((n_xi,nd),Numeric.Float)
    xiArray[0,:] = ebary[0:nd]
    if verbose > 4:
        print 'ebary[0:nd]= ',ebary[0:nd],'\n xiArray= ',xiArray
    #go back and recalculate element mapping info for now, then
    #compare to what's already been calculated in ScalarTransport system
    jacArray    = Numeric.zeros((ne,n_xi,nd,nd),Numeric.Float)
    jacInvArray = Numeric.zeros((ne,n_xi,nd,nd),Numeric.Float)
    jacDetArray = Numeric.zeros((ne,n_xi),Numeric.Float)
    xbary       = Numeric.zeros((ne,n_xi,nd),Numeric.Float)
    
    gradArray   = Numeric.zeros((ne,n_xi,ndofs,nd),Numeric.Float)

    #map to element barycenters in physical space
    uh.femSpace.elementMaps.getValues(xiArray,xbary)
    #get jacobian information for elements
    uh.femSpace.elementMaps.getJacobianValues(xiArray,jacArray,jacInvArray,jacDetArray)
    #shape function gradients in physical space
    uh.femSpace.getBasisGradientValues(xiArray,jacInvArray,gradArray)

    #for local evaluations
    Ah    = Numeric.zeros((nd,nd),Numeric.Float)
    aT    = Numeric.zeros(nd,Numeric.Float)
    fh = 0.0
    bT = 0.0
    #do loop over elements
    for eN in range(ne):
        #get A_h and f_h
        Ah = A(xbary[eN,0,:])
        fh = f(xbary[eN,0,:])
        #get gradient
        gradu = Numeric.zeros(nd,Numeric.Float)
        for j in uh.femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
            J = uh.femSpace.dofMap.l2g[eN,j]
            gradu[:] += uh.dof[J]*gradArray[eN,0,j,:]
        #end grad u calc
        bT = fh/float(nd)
        aT = -Numeric.dot(Ah,gradu) -fh/float(nd)*xbary[eN,0,:]
        fluxdofs[eN,0:nd]= aT
        fluxdofs[eN,nd]  = bT
    #end element loop
#end postprocess flux
def postProcessRT0potential(uh,A,f,nd,fluxdofs,uRT0dof,verbose=0):
    """
     follow example in Chen '94 paper to get RTO potential value.
     Assumes RT0 flux has already been calculated!
     
     here, we basically have to plug in the P^1_{nc} solution,
     and the postprocessed flux into the local Darcy's law expression
     from the mixed hybrid formulation using the test function 
       \vec v_h= \vec x. 
     Then we solve for the pressure unknown.
     
     (\mat{A}_h^{-1}\vec q_h,\vec v_h)_T  - (p_h,\div \vec v_h)_T 
        + (\lambda_h,\vec v_h\cdot n_{T})_{\partial T}
     
     Recall, that \lambda^{j}_h = p^j_h, the P^1_{nc} solution for edge j
    
     To evaluate the weak integrals, I'll use the quadrature formula
      \int_{T} f \dx \approx |T|/3 \sum_{j}f(\vec \bar{x}^j)
    
     for the mass matrix term, and basically exact integration (midpoint rule)
     for the boundary integrals of the Lagrange multipler term since it's
     linear on each edge
    
    """
    #for numerical quadrature, could just use built in quadrature module
    n_xi  = nd+2
    xiArray  = Numeric.zeros((n_xi,nd),Numeric.Float)
    n_bxi = 1
    bxiArray = Numeric.zeros((n_bxi,nd),Numeric.Float)
    qweight = Numeric.ones((n_xi,),Numeric.Float)*1.0/6.0
    nquad   = nd+1
    bqweight= Numeric.ones((n_bxi,),Numeric.Float)
    nbquad  = 1
    refVol  = 0.5
    if nd == 1:
        qa = 0.5*sqrt(3.)/5.
        xiArray[0,:] = Numeric.array([0.5])
        xiArray[1,:] = Numeric.array([0.5-qa])
        xiArray[2,:] = Numeric.array([0.5+qa])

        bxiArray[0,:]= Numeric.array([0.0])
        refVol  = 1.0

    if nd == 2:
        xiArray[0,:] = Numeric.array([1./3.,1./3.])
        
        xiArray[1,:] = Numeric.array([0.5,0.5])
        xiArray[2,:] = Numeric.array([0.,0.5])
        xiArray[3,:] = Numeric.array([0.5,0.])
        
        bxiArray[0,:]= Numeric.array([0.5,0.0])
    elif nd == 3:
        xiArray[0,:] = Numeric.array([1./4.,1./4.,1./4.])

        qa = (5.0-sqrt(5.0))/20.0
        qb = 1.0-3.0*qa
        qweight = Numeric.ones((n_xi,),Numeric.Float)*0.25*1.0/6.0
        refVol  = 1.0/6.0
        xiArray[1,:] = Numeric.array([qa,qa,qa])
        xiArray[2,:] = Numeric.array([qb,qa,qa])
        xiArray[3,:] = Numeric.array([qa,qb,qa])
        xiArray[4,:] = Numeric.array([qa,qa,qb])
        bxiArray[0,:]= Numeric.array([1.0/3.0,1.0/3.0,0.0])
        bqweight     = Numeric.ones((n_bxi,),Numeric.Float)*0.5
    #end 3d
    ne    = uh.femSpace.elementMaps.mesh.nElements_global
    ndofs = uh.femSpace.referenceFiniteElement.localFunctionSpace.dim
    nebLoc= uh.femSpace.elementMaps.mesh.nElementBoundaries_element
    if verbose > 4:
        print 'ebary[0:nd]= ',ebary[0:nd],'\n xiArray= ',xiArray
    #go back and recalculate element mapping info for now, then
    #compare to what's already been calculated in ScalarTransport system
    jacArray    = Numeric.zeros((ne,n_xi,nd,nd),Numeric.Float)
    jacInvArray = Numeric.zeros((ne,n_xi,nd,nd),Numeric.Float)
    jacDetArray = Numeric.zeros((ne,n_xi),Numeric.Float)
    xArray      = Numeric.zeros((ne,n_xi,nd),Numeric.Float)
    #get information for boundaries too
    #indexed like [nelems, nelemBndsLoc,nQpts,nSpaceDim]
    bxArray            = Numeric.zeros((ne,nebLoc,n_bxi,nd),Numeric.Float)
    #indexed like [nelems, nelemBndsLoc,nQpts,nSpaceDim,nSpaceDim]
    bjacInvArray    = Numeric.zeros((ne,nebLoc,n_bxi,nd,nd),Numeric.Float)
    #indexed like [nelems, nelemBndsLoc,nQpts,nSpaceDim-1,nSpaceDim-1]
    metricTensArray    = Numeric.zeros((ne,nebLoc,n_bxi,max(nd-1,1),max(nd-1,1)),
                                       Numeric.Float)
    #indexed like [nelems, nelemBndsLoc,nQpts]
    metricTensDetArray = Numeric.zeros((ne,nebLoc,n_bxi),Numeric.Float)
    #indexed like [nelems, nelemBndsLoc,nQpts,nSpaceDim]
    unitNormalArray    = Numeric.zeros((ne,nebLoc,n_bxi,nd),Numeric.Float)

    #indexed like [nelems, nelemBndsLoc,nQpts,nElementDim]
    bshapeArray        = Numeric.zeros((ne,nebLoc,n_bxi,ndofs),Numeric.Float)


    #map to element barycenters and quadrature points in physical space
    uh.femSpace.elementMaps.getValues(xiArray,xArray)
    #get jacobian information for elements
    uh.femSpace.elementMaps.getJacobianValues(xiArray,jacArray,
                                              jacInvArray,jacDetArray)
    #map to element boundary quadrature points in physical space
    uh.femSpace.elementMaps.getValuesTrace(bxiArray,bxArray)
    #get jacobian information for elements boundary mappints
    uh.femSpace.elementMaps.getJacobianValuesTrace(bxiArray,
                                                   bjacInvArray,
                                                   metricTensArray,
                                                   metricTensDetArray,
                                                   unitNormalArray)
    #get value of shape functions at boundary quadrature points
    uh.femSpace.getBasisValuesTrace(bxiArray,bshapeArray)
    
    #for local evaluations
    Ah    = Numeric.zeros((nd,nd),Numeric.Float)
    AhInv = Numeric.zeros((nd,nd),Numeric.Float)
    qi    = Numeric.zeros(nd,Numeric.Float)
    AqInv  = Numeric.zeros(nd,Numeric.Float)
    

    #do loop over elements
    for eN in range(ne):
        #get A_h and f_h
        Ah = A(xArray[eN,0,:])
        AhInv = inv(Ah)
        #compute (\mat{A}_h^{-1}\vec q_h,\vec v_h)_T
        divsum = 0.0
        #recall that xArray starts with barycenter value
        for k in range(nquad):
            detT     = abs(jacDetArray[eN,k])
            qi[:]    = fluxdofs[eN,0:nd] + fluxdofs[eN,nd]*xArray[eN,1+k,0:nd]
            AqInv[:] = Numeric.dot(AhInv,qi)
            divsum  += Numeric.dot(AqInv[:],xArray[eN,1+k,0:nd])*qweight[k]*detT
            if verbose > 3:
                print 'ppRT0p eN= ',eN,' k= ',k,'\n qi= ',qi,' AqInv= ',AqInv
                print 'detT= ',detT,' divsum= ',divsum
            #end verbose
        #end k loop for div integral

        #could manually compute integral using face-midpoint quadrature which should
        #be exact here and
        #take advantage of correspondence between local ordering of shape functions
        #and faces

        #try testing out the pyadh facilities for doing
        #boundary integrals instead
        bndrySum = 0.0
        for i in range(nebLoc):
            #do not have boundary barycenter included in xi
            for k in range(n_bxi):
                uk = 0.0
                for j in range(ndofs):
                    J = uh.femSpace.dofMap.l2g[eN,j]
                    uk += uh.dof[J]*bshapeArray[eN,i,k,j]
                #end j
                xnik= Numeric.dot(bxArray[eN,i,k,:],unitNormalArray[eN,i,k,:])
                bndrySum += uk*xnik*bqweight[k]*abs(metricTensDetArray[eN,i,k])
                if verbose > 3:
                    print 'ppRT0p uh[eN,i]= ',uh.dof[uh.femSpace.dofMap.l2g[eN,i]]
                    print 'uk= ',uk
                    print 'normal[eN,i,k]= ',unitNormalArray[eN,i,k,:]
                    print 'x.normal[eN,i,k]= ',xnik
                #end verbose
            #end k
        #end i
        #now update potential
        vol = abs(jacDetArray[eN,0])*refVol
        denom = float(nd)*vol
        uRT0dof[eN] = (bndrySum + divsum)/denom
    #end element loop
#end postprocess potential

def postProcessP1flux(uh,A,f,nd,fluxdofs,verbose=0):
    """
    Just compute standard? elemental flux
      \vec q = -\mat{A}\grad u
    for P^1 approximation, so it's constant on an element

    To match what I'm doing with the RT0 flux, I'll write the
    values as
       \vec q_h = \vec a_T  + 0*\vec x
    and store
       the vector a_T and scalar b_T := 0

    assumes that solution holds all degrees of freedom and boundary
    values.

    returns array fluxdofs that's nelements by nd+1, that stores
       fluxdofs[eN,:] = [\vec a_T,b_T]
    """
    ebary = Numeric.array([1./3.,1./3.,0.])
    if nd == 3:
        ebary = Numeric.array([1./4.,1./4.,1./4.])
    #end if
    ne    = uh.femSpace.elementMaps.mesh.nElements_global
    ndofs = uh.femSpace.referenceFiniteElement.localFunctionSpace.dim
    n_xi  = 1
    xiArray = Numeric.zeros((n_xi,nd),Numeric.Float)
    xiArray[0,:] = ebary[0:nd]
    #go back and recalculate element mapping info for now, then
    #compare to what's already been calculated in ScalarTransport system
    jacArray    = Numeric.zeros((ne,n_xi,nd,nd),Numeric.Float)
    jacInvArray = Numeric.zeros((ne,n_xi,nd,nd),Numeric.Float)
    jacDetArray = Numeric.zeros((ne,n_xi),Numeric.Float)
    xbary       = Numeric.zeros((ne,n_xi,nd),Numeric.Float)
    
    gradArray   = Numeric.zeros((ne,n_xi,ndofs,nd),Numeric.Float)

    #map to element barycenters in physical space
    uh.femSpace.elementMaps.getValues(xiArray,xbary)
    #get jacobian information for elements
    uh.femSpace.elementMaps.getJacobianValues(xiArray,jacArray,jacInvArray,jacDetArray)
    #shape function gradients in physical space
    uh.femSpace.getBasisGradientValues(xiArray,jacInvArray,gradArray)

    #for local evaluations
    Ah    = Numeric.zeros((nd,nd),Numeric.Float)
    aT    = Numeric.zeros(nd,Numeric.Float)
    bT = 0.0
    #do loop over elements
    for eN in range(ne):
        #get A_h and f_h
        Ah = A(xbary[eN,0,:])
        #get gradient
        gradu = Numeric.zeros(nd,Numeric.Float)
        for j in uh.femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
            J = uh.femSpace.dofMap.l2g[eN,j]
            gradu[:] += uh.dof[J]*gradArray[eN,0,j,:]
        #end grad u calc
        aT = -Numeric.dot(Ah,gradu)
        fluxdofs[eN,0:nd]= aT
        fluxdofs[eN,nd]  = bT
    #end element loop
    
#end post process P1 flux
def postProcessP1avg(uh,A,f,nd,fluxdofs,uAvgDofs,verbose=0):
    """
    Just compute standard? elemental average for P^1 approximation,
      \bar{u}_h = sum_j=1^d u_h^j\psi_j(\bar{x}) = (sum_j=1^d u_h^j)/d
     so it's constant on an element

    To match what I'm doing with the RT0 flux, I'll write the
    values as
       \vec q_h = \vec a_T  + 0*\vec x
    and store
       the vector a_T and scalar b_T := 0

    assumes that solution holds all degrees of freedom and boundary
    values.

    returns array of averages nelements by 1
    """
    ne    = uh.femSpace.elementMaps.mesh.nElements_global
    denom = 1.0/float(nd+1)
    #do loop over elements
    for eN in range(ne):
        usum = 0.0
        for j in uh.femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
            J = uh.femSpace.dofMap.l2g[eN,j]
            usum += uh.dof[J]
        #end grad u calc
        uAvgDofs[eN]  = usum*denom
    #end element loop
    
#end post process P1 flux

def testLaplacianP1nc(mesh,Lx,Ly,Lz,nx,ny,nz,nd=2,useCG=False,verbose=0):
    """
    try to write some code to define a simple linear steady-state problem and solve it
    using the ScalarTransport interface and CrR elements

    solve
       -\deld u = f   in \Omega = [0,1]^d
              u = u^b   on \partial \Omega
              u = \sum^d_i=1 x_i^2

    """
    # #  setup problem coefficients 
    if verbose > 0:
        print 'definining problem coefficients'
    #end verbose
    if nd == 1:
        fProb   = PoissonTestProblems.fP5    
        AProb   = PoissonTestProblems.AP5    
        A0Prob  = DiagUtils.Ident1
        uexProb = PoissonTestProblems.uexP5  
        qexProb = PoissonTestProblems.qexP5  
        dirProb = PoissonTestProblems.dirBCP5
    elif nd == 3:
        fProb   = PoissonTestProblems.fP4     # fP3    
        AProb   = PoissonTestProblems.AP4     # AP3    
        A0Prob  = DiagUtils.Ident3
        uexProb = PoissonTestProblems.uexP4   # uexP3  
        qexProb = PoissonTestProblems.qexP4   # qexP3  
        dirProb = PoissonTestProblems.dirBCP4 # dirBCP3
    else:
        fProb   = PoissonTestProblems.fP1
        AProb   = PoissonTestProblems.AP1
        A0Prob  = DiagUtils.Ident2
        uexProb = PoissonTestProblems.uexP1
        qexProb = PoissonTestProblems.qexP1
        dirProb = PoissonTestProblems.dirBCP1
    #end else on space dim
    probCoefficients = PoissonTestProblems.PoissonCoefficients(A=A0Prob,S=fProb,nd=nd)
    exactSoln = PoissonTestProblems.poissonExactSoln(uexProb)
    qexactSoln= PoissonTestProblems.poissonExactSoln(qexProb)
    # # setup boundary conditions
    # # setup numerical quadrature (midpoint rules over elements and element boundaries
    if verbose > 0:
        print 'definining quadrature'
    #end verbose
    quadrature = {}
    if not useCG:
        quadratureOrder = 1
    else:
        quadratureOrder = 2
    #end quadrature order
    elemQuadrature = SimplexGaussQuadrature(nd)
    elemQuadrature.setOrder(quadratureOrder)
    for integral in OneLevelScalarTransport.integralKeys:
        quadrature[integral] = elemQuadrature
    #end for
    elementBoundaryQuadrature = {}
    if not useCG:
        bquadratureOrder = 1
    else:
        bquadratureOrder = 2 #don't need two here?
    boundaryQuadrature = SimplexGaussQuadrature(nd-1)
    boundaryQuadrature.setOrder(bquadratureOrder)
    for bintegral in OneLevelScalarTransport.elementBoundaryIntegralKeys:
        elementBoundaryQuadrature[bintegral] = boundaryQuadrature
    # # setup finite element spaces and stuff
    if verbose > 0:
        print 'definining finite element spaces'
    #end verbose
    #try P^1 nc ?
    if useCG:
        FemSpace = C0_AffineLinearOnSimplexWithNodalBasis(mesh,nd)
    else:
        FemSpace = NC_AffineLinearOnSimplexWithNodalBasis(mesh,nd)
    u   = FiniteElementFunction(FemSpace)
    phi = u
    if verbose > 0:
        print 'definining boundary conditions'
    #end verbose
    dirichletBCs=DOFBoundaryConditions(
        FemSpace,dirProb)

    # # setup numerical approximation configuration
    if verbose > 0:
        print 'definining linear algebra'
    #end verbose
    MatType = SparseMat
    matType = 'csr'
    linearSolverType=  'SparseLU'
    # # time integration
    #steady state
    if verbose > 0:
        print 'definining time integration'
    #end verbose
    TimeIntegrationClass = TimeIntegration
    # # flux approximations
    if verbose > 0:
        print 'definining flux approximations'
    #end verbose
    #diffusion only
    conservativeFlux = None
    numericalFlux    = None
    stabilization = None
    shockCapturing= None
    shockCapturingDiffusion= None

    # #create a single level solver
    if verbose > 0:
        print 'creating single level system'
    #end verbose
    system = OneLevelScalarTransport(u,
                                     phi,
                                     FemSpace,
                                     dirichletBCs,
                                     probCoefficients,
                                     quadrature,
                                     elementBoundaryQuadrature,
                                     stabilization,
                                     shockCapturing,
                                     shockCapturingDiffusion,
                                     conservativeFlux,
                                     numericalFlux,
                                     TimeIntegrationClass)
    #need this?
    if verbose > 0:
        print 'femSpace dim= ',FemSpace.dim
        print 'system dim = ',system.dim
    #end if
    y  = Numeric.zeros((system.dim,),Numeric.Float)
    dy = Numeric.zeros((system.dim,),Numeric.Float)
    r  = Numeric.zeros((system.dim,),Numeric.Float)
    system.updateQuadrature()
    system.setFreeDOF(y)
    system.updateCoefficients()
    system.updateQuadrature()
    system.matType = matType
    
    # # create nonlinear system for solving problem
    if verbose > 0:
        print 'creating nonlinear system'
    #end verbose
    
    jacobian = MatType(system.dim,system.dim,7)
    # # create linear system for solving problem
    if verbose > 5:
        print 'creating linear system'
        print 'initial jacobian mat = ',jacobian
    #end verbose

    linearSolver = SparseLU(jacobian)

    nlSolver = Newton(linearSolver,system,jacobian)
    #need intial residual calculated
    system.getResidual(y,r)
    if verbose > 5:
        print 'system y0= ',y
        system.getJacobian(jacobian)
        print 'system residual= ',r
        print 'system jacobian = ',jacobian
    # # solving problem
    if verbose > 0:
        print 'trying to solve nonlinear system'
    #end verbose
    nlSolver.solve(y,r)

        
    if verbose > 0:
        print 'nonlinear solver done'
    #end if

    #create some numerical quadrature points to evaluate errors. Allow for
    #these to be different order than the ones used in solution method

    errorQuadOrder  = 2
    errorQuadrature = SimplexGaussQuadrature(nd)
    errorQuadrature.setOrder(errorQuadOrder)

    errQuadX,errQuadW,errQuadRefX,errQuadRefW = \
      DiagUtils.getQuadraturePhysPointsAndWeights(mesh,FemSpace,errorQuadrature,verbose)
    nquad = Numeric.shape(errQuadX)[1]
    if verbose > 2:
        print 'error calc nquad= ',nquad
        print 'errQuadX.shape= ',Numeric.shape(errQuadX)
        print 'errQuadW.shape= ',Numeric.shape(errQuadW)
    #end if verbose
    if verbose > 3:
        print 'errQuadW = ',errQuadW
        print 'errQuadX = ',errQuadX
    #end if verbose

    uvals = DiagUtils.getFEMvals(u,errQuadRefX,verbose)
    
    fluxdofs = Numeric.zeros((mesh.nElements_global,nd+1),Numeric.Float)
    #should I make this always be 3d?
    qvals     = Numeric.zeros((mesh.nElements_global,nquad,nd),Numeric.Float)
    qexvals   = Numeric.zeros((mesh.nElements_global,nquad,nd),Numeric.Float)
    uAvgDofs  = Numeric.zeros(mesh.nElements_global,Numeric.Float)
    #need to put in correct higher order
    #quadrature and calculate average to higher order too
    uAvgVals  = Numeric.zeros((mesh.nElements_global,nquad),Numeric.Float)
    uexAvgVals= Numeric.zeros((mesh.nElements_global,nquad),Numeric.Float)
    if useCG:
        postProcessP1flux(u,AProb,fProb,nd,fluxdofs,verbose)
        postProcessP1avg(u,AProb,fProb,nd,fluxdofs,uAvgDofs,verbose)
    else:
        postProcessRT0flux(u,AProb,fProb,nd,fluxdofs,verbose)
        postProcessRT0potential(u,AProb,fProb,nd,fluxdofs,uAvgDofs,verbose)
    #end if

    
    for eN in range(mesh.nElements_global):
        uexAvgSum = 0.0
        volSum    = 0.0
        for k in range(nquad):
            qvals[eN,k,:]   = fluxdofs[eN,0:nd]+fluxdofs[eN,nd]*errQuadX[eN,k,0:nd]
            #mwf debug
            #print 'qex.uOfXT ',qexactSoln.uOfXT(errQuadX[eN,k,0:nd],0.0)
            #print 'qexvals ',qexvals[eN,k,:]
            if nd == 1:
                qexvals[eN,k,:] = qexactSoln.uOfXT(errQuadX[eN,k,0:nd],0.0)[0]
            else:
                qexvals[eN,k,:] = qexactSoln.uOfXT(errQuadX[eN,k,0:nd],0.0)
            #end if
            uexAvgSum += exactSoln.uOfXT(errQuadX[eN,k,0:nd],0.0)*errQuadW[eN,k]
            volSum    += errQuadW[eN,k]
        #end k
        for k in range(nquad):
            uAvgVals[eN,k]   = uAvgDofs[eN]
            uexAvgVals[eN,k] = uexAvgSum/volSum
        #end k again
    #end eN
    if verbose > 1:
        qvalb  = Numeric.zeros(nd,Numeric.Float)
        qexvalb= Numeric.zeros(nd,Numeric.Float)
        print 'fluxes at barycenters are '
        for eN in range(mesh.nElements_global):
            #mwf debug
            print 'x_T= ',mesh.elementList[eN].barycenter[0:nd]
            print 'qdofs= ',fluxdofs[eN,0:nd+1]
            qexvalb  = qexactSoln.uOfXT(mesh.elementList[eN].barycenter[0:nd],0.0)
            qvalb[:] = fluxdofs[eN,0:nd] \
                       + fluxdofs[eN,nd]*mesh.elementList[eN].barycenter[0:nd]
            print 'x= ',mesh.elementList[eN].barycenter[0:nd], \
                  ' q= ',qvalb,' qex= ',qexvalb
        #end eN

        print 'average potential values are '
        for eN in range(mesh.nElements_global):
            #mwf debug
            print 'x_T= ',mesh.elementList[eN].barycenter[0:nd]
            uexb = exactSoln.uOfXT(mesh.elementList[eN].barycenter[0:nd],0.0)
            print 'uRT0= ',uAvgVals[eN],' uex(x_T)= ',uexb, \
                  'uexAvg[eN,k]= ',uexAvgVals[eN,:]
        #end eN
    #end verbose
                              
    # # compute some errors
    eL2  = DiagUtils.L2errorSFEMvsAF(exactSoln,errQuadX,errQuadW,
                           uvals,T=0.0)
    print 'L2 error = ',eL2

    qeL2 = DiagUtils.L2errorFEMvsAF(qexactSoln,errQuadX,errQuadW,
                          qvals,T=0.0)
    print 'L2 error in flux = ',qeL2

    eaL2 = 0.0
    for eN in range(mesh.nElements_global):
        for k in range(nquad):
            eaL2 += ((uAvgVals[eN,k]-uexAvgVals[eN,k])**2)*errQuadW[eN,k]
        #end k
    #end eN
    eaL2 = sqrt(eaL2)
    
    print 'L2 error in averages = ',eaL2

    outmesh = """mesh%dd """ % nd
    if nd == 1:
        FemSpace.writeFunctionGnuplot(u,'solution')
    elif nd == 2:
        FemSpace.writeFunctionMatlab(u,'solution',outmesh)

    #go ahead and plot averages for matlab too? should work as before as long
    #as solution vector is same size as number of triangles
    if useCG:
        mesh.writeMeshEnsight(filename='solen')
        FemSpace.writeFunctionEnsight(u,filename='solen',append=False)
    
    # # try to look at solution # #
    import Gnuplot
    if useCG and nd == 2:
        solPlot = Gnuplot.Gnuplot()
        solPlot("set terminal x11")
        
        x = Numeric.arange(nx)/float(nx-1)
        y = Numeric.arange(ny)/float(ny-1)
        
        nSol = Numeric.reshape(system.u.dof,(nx,ny))
        solPlot('set parametric')
        solPlot('set data style lines')
        solPlot('set hidden')
        solPlot('set contour base')
        solPlot.xlabel('x')
        solPlot.ylabel('y')
        solPlot.splot(Gnuplot.GridData(nSol,x,y,binary=0,inline=0,
                                       title='numerical solution'))
    #end if
    if nd == 2:
        aPlot = Gnuplot.Gnuplot()
        aPlot("set terminal x11")
        
        x = Numeric.arange(nx)/float(nx-1)
        y = Numeric.arange(ny)/float(ny-1)
        aSol = Numeric.zeros((nx,ny),Numeric.Float)
        for i in range(nx):
            for j in range(ny):
                aSol[i,j] = exactSoln.uOfXT([x[i],y[j],0.0],0.0)
            #end j
        #end i
        aPlot('set parametric')
        aPlot('set data style lines')
        aPlot('set hidden')
        aPlot('set contour base')
        aPlot.xlabel('x')
        aPlot.ylabel('y')
        aPlot.splot(Gnuplot.GridData(aSol,x,y,binary=0,inline=0,title='analytical solution'))
    #end if 2d
    
    raw_input('Please press return to continue... \n')

## @}

if __name__ == '__main__':
    import sys
    from MeshTools import TriangularMesh
    import DiagUtils
    from optparse import OptionParser
    parser = OptionParser()
    #options controlling simulation behavior
    parser.add_option('-d','--spacedim',
                      default=2,
                      help="""number of spatial dimensions [2]""")
    parser.add_option('--useCG',
                      default=False,action='store_true',
                      help="""use C0 P1 approximation [False]""")
    parser.add_option('--nx',
                      default=11,
                      help="""number of nodes along x axis [11]""")
    parser.add_option('--ny',
                      default=11,
                      help="""number of nodes along y axis [11]""")

    parser.add_option('--nz',
                      default=11,
                      help="""number of nodes along z axis [11]""")

    parser.add_option('-t','--testNumber',
                      default=2,
                      help="""which test to use [2]
                      0 --- testCrRavNodalBasis
                      1 --- testEdgeDOFMap
                      2 --- testLaplacian
                      """)
    parser.add_option('-v','--verbose',
                      default=0,
                      help="""level of verbosity in simulation [0]""")


    #get options
    options, args = parser.parse_args(sys.argv[1:]) #get command line args
    #
    verbose = int(options.verbose)
    useCG   = options.useCG
    #number of space dimensions
    nd = int(options.spacedim)
    nx = int(options.nx)
    ny = int(options.ny)
    nz = int(options.nz)
    testParams = {}
    testParams['number'] = int(options.testNumber)
    #create mesh
    Lx = 1.0
    Ly = 1.0
    Lz = 1.0
    if nd == 1:
        grid1d = RectangularGrid(nx,1,1,Lx,Ly,Lz)
        mesh1d = EdgeMesh()
        mesh1d.rectangularToEdge(grid1d)
        mesh1d.computeGeometricInfo()
        testParams['mesh'] = mesh1d
    elif nd == 2:
        mesh2d = TriangularMesh()
        mesh2d.constructTriangularMeshOnRectangle(Lx,Ly,nx,ny)
        mesh2d.computeGeometricInfo()
        testParams['mesh'] = mesh2d
    else:
        grid3d = RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
        mesh3d = TetrahedralMesh()
        mesh3d.rectangularToTetrahedral(grid3d)
        mesh3d.computeGeometricInfo()
        testParams['mesh'] = mesh3d
    #end if on nd
    if testParams['number'] == 0:
        DiagUtils.testCrRavNodalBasis(nd,verbose)
    elif testParams['number'] == 1:
        DiagUtils.testEdgeDOFMap(testParams['mesh'],nd)
    else:
        testLaplacianP1nc(testParams['mesh'],Lx,Ly,Lz,nx,ny,nz,nd,useCG,verbose)
    #end else
