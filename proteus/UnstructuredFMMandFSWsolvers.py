"""
Fast marching and fast sweeping solvers

.. inheritance-diagram:: proteus.UnstructuredFMMandFSWsolvers
   :parts: 1
"""
import numpy
import math
import sys,atexit
import FemTools,MeshTools,EGeometry
import StupidHeap as SHeap

########################################################################
#solvers
########################################################################
class FMMEikonalSolver:
    """Encapsulate naive implementation of Fast Marching Methods on
unstructured grids for

    .. math::

      \|\grad T\| = 1/F

    :math:`T = 0` on :math:`\Gamma`

    1d local solver is standard upwind approximation

    2d local solver variations: acute triangulations version 1 or
    version 2 from Qian Zhang etal 07 obtuse triangulation not
    implemented

    3d local solver varitions: not fully checked

    For now, the input should be non-negative!

    """
#    TODO:
#       3D version needs to be tested more
    from proteus import cfmmfsw
    def __init__(self,mesh,dofMap,nSpace,localSolverType='QianEtalV2',frontInitType='magnitudeOnly',#'magnitudeOnly',
                 debugLevel=3):
        self.mesh   = mesh
        self.nSpace = nSpace
        self.orderApprox = 1
        self.debugLevel= debugLevel
        #reality check
        assert 1 <= self.nSpace and self.nSpace <= 3, "1d,2d, and 3d only right now"
        assert self.orderApprox == 1, "first order only for now"
        #default speeds for Eikonal equation
        import numpy
        self.unitNodalSpeeds = numpy.ones((self.mesh.nNodes_global,),'d')

        self.frontInitFlag = 1
        if frontInitType == 'magnitudeOnly':
            self.frontInitFlag = 0
        #could pass in frontInit type here
        self.csolver = FMMEikonalSolver.cfmmfsw.FMMEikonalSolver(self.nSpace,self.mesh.cmesh)
        #
        self.localPWLreconstruction = FMMEikonalSolver.cfmmfsw.localPWLreconstruction

    def solve(self,phi0,T,nodalSpeeds=None,zeroTol=1.0e-4,trialTol=1.0e-1,verbose=0):
        """
        Test first order fast marching method algorithm for eikonal equation

        \|\grad T \| = 1, \phi(\vec x) = 0, x \in \Gamma

        assuming \phi_0 describes initial location of interface Gamma and
        has reasonable values (absolute values) for T close to Gamma. Here
        T can be interpreted as the travel time from Gamma.

        Right now assumes global node numbers <--> global dofs but this can be
        fixed easily
        Input


        phi0: dof array from P1 C0 FiniteElementFunction holding initial condition

        T   : dof array from P1 C0 FiniteElementFunction for solution


        Output
        T(\vec x_n)    : travel time from initial front to node (\vec x_n)

        Internal data structures

        Status : status of nodal point (dictionary)
             -1 --> Far
              0 --> Trial
              1 --> Known
        Trial : nodal points adjacent to front tuples (index,val) stored in heap

        TODO
          have return flag
        """
        import numpy

        assert len(T) == len(phi0), "phi0 and T must be same dimensionality"
        assert len(T) == self.mesh.nNodes_global, "FemSpaces must be C0 P1"
        #mwf debug
        #import pdb
        #pdb.set_trace()
        failed = False
        if nodalSpeeds is None:
            speed = self.unitNodalSpeeds
        else:
            speed = nodalSpeeds
        assert len(speed) == self.mesh.nNodes_global, "nodalSpeed dim= %s must be %s " % (len(speed),self.mesh.nNodes_global)

        failed = self.csolver.solve(phi0,speed,T,zeroTol=zeroTol,trialTol=trialTol,
                                    initFlag=self.frontInitFlag,verbose=verbose)
        return bool(failed)
    #solve
#class


class FSWEikonalSolver:
    """Encapsulate naive implementation of Fast Marching Methods on unstructured grids
      for

    .. math::

        \|\grad T\| = 1/F

    :math:`T = 0` on :math:`\Gamma`

    1d local solver is standard upwind approximation

    2d local solver variations: acute triangulations version 1 or
    version 2 from Qian Zhang etal 07 obtuse triangulation not
    implemented

    3d local solver variations: not fully checked


    For now, the input should be non-negative!

    """
#    TODO:
#       3D version needs to be tested more
    from proteus import cfmmfsw

    def __init__(self,mesh,dofMap,nSpace,iterAtol=1.0e-8,iterRtol=0.0,maxIts=100,
                 localSolverType='QianEtalV2',frontInitType='magnitudeOnly',#frontInitType='magnitudeOnly',
                 refPoints=None,
                 orderApprox=1,LARGE=1.234e28,debugLevel=3):
        self.mesh   = mesh
        self.nSpace = nSpace
        self.iterAtol= iterAtol
        self.iterRtol= iterRtol
        self.maxIts = maxIts
        self.orderApprox = 1
        self.LARGE       = LARGE
        self.debugLevel  = debugLevel
        self.xRefOrderingPoints = refPoints
        self.nRefOrderingPoints = None
        if self.xRefOrderingPoints is not None:
            self.nRefOrderingPoints = len(self.xRefOrderingPoints)
        #reality check
        assert 1 <= self.nSpace and self.nSpace <= 3, "1d,2d, and 3d only right now"
        assert self.orderApprox == 1, "first order only for now"
        #default speeds for Eikonal equation
        import numpy
        self.unitNodalSpeeds = numpy.ones((self.mesh.nNodes_global,),'d')

        self.frontInitFlag = 1
        if frontInitType == 'magnitudeOnly':
            self.frontInitFlag = 0
        self.csolver = None

        if self.xRefOrderingPoints is None:
            self.csolver = FSWEikonalSolver.cfmmfsw.FSWEikonalSolver(self.nSpace,self.mesh.cmesh,
                                                                     atol=self.iterAtol,rtol=self.iterRtol,
                                                                     maxIts=self.maxIts,
                                                                     initFlag=self.frontInitFlag)
        else:
            self.csolver = FSWEikonalSolver.cfmmfsw.FSWEikonalSolver(self.nSpace,self.mesh.cmesh,
                                                                     atol=self.iterAtol,rtol=self.iterRtol,
                                                                     maxIts=self.maxIts,
                                                                     initFlag=self.frontInitFlag,
                                                                     nRefPoints=self.nRefOrderingPoints,
                                                                     refPoints=self.xRefOrderingPoints)
        #
        self.localPWLreconstruction = FSWEikonalSolver.cfmmfsw.localPWLreconstruction
    #end init

    def solve(self,phi0,T,nodalSpeeds=None,zeroTol=1.0e-4,trialTol=1.0e-1,verbose=0):
        """
        Test first order fast sweeping method algorithm for eikonal equation

        \|\grad T \| = 1, \phi(\vec x) = 0, x \in \Gamma

        assuming \phi_0 describes initial location of interface Gamma and
        has reasonable values (absolute values) for T close to Gamma. Here
        T can be interpreted as the travel time from Gamma.

        Right now assumes global node numbers <--> global dofs but this can be
        fixed easily
        Input


        phi0: dof array holding P1 C0 FiniteElementFunction holding initial condition

        T   : dof array holding P1 C0 FiniteElementFunction for solution


        Output
        T(\vec x_n)    : travel time from initial front to node (\vec x_n)

        Internal data structures

        Status : status of nodal point (dictionary)
              0 --> Not Known (Trial)
              1 --> Known

        Order  : ordering of points in domain using l_p metric from fixed reference points

        """
        import numpy
        from math import sqrt, fmod
        assert len(T) == len(phi0), "phi0 and T must be same dimensionality"
        assert len(T) == self.mesh.nNodes_global, "FemSpaces must be C0 P1"

        failed = False
        if nodalSpeeds is None:
            speed = self.unitNodalSpeeds
        else:
            speed = nodalSpeeds
        assert len(speed) == len(T), "nodalSpeed dim= %s must be %s " % (len(speed),len(T))


        failed = self.csolver.solve(phi0,speed,T,zeroTol=zeroTol,trialTol=trialTol,
                                    initFlag=self.frontInitFlag,verbose=1)#mwf hack
        return bool(failed)
    #solve
#class

########################################################################
#test codes
########################################################################
def unstructuredEx1d(initFunc,Lx,nx,method='FMM',verbose=0):
    """
    run a couple of redistancing examples in 1d: circle and two circles
    """
    import numpy

    mesh = MeshTools.EdgeMesh()
    mesh.generateEdgeMeshFromRectangularGrid(nx,Lx)

    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh)

    FemPhi0  = FemTools.FiniteElementFunction(femSpace,name="phi0")
    FemPhi0p = FemTools.FiniteElementFunction(femSpace,name="phi0p")
    FemPhi0m = FemTools.FiniteElementFunction(femSpace,name="phi0m")
    FemTp    = FemTools.FiniteElementFunction(femSpace,name="Tp")
    FemTm    = FemTools.FiniteElementFunction(femSpace,name="Tm")

    phi0 = FemPhi0.dof ; phi0p = FemPhi0p.dof ; phi0m = FemPhi0m.dof ;
    Tp   = FemTp.dof; Tm = FemTm.dof


    icout = open("phi0.dat",'w')

    #construct initial level set, short cut assuming dofs <--> node numbers
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]
        phi0[I] = initFunc(x)
        phi0p[I]= max(phi0[I],0.0)
        phi0m[I]= abs(min(phi0[I],0.0))
        icout.write("%g %g \n" % (x,phi0[I]))
    #
    failed = False
    if method == 'FSW':
        solver = FSWEikonalSolver(mesh,FemPhi0.femSpace.dofMap.l2g,1,iterAtol=1.0e-8,maxIts=100)
        print "calling FSWEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,FemTp.dof,verbose=verbose)
        print "back. calling FSWEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,FemTm.dof,verbose=verbose)
        print "back."
    else:
        solver = FMMEikonalSolver(mesh,FemPhi0.femSpace.dofMap.l2g,1)
        print "calling FMMEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,FemTp.dof,verbose=verbose)
        print "back. calling FMMEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,FemTm.dof,verbose=verbose)
        print "back."


    fout = open("T.dat",'w')
    phout= open("phi.dat",'w')
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]
        fout.write("%g %g \n" % (x,Tp[I]))
        phout.write("%g %g \n" % (x,Tp[I]-Tm[I]))

    icout.close()
    fout.close()
    phout.close()

def unstructuredEx2d(initFunc,Lx,Ly,nx,ny,method='FMM',verbose=0):
    """
    run a couple of redistancing examples in 2d:
    """
    import numpy

    mesh = MeshTools.TriangularMesh()
    mesh.generateTriangularMeshFromRectangularGrid(nx,ny,Lx,Ly)

    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh)

    FemPhi0  = FemTools.FiniteElementFunction(femSpace,name="phi0")
    FemPhi0p = FemTools.FiniteElementFunction(femSpace,name="phi0p")
    FemPhi0m = FemTools.FiniteElementFunction(femSpace,name="phi0m")
    FemTp    = FemTools.FiniteElementFunction(femSpace,name="Tp")
    FemTm    = FemTools.FiniteElementFunction(femSpace,name="Tm")

    phi0 = FemPhi0.dof ; phi0p = FemPhi0p.dof ; phi0m = FemPhi0m.dof ;
    Tp   = FemTp.dof; Tm = FemTm.dof


    icout = open("phi0.dat",'w')

    #construct initial level set, short cut assuming dofs <--> node numbers
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]
        phi0[I] = initFunc(x,y)
        phi0p[I]= max(phi0[I],0.0)
        phi0m[I]= abs(min(phi0[I],0.0))
        icout.write("%g %g %g \n" % (x,y,phi0[I]))
    #
    failed = False
    if method == 'FSW':
        #test different ref nodes (say just 3 in middle of domain?
        refNodes = numpy.array([[0.25,0.25,0.0],[0.5,0.5,0.0],[0.75,0.75,0.0]])
        #solver = FSWEikonalSolver(mesh,FemPhi0.femSpace.dofMap.l2g,2,iterAtol=1.0e-8,refPoints=refNodes,maxIts=100)
        solver = FSWEikonalSolver(mesh,FemPhi0.femSpace.dofMap.l2g,2,iterAtol=1.0e-8,maxIts=100)
        print "calling FSWEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,FemTp.dof,verbose=verbose)
        print "back. calling FSWEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,FemTm.dof,verbose=verbose)
        print "back."
    else:
        solver = FMMEikonalSolver(mesh,FemPhi0.femSpace.dofMap.l2g,2)
        print "calling FMMEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,FemTp.dof,verbose=verbose)
        print "back. calling FMMEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,FemTm.dof,verbose=verbose)
        print "back."
    #meth switch

    fout = open("T.dat",'w')
    phout= open("phi.dat",'w')
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]
        fout.write("%g %g %g \n" % (x,y,Tp[I]))
        phout.write("%g %g %g \n" % (x,y,Tp[I]-Tm[I]))

    icout.close()
    fout.close()
    phout.close()

def unstructuredEx3d(initFunc,Lx,Ly,Lz,nx,ny,nz,method='FMM',verbose=0):
    """
    run a redistancing example in 3d:
    """
    import numpy

    mesh = MeshTools.TetrahedralMesh()
    mesh.generateTetrahedralMeshFromRectangularGrid(nx,ny,nz,Lx,Ly,Lz)

    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh)

    FemPhi0  = FemTools.FiniteElementFunction(femSpace,name="phi0")
    FemPhi0p = FemTools.FiniteElementFunction(femSpace,name="phi0p")
    FemPhi0m = FemTools.FiniteElementFunction(femSpace,name="phi0m")
    FemTp    = FemTools.FiniteElementFunction(femSpace,name="Tp")
    FemTm    = FemTools.FiniteElementFunction(femSpace,name="Tm")

    phi0 = FemPhi0.dof ; phi0p = FemPhi0p.dof ; phi0m = FemPhi0m.dof ;
    Tp   = FemTp.dof; Tm = FemTm.dof


    icout = open("phi0.dat",'w')

    #construct initial level set, short cut assuming dofs <--> node numbers
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]; z=mesh.nodeArray[I,2]
        phi0[I] = initFunc(x,y,z)
        phi0p[I]= max(phi0[I],0.0)
        phi0m[I]= abs(min(phi0[I],0.0))
        icout.write("%g %g %g %g \n" % (x,y,z,phi0[I]))
    #
    failed = False
    if method == 'FSW':
        solver = FSWEikonalSolver(mesh,FemPhi0.femSpace.dofMap.l2g,3,iterAtol=1.0e-8,maxIts=100)
        print "calling FSWEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,FemTp.dof,verbose=verbose)
        print "back. calling FSWEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,FemTm.dof,verbose=verbose)
        print "back."
    else:
        solver = FMMEikonalSolver(mesh,FemPhi0.femSpace.dofMap.l2g,3)
        print "calling FMMEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,FemTp.dof,verbose=verbose)
        print "back. calling FMMEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,FemTm.dof,verbose=verbose)
        print "back."
    #method switch

    fout = open("T.dat",'w')
    phout= open("phi.dat",'w')
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]; z=mesh.nodeArray[I,2]
        fout.write("%g %g %g %g \n" % (x,y,z,Tp[I]))
        phout.write("%g %g %g %g \n" % (x,y,z,Tp[I]-Tm[I]))

    icout.close()
    fout.close()
    phout.close()



########################################################################
#try to test out 3d versions
def test3dLocalSolver(verbose=0):
    #try some simple configurations that I can back out soln for
    import numpy,math
    nNodes=4; nSpace=3;
    nodes = numpy.zeros((nNodes,nSpace),'d')
    #reference tet
    nodes[1,:]=[1.0,0.0,0.0]; nodes[2,:]=[0.0,1.0,0.0]; nodes[3,:]=[0.0,0.0,1.0]

    T = numpy.zeros((nNodes,),'d')

    sqrt3 = math.sqrt(3.)
    waveNormal = 1.0/sqrt3*numpy.array([-1.,1.,1.])
    eikSpeed=1.0

    eN = 0;
    #nodes with causal ordering
    N_A = 1; N_B=0; N_C=2; N_D=3
    #generic node numbering
    N = [0,1,2];
    T[N_A]=0; T[N_B]=sqrt3/3.0; T[N_C]=2.0*sqrt3/3.0

    print "calling qianZhangLocalSolver\n\t nodes=%s \n N=%s \n\t T=%s " % (nodes,N,T)
    T_D = qianZhangLocalSolver3d(eN,N_D,N[0],N[1],N[2],nodes,T,eikSpeed,verbose=verbose)
    print "T_D= %s " % T_D

def unstructuredEx1dInCpp(initFunc,Lx,nx,method='FMM',verbose=0):
    """
    run a couple of redistancing examples in 1d: circle and two circles
    use c++ interface
    """
    import numpy
    from proteus import cfmmfsw

    mesh = MeshTools.EdgeMesh()
    mesh.generateEdgeMeshFromRectangularGrid(nx,Lx)

    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh)

    FemPhi0  = FemTools.FiniteElementFunction(femSpace,name="phi0")
    FemPhi0p = FemTools.FiniteElementFunction(femSpace,name="phi0p")
    FemPhi0m = FemTools.FiniteElementFunction(femSpace,name="phi0m")
    FemTp    = FemTools.FiniteElementFunction(femSpace,name="Tp")
    FemTm    = FemTools.FiniteElementFunction(femSpace,name="Tm")

    phi0 = FemPhi0.dof ; phi0p = FemPhi0p.dof ; phi0m = FemPhi0m.dof ;
    Tp   = FemTp.dof; Tm = FemTm.dof


    icout = open("phi0.dat",'w')

    #construct initial level set, short cut assuming dofs <--> node numbers
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]
        phi0[I] = initFunc(x)
        phi0p[I]= max(phi0[I],0.0)
        phi0m[I]= abs(min(phi0[I],0.0))
        icout.write("%g %g \n" % (x,phi0[I]))
    #
    failed = False
    nd = 1
    nodalSpeeds = numpy.ones((mesh.nNodes_global,),'d')
    if method == 'FSW':
        solver = cfmmfsw.FSWEikonalSolver(nd,mesh.cmesh,atol=1.0e-8,rtol=1.0e-8,maxIts=100,
                                          initFlag=0)
        print "calling FSWEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,nodalSpeeds,FemTp.dof,initFlag=0,verbose=verbose)
        print "back. calling FSWEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,nodalSpeeds,FemTm.dof,initFlag=0,verbose=verbose)
        print "back."
    else:
        solver = cfmmfsw.FMMEikonalSolver(nd,mesh.cmesh)
        print "calling FMMEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,nodalSpeeds,FemTp.dof,initFlag=0,verbose=verbose)
        print "back. calling FMMEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,nodalSpeeds,FemTm.dof,initFlag=0,verbose=verbose)
        print "back."


    fout = open("T.dat",'w')
    phout= open("phi.dat",'w')
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]
        fout.write("%g %g \n" % (x,Tp[I]))
        phout.write("%g %g \n" % (x,Tp[I]-Tm[I]))

    icout.close()
    fout.close()
    phout.close()

def unstructuredEx2dInCpp(initFunc,Lx,Ly,nx,ny,method='FMM',verbose=0):
    """
    run a couple of redistancing examples in 2d:
    """
    import numpy
    from proteus import cfmmfsw

    mesh = MeshTools.TriangularMesh()
    mesh.generateTriangularMeshFromRectangularGrid(nx,ny,Lx,Ly)

    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh)

    FemPhi0  = FemTools.FiniteElementFunction(femSpace,name="phi0")
    FemPhi0p = FemTools.FiniteElementFunction(femSpace,name="phi0p")
    FemPhi0m = FemTools.FiniteElementFunction(femSpace,name="phi0m")
    FemTp    = FemTools.FiniteElementFunction(femSpace,name="Tp")
    FemTm    = FemTools.FiniteElementFunction(femSpace,name="Tm")

    phi0 = FemPhi0.dof ; phi0p = FemPhi0p.dof ; phi0m = FemPhi0m.dof ;
    Tp   = FemTp.dof; Tm = FemTm.dof


    icout = open("phi0.dat",'w')

    #construct initial level set, short cut assuming dofs <--> node numbers
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]
        phi0[I] = initFunc(x,y)
        phi0p[I]= max(phi0[I],0.0)
        phi0m[I]= abs(min(phi0[I],0.0))
        icout.write("%g %g %g \n" % (x,y,phi0[I]))
    #
    failed = False
    nd = 2
    nodalSpeeds = numpy.ones((mesh.nNodes_global,),'d')
    if method == 'FSW':
        solver = cfmmfsw.FSWEikonalSolver(nd,mesh.cmesh,atol=1.0e-8,rtol=1.0e-8,maxIts=100,
                                          initFlag=0)
        print "calling FSWEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,nodalSpeeds,FemTp.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
        print "back. calling FSWEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,nodalSpeeds,FemTm.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
        print "back. failed= %s" % failed
    else:
        solver = cfmmfsw.FMMEikonalSolver(nd,mesh.cmesh)
        print "calling FMMEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,nodalSpeeds,FemTp.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
        print "back. calling FMMEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,nodalSpeeds,FemTm.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
        print "back."
    #meth switch

    fout = open("T.dat",'w')
    phout= open("phi.dat",'w')
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]
        fout.write("%g %g %g \n" % (x,y,Tp[I]))
        phout.write("%g %g %g \n" % (x,y,Tp[I]-Tm[I]))

    icout.close()
    fout.close()
    phout.close()

def unstructuredEx3dinCpp(initFunc,Lx,Ly,Lz,nx,ny,nz,method='FMM',verbose=0):
    """
    run a redistancing example in 3d:
    """
    import numpy
    from proteus import cfmmfsw


    mesh = MeshTools.TetrahedralMesh()
    mesh.generateTetrahedralMeshFromRectangularGrid(nx,ny,nz,Lx,Ly,Lz)

    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh)

    FemPhi0  = FemTools.FiniteElementFunction(femSpace,name="phi0")
    FemPhi0p = FemTools.FiniteElementFunction(femSpace,name="phi0p")
    FemPhi0m = FemTools.FiniteElementFunction(femSpace,name="phi0m")
    FemTp    = FemTools.FiniteElementFunction(femSpace,name="Tp")
    FemTm    = FemTools.FiniteElementFunction(femSpace,name="Tm")

    phi0 = FemPhi0.dof ; phi0p = FemPhi0p.dof ; phi0m = FemPhi0m.dof ;
    Tp   = FemTp.dof; Tm = FemTm.dof


    icout = open("phi0.dat",'w')

    #construct initial level set, short cut assuming dofs <--> node numbers
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]; z=mesh.nodeArray[I,2]
        phi0[I] = initFunc(x,y,z)
        phi0p[I]= max(phi0[I],0.0)
        phi0m[I]= abs(min(phi0[I],0.0))
        icout.write("%g %g %g %g \n" % (x,y,z,phi0[I]))
    #
    failed = False
    nd = 3
    nodalSpeeds = numpy.ones((mesh.nNodes_global,),'d')
    if method == 'FSW':
        solver = cfmmfsw.FSWEikonalSolver(nd,mesh.cmesh,atol=1.0e-8,rtol=1.0e-8,maxIts=100,
                                          initFlag=0)
        print "calling FSWEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,nodalSpeeds,FemTp.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
        print "back. failed= %s  calling FSWEikonalSolver.solve for - ..."  % failed
        failed = solver.solve(FemPhi0m.dof,nodalSpeeds,FemTm.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
        print "back. failed= %s" % failed
    else:
        solver = cfmmfsw.FMMEikonalSolver(nd,mesh.cmesh)
        print "calling FMMEikonalSolver.solve for + ..."
        failed = solver.solve(FemPhi0p.dof,nodalSpeeds,FemTp.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
        print "back. calling FMMEikonalSolver.solve for - ..."
        failed = solver.solve(FemPhi0m.dof,nodalSpeeds,FemTm.dof,zeroTol=1.0e-4,trialTol=1.0e-1,
                              initFlag=0,verbose=verbose)
    print "back."
    #method switch

    fout = open("T.dat",'w')
    phout= open("phi.dat",'w')
    for I in range(mesh.nNodes_global):
        x = mesh.nodeArray[I,0]; y = mesh.nodeArray[I,1]; z=mesh.nodeArray[I,2]
        fout.write("%g %g %g %g \n" % (x,y,z,Tp[I]))
        phout.write("%g %g %g %g \n" % (x,y,z,Tp[I]-Tm[I]))

    icout.close()
    fout.close()
    phout.close()



if __name__ == "__main__":
    import math
    #method = 'FMM'
    method = 'FSW'
    dim = 2
    #now need for mpi
    from proteus import Comm
    comm = proteus.Comm.get()

    def circle1d(x):
        return (x-0.5)**2 - 0.2**2
    def twoCircle1d(x):
        return min((x-0.25)**2 - 0.1**2,(x-0.75)**2 - 0.1**2)
    #
    def circle2d(x,y):
        return (x-0.5)**2 + (y-0.5)**2 - 0.2**2
    def fourPetal(x,y):
        r0 = 0.25; a = 40; b = 4;
        tx = x-0.5; ty = y-0.5
        r = math.sqrt(tx**2 + ty**2); th = math.atan2(tx,ty)
        pr = 0.5*(r0 + math.cos(b*th)/(a*r0))
        return r**2 - pr**2
    def twoCircle2d(x,y):
        r0 = 0.15; r1 = 0.15; c0 = (0.25,0.25); c1=(0.75,0.75)
        d20= (x-c0[0])**2 + (y-c0[1])**2 - r0**2; d21 =  (x-c1[0])**2 + (y-c1[1])**2 - r1**2
        return min(d20,d21)
    #
    def sphere3d(x,y,z):
        return (x-0.5)**2 + (y-0.5)**2 + (z-0.5)**2 - 0.2**2
    def twoSphere3d(x,y,z):
        return min((x-0.25)**2 + (y-0.25)**2 + (z-0.25)**2 - 0.1**2,
                   (x-0.75)**2 + (y-0.75)**2 + (z-0.75)**2 - 0.1**2)

    Lx = 1.; Ly = 1.; Lz = 1.

    if dim == 2:
        nx=21; ny=21
        #testFunc= circle2d
        #testFunc= fourPetal
        testFunc= twoCircle2d
        unstructuredEx2d(testFunc,Lx,Ly,nx,ny,method=method,verbose=0)
        #unstructuredEx2dInCpp(testFunc,Lx,Ly,nx,ny,method=method,verbose=1)

    elif dim == 1:
        #nx=11
        #testFunc= circle1d
        nx=41
        testFunc= twoCircle1d
        unstructuredEx1d(testFunc,Lx,nx,method=method,verbose=9)
        #unstructuredEx1dInCpp(testFunc,Lx,nx,method=method,verbose=9)
    else:
        #test3dLocalSolver(verbose=10)
        nx=21; ny = 21; nz=21
        testFunc= sphere3d
        #testFunc= twoSphere3d
        #unstructuredEx3d(testFunc,Lx,Ly,Lz,nx,ny,nz,method=method,verbose=0)
        unstructuredEx3dinCpp(testFunc,Lx,Ly,Lz,nx,ny,nz,method=method,verbose=1)
