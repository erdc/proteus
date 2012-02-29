from math import *
import numpy
import Quadrature
import Profiling
log = Profiling.logEvent


import cellam,cfemIntegrals,ctracking

"""
Tools to implement an ELLAM approximation within Proteus

TODO

"""

class ELLAMdiscretization:
    """
    implements the integral approximations for ELLAM
    """
    def __init__(self,transport,options):
        self.transport = transport

        #mesh data structures needed for now to do physical space evaluation of shape functions
        self.elementBoundaryOuterNormalsArray = numpy.zeros((self.transport.mesh.nElements_global,self.transport.mesh.nElementBoundaries_element,self.transport.nSpace_global),'d')

        ###ellam specific options with defauls here
        self.particleTrackingType = None
        self.useBackwardTrackingForOldMass = False
        self.slumpingFlag = 0 #0 -- none, 1 -- Russell, Binning (1d any way), 2 -- Berzins (or something close), 3 -- FCT approach
        self.SSIPflag = 0 #use strategic spatial integration points


        #grab options from user if available
        #todo clean this up, add set from options for zeroTol etc
        assert 'particleTracking' in dir(options), "ELLAM requires particleTracking type to be set in n file"
        self.particleTrackingType = options.particleTracking

        if 'useBackwardTrackingForOldMass' in dir(options):
            self.useBackwardTrackingForOldMass = options.useBackwardTrackingForOldMass

        if 'slumpingFlag' in dir(options):
            self.slumpingFlag = options.slumpingFlag

        if 'SSIPflag' in dir(options):
            self.SSIPflag = options.SSIPflag

        ##determine algorithm behaviors based on user options
        self.needToBackTrackSolution = False
        if self.slumpingFlag in [2,3,4] or self.SSIPflag > 0:
            self.needToBackTrackSolution = True

        ###for tracking
        #quadrature points
        self.q_x_track = {}; self.q_t_track={}; self.q_t_depart={}; self.q_dt_track={}; self.q_flag_track={}; self.q_element_track={}
        for ci in range(self.transport.nc):
            self.q_x_track[ci] = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nQuadraturePoints_element,3),'d')
            self.q_t_track[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nQuadraturePoints_element),'d')
            self.q_t_depart[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nQuadraturePoints_element),'d')
            self.q_dt_track[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nQuadraturePoints_element),'d')
            self.q_flag_track[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nQuadraturePoints_element),'i')
            self.q_flag_track[ci].fill(-1)
            self.q_element_track[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nQuadraturePoints_element),'i')
        #interpolation points for solution values
        self.x_track_ip = {}; self.t_track_ip={}; self.t_depart_ip={}; self.u_track_ip={}; self.flag_track_ip={}; self.element_track_ip={}
        self.u_dof_track = {}; self.u_dof_track_tmp = {} ;
        for ci in range(self.transport.nc):
            self.x_track_ip[ci] = numpy.zeros((self.transport.mesh.nElements_global,self.transport.n_phi_ip_element[ci],3),'d')
            self.t_track_ip[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.n_phi_ip_element[ci]),'d')
            self.t_depart_ip[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.n_phi_ip_element[ci]),'d')
            self.u_track_ip[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.n_phi_ip_element[ci]),'d')
            self.flag_track_ip[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.n_phi_ip_element[ci]),'i')
            self.flag_track_ip[ci].fill(-1)
            self.element_track_ip[ci]   = numpy.zeros((self.transport.mesh.nElements_global,self.transport.n_phi_ip_element[ci]),'i')
            self.u_dof_track[ci] = numpy.copy(self.transport.u[ci].dof)
            self.u_dof_track_tmp[ci] = numpy.copy(self.transport.u[ci].dof)

            #manually insert value in phi_ip for now
            if not self.transport.phi_ip.has_key(('u',ci)):
                self.transport.phi_ip[('u',ci)] = numpy.zeros((self.transport.mesh.nElements_global,self.transport.n_phi_ip_element[ci]),'d')


        ##'global tracked points arrays' when have variable quadrature rules (e.g., SSIPs)
        self.gq_x_track_offsets=None; self.gq_x_track=None; self.gq_t_track=None; self.gq_t_depart=None; self.gq_dt_track=None;
        self.gq_flag_track=None; self.gq_element_track=None; self.gq_dV=None; self.gq=None; self.gq_last=None
        self.gq_x_depart=None; self.gq_element_depart=None;
        #not really needed except for evaluateSolutionAtTrackedPoints convention
        #todo get rid of
        self.gq_flag_depart=None;

        ##particle tracker setup
        self.particle_tracker = options.particleTracking(self.transport.mesh,self.transport.nSpace_global,
                                                         activeComponentList=range(self.transport.nc))

        self.particle_tracker.setFromOptions(options)
        self.particle_tracker.updateTransportInformation(self)
        self.zeroSolutionTol_track = {}
        if 'zeroSolutionTol_track' in dir(options):
            for ci in range(self.transport.nc):
                self.zeroSolutionTol_track[ci]=options.zeroSolutionTol_track[ci]
        else:
            for ci in range(self.transport.nc):
                self.zeroSolutionTol_track[ci]=1.0e-8


        #need to be able to evaluate solution at old and new time levels in some cases
        #could make this a shallow copy otherwise to save memory
        self.u_dof_last = {}
        for ci in range(self.transport.nc):
            self.u_dof_last[ci] = numpy.copy(self.transport.u[ci].dof)


        if self.useBackwardTrackingForOldMass:
            #need this to evaluate coefficients at tracked points
            self.q_backtrack = {}
            #deep_keys = [('u',0),('m',0),('dm',0,0),('f',0),('df',0,0),('a',0,0),('velocity',0)]
            deep_keys = set([('u',ci) for ci in range(self.transport.nc)])
            deep_keys |= set([('m',ci) for ci in self.transport.coefficients.mass.keys()])
            #no longer need to evaluate these if calling evaluateMassOnly
            #deep_keys |= set([('f',ci) for ci in self.transport.coefficients.advection.keys()])
            #deep_keys |= set([('velocity',ci) for ci in range(self.transport.nc)])
            for ci,cjDict in self.transport.coefficients.mass.iteritems():
                deep_keys |= set([('dm',ci,cj) for cj in cjDict.keys()])
            #for ci,cjDict in self.transport.coefficients.advection.iteritems():
            #    deep_keys |= set([('df',ci,cj) for cj in cjDict.keys()])
            #for ci,ckDict in self.transport.coefficients.diffusion.iteritems():
            #    deep_keys |= set([('a',ci,ck) for ck in ckDict.keys()])

            shallow_keys=set(['x'])
            for k in self.transport.q.keys():
                if k in deep_keys:
                    self.q_backtrack[k] = numpy.copy(self.transport.q[k])
                elif k in shallow_keys:
                    self.q_backtrack[k] = self.transport.q[k]

        else:
            self.q_backtrack=self.transport.q #could only grab shallow copies of specific keys

        #boundary point tracking (inflow approx)
        self.ebqe_x_track = {}; self.ebqe_t_track={}; self.ebqe_t_depart={};  self.ebqe_flag_track={}; self.ebqe_element_track={}
        for ci in range(self.transport.nc):
            self.ebqe_x_track[ci] = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
            self.ebqe_t_track[ci]   = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.ebqe_t_depart[ci]   = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.ebqe_flag_track[ci]   = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'i')
            self.ebqe_flag_track[ci].fill(-1)
            self.ebqe_element_track[ci]   = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'i')

        #keep track of when to update tracking
        self.needToTrackPoints = True;
        self.tForLastTrackingStep = 0.0

        #outflow boundary approximation via trapezoidal rule
        for ci in range(self.transport.nc):
            self.transport.ebqe[('outflow_flux',ci)] = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.transport.ebqe[('outflow_flux_last',ci)] = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'d')

            self.transport.ebqe[('u',ci)] = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.transport.ebqe[('grad(u)',ci)] = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary,self.transport.nSpace_global),'d')


            self.transport.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'i')
            self.transport.ebqe[('advectiveFlux_bc',ci)] = numpy.zeros((self.transport.mesh.nExteriorElementBoundaries_global,self.transport.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        #

        ##data structures for slumping
        self.rightHandSideForLimiting = {}; self.elementResidualTmp = {}; self.elementModifiedMassMatrixCorrection = {}; self.elementSlumpingParameter = {}
        for ci in range(self.transport.nc):
            self.rightHandSideForLimiting[ci]= numpy.zeros((self.transport.nFreeDOF_global[ci],),'d')
            self.elementResidualTmp[ci] = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nDOF_test_element[ci]),'d')
            self.elementModifiedMassMatrixCorrection[ci] = numpy.zeros((self.transport.mesh.nElements_global,self.transport.nDOF_test_element[ci],self.transport.nDOF_test_element[ci]),'d')
            self.elementSlumpingParameter[ci] = numpy.zeros((self.transport.mesh.nElements_global,),'d')

        if self.slumpingFlag in [3,4]: #FCT approach, assumes C0 P1 for now
            #TODO only works for 1 component right now!!
            assert ci == 0, "slumpingFlag == 3,4 only works for 1 component right now, fix rowptr colind info"
            assert self.transport.nFreeVDOF_global == self.transport.mesh.nNodes_global, "slumpingFlag == 3,4 only works for no hardwired dirichlet bcs"
            self.globalEdgeLimiter = {}; self.consistentMassMatrix= {}
            self.FCT_Rim = numpy.zeros((self.transport.mesh.nNodes_global,),'d')
            self.FCT_Rip = numpy.zeros((self.transport.mesh.nNodes_global,),'d')
            for ci in range(self.transport.nc):
                self.globalEdgeLimiter[ci]   = numpy.zeros((self.transport.mesh.nodeStarArray.shape[0]+self.transport.mesh.nNodes_global,),'d')
                self.consistentMassMatrix[ci]= numpy.copy(self.globalEdgeLimiter[ci])

        ###variables for debugging
        self.totalInflowFlux_cur = numpy.zeros((self.transport.nc,),'d'); self.totalOutflowFlux_cur = numpy.zeros((self.transport.nc,),'d')
        self.totalMassOld_cur = numpy.zeros((self.transport.nc,),'d'); self.totalMassNew_cur = numpy.zeros((self.transport.nc,),'d')
    ### routines for tracking
    def trackQuadraturePoints(self,q):
        """
        track quadrature points in q['x'] backward from t^{n+1} --> t^{n},
          loads
             x_track[0]      : location of point at end of tracking
             t_track[0]      : time tracking ended
             flag_track[0]   : -1  -- point in interior at tOut
                               -2  -- point exited domain somewhere in (tIn,tOut)
                               -3  -- did not track (e.g., v = 0 or u = 0)
             element_track[0]     : element containing point at end of tracking
        save time steps for domain in
             dt_track[0] = t^{n+1} - t_track[0]
        Then
        track quadrature points in q['x'] forward from t^n --> t^{n+1},
          save
             'x_track[0]     : location of point at end of tracking
             t_track[0]      : time tracking ended
             flag_track[0]   : -1  -- point in interior at tOut
                               -2  -- point exited domain somewhere in (tIn,tOut)
                               -3  -- did not track (e.g., v = 0 or u = 0)
             element_track[0]     : element containing point at end of tracking

        """
        import pdb
        timeToTrackPoints = (self.transport.timeIntegration.t > self.transport.timeIntegration.tLast + 1.0e-8 or
                             abs(self.tForLastTrackingStep-self.transport.timeIntegration.t) > 1.0e-8)

        #by default, tracking element quadrature points only (q array)
        x_depart = {}
        nPoints_track  = {}
        for ci in range(self.transport.nc):
            x_depart[ci] = q['x']
            nPoints_track[ci] = self.transport.mesh.nElements_global*self.transport.nQuadraturePoints_element

        def setupInitialElementLocations(ci,q_e):
            for k in range(q_e[ci].shape[1]):
                q_e[ci][:,k] = numpy.arange(self.transport.mesh.nElements_global,dtype='i')
        #todo need to allow skipping nonzero points with q or gq

        #first generate SSIPs if needed
        #todo this could be turned into a data member
        #0 -- not backtracked at all
        #1 -- backtracked only nonzero solution points
        #2 -- backtracked everything
        #mwf debug
        #import pdb
        #pdb.set_trace()
        solutionBackTrackedFlag = 0
        if self.needToTrackPoints and timeToTrackPoints and self.SSIPflag > 0:
            self.trackSolutionBackwards(skipPointsWithZeroSolution=True)
            self.generateSSIPs()
            solutionBackTrackedFlag = 1
            self.trackSSIPs()
        if self.needToTrackPoints and timeToTrackPoints:
            #mwf debug
            #pdb.set_trace()
            #update velocity fields for particle tracking
            for ci in range(self.transport.nc):
                self.particle_tracker.setTrackingVelocity(self.transport.coefficients.adjoint_velocity_dofs_last[ci],ci,
                                                          self.transport.coefficients.adjoint_velocity_times_last[ci],
                                                          timeLevel=0,
                                                          trackingVelocity_l2g=self.transport.coefficients.adjoint_velocity_l2g[ci])
                self.particle_tracker.setTrackingVelocity(self.transport.coefficients.adjoint_velocity_dofs[ci],ci,
                                                          self.transport.coefficients.adjoint_velocity_times[ci],
                                                          timeLevel=1)


                log(" LADRellam tracking integration points backward ci=%s" % ci,level=2)
                self.q_t_depart[ci].fill(self.transport.timeIntegration.t)
                #in desired output time, out actual time
                self.q_t_track[ci].fill(self.transport.timeIntegration.tLast)
                #try all points, now set to -1 to try, -3 to skip, 0 or greater if a node of the mesh
                self.q_flag_track[ci].fill(-1)
                #assign ownership of quadrature points to elements
                setupInitialElementLocations(ci,self.q_element_track)

            #todo make sure activeComponents set explicitly?
            #mwf debug just play with forwardTrack call, normally backward tracking
            self.particle_tracker.backwardTrack(self.q_t_depart,
                                                self.q_t_track,
                                                nPoints_track,
                                                x_depart,
                                                self.q_element_track,
                                                self.q_x_track,
                                                self.q_flag_track)


            #mwf debug
            #pdb.set_trace()
            for ci in range(self.transport.nc):
                self.q_dt_track[ci]  = numpy.copy(self.q_t_depart[ci])
                self.q_dt_track[ci] -= self.q_t_track[ci]

            if not self.useBackwardTrackingForOldMass:
                for ci in range(self.transport.nc):
                    log(" LADRellam tracking integration points forward ci=%s " % ci,level=2)
                    #forward
                    self.q_t_depart[ci].fill(self.transport.timeIntegration.tLast)
                    self.q_t_track[ci].fill(self.transport.timeIntegration.t)
                    #todo setup so can skip points with zero solution using q or gq, need to evaluate u at gq
                    #try all points, now set to -1 to try, -3 to skip, 0 or greater if a node of the mesh
                    self.q_flag_track[ci].fill(-1)
                    #assign ownership of quadrature points to elements
                    setupInitialElementLocations(ci,self.q_element_track)


                #todo make sure activeComponents set explicitly?
                self.particle_tracker.forwardTrack(self.q_t_depart,
                                                   self.q_t_track,
                                                   nPoints_track,
                                                   x_depart,
                                                   self.q_element_track,
                                                   self.q_x_track,
                                                   self.q_flag_track)


            if self.needToBackTrackSolution and solutionBackTrackedFlag < 1:
                self.trackSolutionBackwards(skipPointsWithZeroSolution=False)

            #end tracking interpolation points
            self.needToTrackPoints = False
            self.tForLastTrackingStep=self.transport.timeIntegration.t
            #mwf debug
            #pdb.set_trace()
        #end need to track integration points

    #
    def trackSolutionBackwards(self,skipPointsWithZeroSolution=False):
        """
        track interpolation points backwards to get an approximate solution at new time level
        """
        x_depart_ip = {}
        nPoints_track_ip  = {}
        #mwf debug
        #import pdb
        #pdb.set_trace()
        for ci in range(self.transport.nc):
            #todo switch this to characteristic velocity, need _last values!
            self.particle_tracker.setTrackingVelocity(self.transport.coefficients.adjoint_velocity_dofs_last[ci],ci,
                                                      self.transport.coefficients.adjoint_velocity_times_last[ci],
                                                      timeLevel=0,
                                                      trackingVelocity_l2g=self.transport.coefficients.adjoint_velocity_l2g[ci])
            self.particle_tracker.setTrackingVelocity(self.transport.coefficients.adjoint_velocity_dofs[ci],ci,
                                                      self.transport.coefficients.adjoint_velocity_times[ci],
                                                      timeLevel=1)


            log(" LADRellam tracking integration points backward ci=%s" % ci,level=2)
            self.t_depart_ip[ci].fill(self.transport.timeIntegration.t)
            #in desired output time, out actual time
            self.t_track_ip[ci].fill(self.transport.timeIntegration.tLast)
            #try all points, now set to -1 to try, -3 to skip, 0 or greater if a node of the mesh
            self.flag_track_ip[ci].fill(-1)

            for k in range(self.element_track_ip[ci].shape[1]):
                self.element_track_ip[ci][:,k] = numpy.arange(self.transport.mesh.nElements_global,dtype='i')

            x_depart_ip[ci] = self.transport.u[ci].femSpace.interpolationPoints
            nPoints_track_ip[ci] = self.transport.mesh.nElements_global*x_depart_ip[ci].shape[1]
            if skipPointsWithZeroSolution:
                #could use normal FemSpace machinery if all of the trial functions have been built etc
                cellam.evaluateSolutionAtTrackedPoints(self.transport.nSpace_global,
                                                       self.transport.nDOF_trial_element[ci],
                                                       nPoints_track_ip[ci],
                                                       self.transport.mesh.nElements_global,
                                                       self.transport.mesh.nNodes_global,
                                                       self.transport.mesh.nNodes_element,
                                                       self.transport.mesh.nElementBoundaries_element,
                                                       self.transport.mesh.nodeArray,
                                                       self.transport.mesh.elementNodesArray,
                                                       self.transport.mesh.elementNeighborsArray,
                                                       self.elementBoundaryOuterNormalsArray,
                                                       x_depart_ip[ci],
                                                       self.t_depart_ip[ci],
                                                       self.element_track_ip[ci],
                                                       self.flag_track_ip[ci],
                                                       self.transport.u[ci].femSpace.dofMap.l2g,
                                                       self.transport.u[ci].dof,
                                                       self.transport.phi_ip[('u',ci)])
                cellam.tagNegligibleIntegrationPoints(nPoints_track_ip[ci],
                                                      self.zeroSolutionTol_track[ci],
                                                      x_depart_ip[ci],
                                                      self.transport.phi_ip[('u',ci)],
                                                      self.flag_track_ip[ci])
                #mwf debug
                #import pdb
                #pdb.set_trace()
        #todo make sure activeComponents set explicitly?

        self.particle_tracker.backwardTrack(self.t_depart_ip,
                                            self.t_track_ip,
                                            nPoints_track_ip,
                                            x_depart_ip,
                                            self.element_track_ip,
                                            self.x_track_ip,
                                            self.flag_track_ip)

        for ci in range(self.transport.nc):
            #evaluate solution at tracked interpolation points using old time level dofs
            cellam.evaluateSolutionAtTrackedPoints(self.transport.nSpace_global,
                                                   self.transport.nDOF_trial_element[ci],
                                                   nPoints_track_ip[ci],
                                                   self.transport.mesh.nElements_global,
                                                   self.transport.mesh.nNodes_global,
                                                   self.transport.mesh.nNodes_element,
                                                   self.transport.mesh.nElementBoundaries_element,
                                                   self.transport.mesh.nodeArray,
                                                   self.transport.mesh.elementNodesArray,
                                                   self.transport.mesh.elementNeighborsArray,
                                                   self.elementBoundaryOuterNormalsArray,
                                                   self.x_track_ip[ci],
                                                   self.t_track_ip[ci],
                                                   self.element_track_ip[ci],
                                                   self.flag_track_ip[ci],
                                                   self.transport.u[ci].femSpace.dofMap.l2g,
                                                   self.u_dof_last[ci],#self.u_dof_lim_last[ci],#todo put this in time integration?
                                                   self.u_track_ip[ci])
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #use finite element machinery to project solution at new time level
            self.u_dof_track_tmp[ci][:] = self.transport.u[ci].dof
            self.transport.u[ci].projectFromInterpolationConditions(self.u_track_ip[ci])
            self.u_dof_track[ci][:] = self.transport.u[ci].dof
            self.transport.u[ci].dof[:] = self.u_dof_track_tmp[ci]
        #ci
    #def


    ###routines for implementing ellam integrals
    def updateElementResidual(self,elementResidual):
        """
        separate out LHS and RHS assembly for slumping
        """
        #setup points and time steps for tracking
        self.trackQuadraturePoints(self.transport.q)
        self.updateElementResidualRHS(elementResidual)
        self.updateElementResidualLHS(elementResidual)
        #mwf debug
        for ci in range(self.transport.nc):
            balance = self.totalMassNew_cur[ci] - self.totalMassOld_cur[ci] + self.totalInflowFlux_cur[ci] + self.totalOutflowFlux_cur[ci]
            print "ELLAMtools in updateElementResidual ci= %s newMassCur= %s ; oldMassCur= %s ; inflow= %s ; outflow = %s ; balance = %s " % (ci,self.totalMassNew_cur[ci],self.totalMassOld_cur[ci],self.totalInflowFlux_cur[ci],self.totalOutflowFlux_cur[ci],balance)
    def updateElementResidualRHS(self,elementResidual):
        """
        accumulate ELLAM approximations in element residual that go on right hand side (at least conceptually)

        """
        for ci in self.transport.coefficients.reaction.keys():
            #weight by time step size
            self.transport.q[('dt*w*dV_r',ci)][:] = self.transport.q[('w*dV_r',ci)]
            #todo need faster loop
            for j in range(self.transport.nDOF_trial_element[0]):
                self.transport.q[('dt*w*dV_r',ci)][:,:,j]   *= self.q_dt_track[ci]
            cfemIntegrals.updateReaction_weak(self.transport.q[('r',ci)],
                                              self.transport.q[('dt*w*dV_r',ci)],
                                              elementResidual[ci])

        #
        # (m^{n},w^{n+1})
        self.approximateOldMassIntegral(elementResidual)
        #inflow
        self.approximateInflowBoundaryIntegral(elementResidual)
        #outflow
        self.approximateOutflowBoundaryIntegral(elementResidual)
        if self.slumpingFlag == 1:
            for ci in range(self.transport.nc):
                self.elementResidualTmp[ci].fill(0.0)
                self.elementResidualTmp[ci] -= elementResidual[ci]

    def updateElementResidualLHS(self,elementResidual):
        """
        accumulate ELLAM approximations in element residual that go on left hand side (at least conceptually)

        """
        # (m^{n+1,w^{n+1}) + (\Delta t(x) a\grad u, grad w^{n+1}) + (\Delta t(x) r,w^{n+1})
        for ci,ckDict in self.transport.coefficients.diffusion.iteritems():
            for ck in ckDict.keys():
                #weight by time step size
                self.transport.q[('dt*grad(w)*dV_a',ck,ci)][:] = self.transport.q[('grad(w)*dV_a',ck,ci)]
                #todo need faster loop
                for j in range(self.transport.nDOF_trial_element[0]):
                    for I in range(self.transport.nSpace_global):
                        self.transport.q[('dt*grad(w)*dV_a',ck,ci)][:,:,j,I] *= self.q_dt_track[ci]
                if self.transport.sd:
                    cfemIntegrals.updateDiffusion_weak_sd(self.transport.coefficients.sdInfo[(ci,ck)][0],self.transport.coefficients.sdInfo[(ci,ck)][1],
                                                          self.transport.q[('a',ci,ck)],
                                                          self.transport.q[('grad(phi)',ck)],
                                                          self.transport.q[('dt*grad(w)*dV_a',ck,ci)],
                                                          elementResidual[ci])
                else:
                    cfemIntegrals.updateDiffusion_weak_lowmem(self.transport.q[('a',ci,ck)],
                                                              self.transport.q[('grad(phi)',ck)],
                                                              self.transport.q[('dt*grad(w)*dV_a',ck,ci)],
                                                              elementResidual[ci])



        if False and self.SSIPflag > 0 and self.gq_x_depart != None:#todo come up with a better way to handle uninitialized cases (first step)
            self.approximateNewMassIntegralUsingSSIPs(elementResidual)
        else:
            for ci in self.transport.coefficients.mass.keys():
                #note not dm/dt but just m
                #cfemIntegrals.updateMass_weak(self.transport.q[('m',ci)],
                #                              self.transport.q[('w*dV_m',ci)],
                #                              elementResidual[ci])
                #self.totalMassNew_cur[ci] = numpy.sum(self.transport.q[('m',ci)]*self.transport.q[('dV_u',ci)])
                self.totalMassNew_cur[ci] = cellam.updateNewMass_weak(self.transport.nSpace_global,
                                                                      self.transport.nDOF_test_element[ci],
                                                                      self.transport.mesh.nElements_global,
                                                                      self.transport.mesh.nNodes_global,
                                                                      self.transport.mesh.nNodes_element,
                                                                      self.transport.mesh.nElementBoundaries_element,
                                                                      self.transport.nQuadraturePoints_element,
                                                                      self.transport.mesh.nodeArray,
                                                                      self.transport.mesh.elementNodesArray,
                                                                      self.transport.mesh.elementNeighborsArray,
                                                                      self.elementBoundaryOuterNormalsArray,
                                                                      self.transport.q[('dV_u',ci)],#self.transport.q['dV'],
                                                                      self.transport.q['x'][ci],
                                                                      self.transport.u[ci].femSpace.dofMap.l2g,
                                                                      self.transport.q[('m',ci)],
                                                                      elementResidual[ci])
        #mwf debug
        #pdb.set_trace()
        if self.slumpingFlag == 1:
            for ci in range(self.transport.nc):
                #assemble right hand side vector
                self.rightHandSideForLimiting[ci].fill(0.)
                cfemIntegrals.updateGlobalResidualFromElementResidual(self.transport.offset[ci],
                                                                      self.transport.stride[ci],
                                                                      self.transport.l2g[ci]['nFreeDOF'],
                                                                      self.transport.l2g[ci]['freeLocal'],
                                                                      self.transport.l2g[ci]['freeGlobal'],
                                                                      self.elementResidualTmp[ci],
                                                                      self.rightHandSideForLimiting[ci]);
            #calculate element level lumping parameters and
            #subtract off element level mass correction from residual
            #mwf hack test what happens in 1d with a local slumping condition
            if self.transport.nSpace_global == 1:
                testLocalApproximation = False
                if testLocalApproximation:
                    #gives i-1 biased solution that has overshoot (i-) side and under shoot (i+) side
                    cellam.calculateSlumpedMassApproximation1d_local(self.transport.u[ci].femSpace.dofMap.l2g,
                                                                     self.transport.mesh.elementNeighborsArray,
                                                                     self.transport.u[ci].dof,self.transport.u[ci].dof,
                                                                     self.transport.q[('dm',ci,ci)],
                                                                     self.transport.q[('w',ci)],
                                                                     self.transport.q[('v',ci)],
                                                                     self.transport.q[('dV_u',ci)],
                                                                     self.rightHandSideForLimiting[ci],
                                                                     elementResidual[ci],
                                                                     self.elementSlumpingParameter[ci],
                                                                     self.elementModifiedMassMatrixCorrection[ci])
                else:
                    cellam.calculateSlumpedMassApproximation1d(self.transport.u[ci].femSpace.dofMap.l2g,
                                                               self.transport.mesh.elementNeighborsArray,
                                                               self.transport.u[ci].dof,self.transport.u[ci].dof,
                                                               self.transport.q[('dm',ci,ci)],
                                                               self.transport.q[('w',ci)],
                                                               self.transport.q[('v',ci)],
                                                               self.transport.q[('dV_u',ci)],
                                                               self.rightHandSideForLimiting[ci],
                                                               elementResidual[ci],
                                                               self.elementSlumpingParameter[ci],
                                                               self.elementModifiedMassMatrixCorrection[ci])

            elif self.transport.nSpace_global == 2:
                tryLocalUpwind = False
                if tryLocalUpwind:
                    cellam.calculateSlumpedMassApproximation2d_upwind(self.transport.mesh.nodeArray,
                                                                      self.transport.mesh.elementNodesArray,
                                                                      self.transport.mesh.elementNeighborsArray,
                                                                      self.transport.mesh.nodeStarOffsets,
                                                                      self.transport.mesh.nodeStarArray,
                                                                      self.elementBoundaryOuterNormalsArray,
                                                                      self.transport.u[ci].femSpace.dofMap.l2g,
                                                                      self.transport.u[ci].dof,self.transport.u[ci].dof,
                                                                      self.transport.q[('dm',ci,ci)],
                                                                      self.transport.q[('df',ci,ci)],
                                                                      self.transport.q[('w',ci)],
                                                                      self.transport.q[('v',ci)],
                                                                      self.transport.q[('dV_u',ci)],
                                                                      self.rightHandSideForLimiting[ci],
                                                                      elementResidual[ci],
                                                                      self.elementSlumpingParameter[ci],
                                                                      self.elementModifiedMassMatrixCorrection[ci])

                else:
                    #test adjusting local slumping criterion?
                    adjustFactor = 1.0#some overshoot, looks pretty good over long term? 1.0/2.0
                    cellam.calculateSlumpedMassApproximation2d(self.transport.u[ci].femSpace.dofMap.l2g,
                                                               self.transport.mesh.elementNeighborsArray,
                                                               self.transport.u[ci].dof,self.transport.u[ci].dof,
                                                               self.transport.q[('dm',ci,ci)],
                                                               self.transport.q[('w',ci)],
                                                               self.transport.q[('v',ci)],
                                                               self.transport.q[('dV_u',ci)],
                                                               self.rightHandSideForLimiting[ci],
                                                               elementResidual[ci],
                                                               self.elementSlumpingParameter[ci],
                                                               self.elementModifiedMassMatrixCorrection[ci],
                                                               adjustFactor)


        elif self.slumpingFlag == 2:
            #start by using current solution to do limiting, then try back tracking
            if self.transport.nSpace_global == 1:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                cellam.calculateBerzinsSlumpedMassApproximation1d(self.transport.u[ci].femSpace.dofMap.l2g,
                                                                  self.transport.mesh.elementNeighborsArray,
                                                                  self.transport.u[ci].dof,self.u_dof_track[ci],
                                                                  self.transport.q[('dm',ci,ci)],
                                                                  self.transport.q[('w',ci)],
                                                                  self.transport.q[('v',ci)],
                                                                  self.transport.q[('dV_u',ci)],
                                                                  self.rightHandSideForLimiting[ci],
                                                                  elementResidual[ci],
                                                                  self.elementModifiedMassMatrixCorrection[ci])
            elif self.transport.nSpace_global == 2:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                cellam.calculateBerzinsSlumpedMassApproximation2d(self.transport.u[ci].femSpace.dofMap.l2g,
                                                                  self.transport.mesh.elementNeighborsArray,
                                                                  self.transport.u[ci].dof,self.u_dof_track[ci],
                                                                  self.transport.q[('dm',ci,ci)],
                                                                  self.transport.q[('w',ci)],
                                                                  self.transport.q[('v',ci)],
                                                                  self.transport.q[('dV_u',ci)],
                                                                  self.rightHandSideForLimiting[ci],
                                                                  elementResidual[ci],
                                                                  self.elementModifiedMassMatrixCorrection[ci])
        elif self.slumpingFlag == 3:
            #TODO move this somewhere else? what if in parallel? just ignore off processor coupling
            #TODO only works for 1 component right now!!
            assert ci == 0, "slumpingFlag == 3 only works for 1 component right now, fix rowptr colind info"
            assert self.transport.nFreeVDOF_global == self.transport.mesh.nNodes_global, "slumpingFlag == 3 only works for no hardwired dirichlet bcs"
            #manualy assemble the global mass matrix
            self.consistentMassMatrix[ci].fill(0.0)
            cellam.manuallyUpdateGlobalMassMatrix(self.transport.rowptr,
                                                  self.transport.colind,
                                                  self.transport.u[ci].femSpace.dofMap.l2g,
                                                  self.transport.u[ci].dof,
                                                  self.transport.q[('dm',ci,ci)],
                                                  self.transport.q[('w',ci)],
                                                  self.transport.q[('v',ci)],
                                                  self.transport.q[('dV_u',ci)],
                                                  self.consistentMassMatrix[ci])

            #assumes C0P1
            cellam.computeSlumpingParametersFCT_KuzminMoeller10(self.transport.rowptr,
                                                                self.transport.colind,
                                                                self.transport.u[ci].dof,
                                                                self.u_dof_track[ci],
                                                                self.consistentMassMatrix[ci],
                                                                self.FCT_Rip,
                                                                self.FCT_Rim,
                                                                self.globalEdgeLimiter[ci])


            cellam.calculateElementSlumpedMassApproximationFromGlobalEdgeLimiter(self.transport.rowptr,
                                                                                 self.transport.colind,
                                                                                 self.transport.u[ci].femSpace.dofMap.l2g,
                                                                                 self.transport.u[ci].dof,
                                                                                 self.transport.q[('dm',ci,ci)],
                                                                                 self.transport.q[('w',ci)],
                                                                                 self.transport.q[('v',ci)],
                                                                                 self.transport.q[('dV_u',ci)],
                                                                                 self.globalEdgeLimiter[ci],
                                                                                 elementResidual[ci],
                                                                                 self.elementModifiedMassMatrixCorrection[ci])

        elif self.slumpingFlag == 4:
            #TODO move this somewhere else? what if in parallel? just ignore off processor coupling
            #TODO only works for 1 component right now!!
            assert ci == 0, "slumpingFlag == 4 only works for 1 component right now, fix rowptr colind info"
            assert self.transport.nFreeVDOF_global == self.transport.mesh.nNodes_global, "slumpingFlag == 3 only works for no hardwired dirichlet bcs"
            #manualy assemble the global mass matrix
            self.consistentMassMatrix[ci].fill(0.0)
            cellam.manuallyUpdateGlobalMassMatrix(self.transport.rowptr,
                                                  self.transport.colind,
                                                  self.transport.u[ci].femSpace.dofMap.l2g,
                                                  self.transport.u[ci].dof,
                                                  self.transport.q[('dm',ci,ci)],
                                                  self.transport.q[('w',ci)],
                                                  self.transport.q[('v',ci)],
                                                  self.transport.q[('dV_u',ci)],
                                                  self.consistentMassMatrix[ci])

            #assumes C0P1
            #mwf hack try using mass-lumping as low-order solution
            #requires 2 nonlinear iterations always, on the first call force lumping
            #also use current solution to do limiting
            #mwf debug
            #import pdb
            #pdb.set_trace()
            if self.transport.timeIntegration.low_order_step == True:
                #skip limiting and lump
                self.globalEdgeLimiter[ci].fill(0)
            elif self.transport.timeIntegration.low_order_step == False:
                cellam.computeSlumpingParametersFCT_KuzminMoeller10(self.transport.rowptr,
                                                                    self.transport.colind,
                                                                    self.transport.u[ci].dof,
                                                                    self.transport.timeIntegration.u_dof_low_order[ci],
                                                                    self.consistentMassMatrix[ci],
                                                                    self.FCT_Rip,
                                                                    self.FCT_Rim,
                                                                    self.globalEdgeLimiter[ci])

            cellam.calculateElementSlumpedMassApproximationFromGlobalEdgeLimiter(self.transport.rowptr,
                                                                                 self.transport.colind,
                                                                                 self.transport.u[ci].femSpace.dofMap.l2g,
                                                                                 self.transport.u[ci].dof,
                                                                                 self.transport.q[('dm',ci,ci)],
                                                                                 self.transport.q[('w',ci)],
                                                                                 self.transport.q[('v',ci)],
                                                                                 self.transport.q[('dV_u',ci)],
                                                                                 self.globalEdgeLimiter[ci],
                                                                                 elementResidual[ci],
                                                                                 self.elementModifiedMassMatrixCorrection[ci])


    #
    def updateElementJacobian(self,elementJacobian):
        """
        accumulate ELLAM approximations in element residual
        """
        for ci,ckDict in self.transport.coefficients.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                for cj in set(cjDict.keys()+self.transport.coefficients.potential[ck].keys()):
                    #assume dt weighting has been set already
                    if self.transport.sd:
                        cfemIntegrals.updateDiffusionJacobian_weak_sd(self.transport.coefficients.sdInfo[(ci,ck)][0],self.transport.coefficients.sdInfo[(ci,ck)][1],
                                                                      self.transport.phi[ck].femSpace.dofMap.l2g,
                                                                      self.transport.q[('a',ci,ck)],
                                                                      self.transport.q[('da',ci,ck,cj)],
                                                                      self.transport.q[('grad(phi)',ck)],
                                                                      self.transport.q[('dt*grad(w)*dV_a',ck,ci)],
                                                                      self.transport.dphi[(ck,cj)].dof,
                                                                      self.transport.q[('v',cj)],
                                                                      self.transport.q[('grad(v)',cj)],
                                                                      elementJacobian[ci][cj])
                    else:
                        cfemIntegrals.updateDiffusionJacobian_weak_lowmem(self.transport.phi[ck].femSpace.dofMap.l2g,
                                                                          self.transport.q[('a',ci,ck)],
                                                                          self.transport.q[('da',ci,ck,cj)],
                                                                          self.transport.q[('grad(phi)',ck)],
                                                                          self.transport.q[('dt*grad(w)*dV_a',ck,ci)],
                                                                          self.transport.dphi[(ck,cj)].dof,
                                                                          self.transport.q[('v',cj)],
                                                                          self.transport.q[('grad(v)',cj)],
                                                                          elementJacobian[ci][cj])
        for ci,cjDict in self.transport.coefficients.reaction.iteritems():
            for cj in cjDict:
                #assume dt weighting has been set already
                cfemIntegrals.updateReactionJacobian_weak_lowmem(self.transport.q[('dr',ci,cj)],
                                                                 self.transport.q[('v',cj)],
                                                                 self.transport.q[('dt*w*dV_r',ci)],
                                                                 elementJacobian[ci][cj])
        #todo handle Jacobian when using SSIPs even though shouldn't matter really for linear problem?
        for ci,cjDict in self.transport.coefficients.mass.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateMassJacobian_weak_lowmem(self.transport.q[('dm',ci,cj)],
                                                             self.transport.q[('v',cj)],
                                                             self.transport.q[('w*dV_m',ci)],
                                                             elementJacobian[ci][cj])

        #TODO unify so that all slumping approaches use same correction
        if self.slumpingFlag == 1:
            #mwf debug
            #import pdb
            #pdb.set_trace()
            useC = True
            for ci,cjDict in self.transport.coefficients.mass.iteritems():
                for cj in cjDict:
                    if useC:
                        cellam.updateElementJacobianWithSlumpedMassApproximation(self.elementSlumpingParameter[ci],
                                                                                 elementJacobian[ci][cj])
                    else:
                        for eN in range(self.transport.mesh.nElements_global):
                            for i in range(self.transport.nDOF_test_element[ci]):
                                self.elementJacobian[ci][cj][eN,i,i] += (self.transport.nDOF_trial_element[cj]-1)*self.elementSlumpingParameter[ci][eN]
                                for j in range(i):
                                    elementJacobian[ci][cj][eN,i,j] -= self.elementSlumpingParameter[ci][eN]
                                for j in range(i+1,self.transport.nDOF_trial_element[cj]):
                                    elementJacobian[ci][cj][eN,i,j] -= self.elementSlumpingParameter[ci][eN]

        elif self.slumpingFlag in [2,3,4]:
            #mwf debug
            #import pdb
            #pdb.set_trace()
            useC = True
            for ci,cjDict in self.transport.coefficients.mass.iteritems():
                for cj in cjDict:
                    if useC:
                        cellam.updateElementJacobianWithSlumpedMassCorrection(self.elementModifiedMassMatrixCorrection[ci],
                                                                              elementJacobian[ci][cj])
                    else:
                        for eN in range(self.transport.mesh.nElements_global):
                            for i in range(self.transport.nDOF_test_element[ci]):
                                for j in range(self.transport.nDOF_trial_element[cj]):
                                    elementJacobian[ci][cj][eN,i,j] += self.elementModifiedMassMatrixCorrection[ci][eN,i,j]


    def updateExteriorElementBoundaryJacobian(self,fluxJacobian_exterior):
        self.approximateOutflowBoundaryIntegralJacobian(fluxJacobian_exterior)

    def approximateOldMassIntegral(self,elementRes):
        """
        approximate weak integral
        \int_{\Omega} m^{n} w^{n+1} \dV
        """
        #by default, using just element quadrature array points (g)

        if self.useBackwardTrackingForOldMass:
            return self.approximateOldMassIntegralWithBackwardTracking(elementRes)
        if self.SSIPflag > 0 and self.gq_x_depart != None: #todo come up with a better way to avoid unitialized cases (first step)
            return self.approximateOldMassIntegralUsingSSIPs(elementRes)
        else:
            log("LADRellam evaluating old mass integral with q and forwardtracking",level=2)
            #mwf debug
            #import pdb
            #pdb.set_trace()
            for ci in range(self.transport.nc):
                self.totalMassOld_cur[ci] = cellam.updateOldMass_weak(self.transport.nSpace_global,
                                                                  self.transport.nDOF_test_element[ci],
                                                                  self.transport.mesh.nElements_global,
                                                                  self.transport.mesh.nNodes_global,
                                                                  self.transport.mesh.nNodes_element,
                                                                  self.transport.mesh.nElementBoundaries_element,
                                                                  self.transport.nQuadraturePoints_element,
                                                                  self.transport.mesh.nodeArray,
                                                                  self.transport.mesh.elementNodesArray,
                                                                  self.transport.mesh.elementNeighborsArray,
                                                                  self.elementBoundaryOuterNormalsArray,
                                                                  self.transport.q[('dV_u',ci)],#self.transport.q['dV'],
                                                                  self.q_x_track[ci],
                                                                  self.q_t_track[ci],
                                                                  self.q_element_track[ci],
                                                                  self.q_flag_track[ci],
                                                                  self.transport.u[ci].femSpace.dofMap.l2g,
                                                                  self.transport.timeIntegration.m_last[ci],
                                                                  elementRes[ci])

                # #mwf hack debug
                # dummyRes = numpy.copy(elementRes[ci]); dummyRes.fill(0.0);
                # totalMassOld_cur_mnew = cellam.updateOldMass_weak(self.transport.nSpace_global,
                #                                                   self.transport.nDOF_test_element[ci],
                #                                                   self.transport.mesh.nElements_global,
                #                                                   self.transport.mesh.nNodes_global,
                #                                                   self.transport.mesh.nNodes_element,
                #                                                   self.transport.mesh.nElementBoundaries_element,
                #                                                   self.transport.nQuadraturePoints_element,
                #                                                   self.transport.mesh.nodeArray,
                #                                                   self.transport.mesh.elementNodesArray,
                #                                                   self.transport.mesh.elementNeighborsArray,
                #                                                   self.elementBoundaryOuterNormalsArray,
                #                                                   self.transport.q[('dV_u',ci)],#self.transport.q['dV'],
                #                                                   self.q_x_track[ci],
                #                                                   self.q_t_track[ci],
                #                                                   self.q_element_track[ci],
                #                                                   self.q_flag_track[ci],
                #                                                   self.transport.u[ci].femSpace.dofMap.l2g,
                #                                                   self.transport.q[('m',ci)],
                #                                                   dummyRes)
                # print "UpdateOldMass totalMassOld_cur = %s with new solution = %s " % (self.totalMassOld_cur[ci],totalMassOld_cur_mnew)

    def approximateOldMassIntegralWithBackwardTracking(self,elementRes):
        """
        approximate weak integral
        \int_{\Omega} m^{n} w^{n+1} \dV using backward tracking
        """
        assert self.useBackwardTrackingForOldMass
        log("LADRellam evaluating old mass integral with backtracking",level=2)
        if self.SSIPflag > 0:
            assert False, "need to handle backtracking for old mass with SSIPs"
        #assumes that x_track, t_track etc correctly set from backtracking step in trackQuadraturePoints
        for ci in range(self.transport.nc):
            cellam.evaluateSolutionAtTrackedPoints(self.transport.nSpace_global,
                                                   self.transport.nDOF_trial_element[ci],
                                                   self.transport.mesh.nElements_global*self.transport.nQuadraturePoints_element,
                                                   self.transport.mesh.nElements_global,
                                                   self.transport.mesh.nNodes_global,
                                                   self.transport.mesh.nNodes_element,
                                                   self.transport.mesh.nElementBoundaries_element,
                                                   self.transport.mesh.nodeArray,
                                                   self.transport.mesh.elementNodesArray,
                                                   self.transport.mesh.elementNeighborsArray,
                                                   self.elementBoundaryOuterNormalsArray,
                                                   self.q_x_track[ci],
                                                   self.q_t_track[ci],
                                                   self.q_element_track[ci],
                                                   self.q_flag_track[ci],
                                                   self.transport.u[ci].femSpace.dofMap.l2g,
                                                   self.u_dof_last[ci],#todo put this in time integration?
                                                   self.q_backtrack[('u',ci)])

        #now evaluate as a standard mass integral
        #todo get rid of all of this, just want mass
        self.q_backtrack['dV']= self.transport.q['dV']
        #mwf debug
        #import pdb
        #pdb.set_trace()
        #if call full evaluate here, need to 'undo' to get velocity straight
        self.transport.coefficients.evaluateMassOnly(self.transport.timeIntegration.tLast,self.q_backtrack)

        for ci in range(self.transport.nc):
            #have to scale by -1
            self.q_backtrack[('m',ci)] *= -1.
            cfemIntegrals.updateMass_weak(self.q_backtrack[('m',ci)],
                                          self.transport.q[('w*dV_m',ci)],
                                          elementRes[ci])
            self.q_backtrack[('m',ci)] *= -1.
        #mwf debug
        #import pdb
        #pdb.set_trace()

    def approximateInflowBoundaryIntegral(self,elementRes):
        """
        approximate term

         \int_{t^n}^{t^{n+1}}  \int_{\Gamma_{I}\sigma^b w \dS \dt

        numerically using composite trapezoidal rule in time (and space too)

          \sum_{p=1}^{NT}\sum_{q=1}^{N_{q,b}}\Delta t^{p}\sigma^b(x_{q},t^p)w^{n+1}_{i}(\tilde{x}_q,t^{n+1})} W_q

        Here (x_q,t^p) tracks forward to  (\tilde{x}_q,t^{n+1}) and w^{n+1}_{i} is any test function with support
          covering (\tilde{x}_q,t^{n+1})

        only points on inflow boundary are tracked
        """
        if self.transport.timeIntegration.t > self.transport.timeIntegration.tLast + 1.0e-8:
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #update velocity fields for particle tracking
            ebqe_x_depart = {}
            ebqe_nPoints_track  = {}
            for ci in range(self.transport.nc):
                self.particle_tracker.setTrackingVelocity(self.transport.coefficients.adjoint_velocity_dofs_last[ci],ci,
                                                          self.transport.coefficients.adjoint_velocity_times_last[ci],
                                                          timeLevel=0,
                                                          trackingVelocity_l2g=self.transport.coefficients.adjoint_velocity_l2g[ci])
                self.particle_tracker.setTrackingVelocity(self.transport.coefficients.adjoint_velocity_dofs[ci],ci,
                                                          self.transport.coefficients.adjoint_velocity_times[ci],
                                                          timeLevel=1)
                ebqe_nPoints_track[ci]=self.transport.mesh.nExteriorElementBoundaries_global*self.transport.nElementBoundaryQuadraturePoints_elementBoundary
                ebqe_x_depart[ci] = self.transport.ebqe['x']
            self.NT = max(2,4*int(ceil(self.transport.timeIntegration.runCFL)))
            dtp = (self.transport.timeIntegration.t-self.transport.timeIntegration.tLast)/float(self.NT)
            integrationTimes = numpy.arange(self.NT+1,dtype='d')*dtp + self.transport.timeIntegration.tLast
            integrationTimeWeights=numpy.zeros(self.NT+1,'d'); integrationTimeWeights.fill(dtp)
            integrationTimeWeights[0] *= 0.5; integrationTimeWeights[-1] *= 0.5

            for tpi,dtpi in zip(integrationTimes,integrationTimeWeights):
                for ci in range(self.transport.nc):
                    #figure out which points on inflow need to be tracked
                    cellam.markInflowBoundaryPoints(self.transport.nSpace_global,
                                                    self.transport.timeIntegration.tLast,
                                                    self.transport.timeIntegration.t,
                                                    tpi,
                                                    self.transport.mesh.nExteriorElementBoundaries_global,
                                                    self.transport.nElementBoundaryQuadraturePoints_elementBoundary,
                                                    self.transport.mesh.exteriorElementBoundariesArray,
                                                    self.transport.mesh.elementBoundaryElementsArray,
                                                    self.transport.mesh.elementBoundaryLocalElementBoundariesArray,
                                                    self.transport.ebqe['x'],
                                                    self.transport.ebqe['n'],
                                                    self.transport.coefficients.ebqe[('velocity',ci)],#need to have time varying v
                                                    self.transport.coefficients.ebqe[('velocity',ci)],
                                                    self.transport.numericalFlux.isDOFBoundary[ci],
                                                    self.transport.ebqe[('advectiveFlux_bc_flag',ci)],
                                                    self.ebqe_element_track[ci],
                                                    self.ebqe_flag_track[ci])

                    #track forward
                    self.ebqe_t_depart[ci].fill(tpi)
                    self.ebqe_t_track[ci].fill(self.transport.timeIntegration.t)
                    #need to skip points with small boundary flux when tracking inflow boundary
                    skipPointsWithZeroSolution = 1

                    if skipPointsWithZeroSolution:
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()
                        cellam.tagNegligibleIntegrationPoints(ebqe_nPoints_track[ci],
                                                              self.zeroSolutionTol_track[ci],
                                                              ebqe_x_depart[ci],
                                                              self.transport.ebqe[('advectiveFlux_bc',ci)],
                                                              self.ebqe_flag_track[ci])

                direction = 1.0 #forward tracking
                if self.transport.timeIntegration.t > tpi + 1.0e-8:
                    self.particle_tracker.forwardTrack(self.ebqe_t_depart,
                                                       self.ebqe_t_track,
                                                       ebqe_nPoints_track,
                                                       ebqe_x_depart,
                                                       self.ebqe_element_track,
                                                       self.ebqe_x_track,
                                                       self.ebqe_flag_track)

                    for ci in range(self.transport.nc):
                        #accumulate into correct locations in residual
                        self.totalInflowFlux_cur[ci] = cellam.accumulateInflowFlux(self.transport.nSpace_global,
                                                                               self.transport.nDOF_test_element[ci],
                                                                               self.transport.mesh.nElements_global,
                                                                               self.transport.mesh.nNodes_global,
                                                                               self.transport.mesh.nNodes_element,
                                                                               self.transport.mesh.nElementBoundaries_element,
                                                                               self.transport.mesh.nExteriorElementBoundaries_global,
                                                                               self.transport.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                               self.transport.mesh.nodeArray,
                                                                               self.transport.mesh.elementNodesArray,
                                                                               self.transport.mesh.elementNeighborsArray,
                                                                               self.transport.mesh.exteriorElementBoundariesArray,
                                                                               self.transport.mesh.elementBoundaryElementsArray,
                                                                               self.transport.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                                    self.elementBoundaryOuterNormalsArray,
                                                                               tpi,
                                                                               dtpi,
                                                                               self.transport.ebqe['dS'],
                                                                               self.ebqe_x_track[ci],
                                                                               self.ebqe_t_track[ci],
                                                                               self.ebqe_element_track[ci],
                                                                               self.ebqe_flag_track[ci],
                                                                               self.transport.u[ci].femSpace.dofMap.l2g,
                                                                               self.transport.u[ci].dof,
                                                                               elementRes[ci],
                                                                               self.transport.coefficients.sdInfo[(ci,ci)][0], #todo fix
                                                                               self.transport.coefficients.sdInfo[(ci,ci)][1],
                                                                               self.transport.ebqe[('advectiveFlux_bc_flag',ci)],
                                                                               self.transport.ebqe[('advectiveFlux_bc',ci)])


    def approximateOutflowBoundaryIntegral(self,elementRes):
        """
        approximate term

         \int_{t^n}^{t^{n+1}}  \int_{\Gamma_{O}} f w \dS \dt

        numerically using trapezoidal rule in time


        """
        for ci in range(self.transport.nc):
            #accumulate into correct locations in residual
            self.totalOutflowFlux_cur[ci] = cellam.updateExteriorOutflowBoundaryFlux(self.transport.timeIntegration.t-self.transport.timeIntegration.tLast,
                                                                                 self.transport.nSpace_global,
                                                                                 self.transport.nDOF_test_element[ci],
                                                                                 self.transport.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                 self.transport.mesh.nExteriorElementBoundaries_global,
                                                                                 self.transport.mesh.exteriorElementBoundariesArray,
                                                                                 self.transport.mesh.elementBoundaryElementsArray,
                                                                                 self.transport.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                 self.transport.coefficients.ebqe[('velocity',ci)],
                                                                                 self.transport.ebqe[('n')],
                                                                                 self.transport.ebqe[('outflow_flux_last',ci)],
                                                                                 self.transport.ebqe[('w*dS_f',ci)],
                                                                                 self.transport.ebqe[('u',ci)],
                                                                                 self.transport.u[ci].femSpace.dofMap.l2g,
                                                                                 self.transport.ebqe[('outflow_flux',ci)],
                                                                                 elementRes[ci])


    #
    def approximateOutflowBoundaryIntegralJacobian(self,fluxJacobian_exterior):
        """
        approximate jacobian for term

         \int_{t^n}^{t^{n+1}}  \int_{\Gamma_{O}} f w \dS \dt

        numerically using trapezoidal rule in time


        """
        for ci in range(self.transport.nc):
            cellam.updateExteriorOutflowBoundaryFluxJacobian(self.transport.timeIntegration.t-self.transport.timeIntegration.tLast,
                                                             self.transport.nSpace_global,
                                                             self.transport.nDOF_test_element[ci],
                                                             self.transport.nDOF_trial_element[ci],
                                                             self.transport.nElementBoundaryQuadraturePoints_elementBoundary,
                                                             self.transport.mesh.nExteriorElementBoundaries_global,
                                                             self.transport.mesh.exteriorElementBoundariesArray,
                                                             self.transport.mesh.elementBoundaryElementsArray,
                                                             self.transport.mesh.elementBoundaryLocalElementBoundariesArray,
                                                             self.transport.coefficients.ebqe[('velocity',ci)],
                                                             self.transport.ebqe[('n')],
                                                             self.transport.ebqe[('outflow_flux_last',ci)],
                                                             self.transport.ebqe[('w*dS_f',ci)],
                                                             self.transport.ebqe[('u',ci)],
                                                             self.transport.ebqe[('v',ci)],
                                                             fluxJacobian_exterior[ci][ci])


    #
    ###time step management
    def updateTimeHistory(self,T,resetFromDOF=False):
        """
        todo find a better place to make sure know when a step is done
        because if step failes need to retrack
        """
        self.needToTrackPoints = True
        for ci in range(self.transport.nc):
            self.transport.ebqe[('outflow_flux_last',ci)].flat[:] = self.transport.ebqe[('outflow_flux',ci)].flat
        #todo put this in time integration
        #don't always need deep copy but go ahead and keep for now
        for ci in range(self.transport.nc):
            self.u_dof_last[ci].flat[:] = self.transport.u[ci].dof.flat

    def setInitialConditions(self,getInitialConditionsDict,T=0.0):
        #dont always need a deep copy but go ahead for now and keep
        for ci in range(self.transport.nc):
            self.u_dof_last[ci].flat[:] = self.transport.u[ci].dof.flat
            #go ahead and set tracking variable to
            self.u_dof_track[ci].flat[:] = self.transport.u[ci].dof.flat

    ###setting up geometry, quadrature point specific information
    def updateElementQuadrature(self,q):
        #extra boundary normal information for 2d, 3d to save need for ebq array
        boundaryNormals = numpy.array(self.transport.testSpace[0].elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
        ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                                q['inverse(J)'],
                                                self.elementBoundaryOuterNormalsArray)

        #mwf hack
        #TODO make sure coefficients has access to quadrature points for velocity evaluation??
        self.transport.coefficients.elementQuadraturePoints = self.transport.elementQuadraturePoints
    def calculateExteriorElementBoundaryQuadrature(self,ebqe):
        #mwf TODO cleanup make sure coefficients has access to quadrature points for velocity evaluation??
        self.transport.coefficients.elementBoundaryQuadraturePoints = self.transport.elementBoundaryQuadraturePoints

    ###SSIP routines
    def trackSSIPs(self):
        """
        track special integration points
        """
        x_depart = {}
        nPoints_track  = {}
        #todo get this loop out of python
        #def setupInitialElementLocations(ci,q_e):
        #    for eN in range(self.mesh.nElements_global):
        #        start = self.gq_x_track_offsets[ci][eN]; finish = self.gq_x_track_offsets[ci][eN+1]
        #        q_e[ci][start:finish] = eN

        #only forward track SSIPs
        assert not self.useBackwardTrackingForOldMass, "no need to use SSIPs with backtracking for mass"
        for ci in range(self.transport.nc):
            log(" LADRellam tracking SSIP  points forward ci=%s " % ci,level=2)
            nPoints_track[ci] = self.gq_x_track[ci].shape[0]

            self.gq_t_depart[ci].fill(self.transport.timeIntegration.tLast)
            self.gq_t_track[ci].fill(self.transport.timeIntegration.t)
            #todo setup so can skip points with zero solution using q or gq, need to evaluate u at gq
            #try all points, now set to -1 to try, -3 to skip, 0 or greater if a node of the mesh
            self.gq_flag_track[ci].fill(-1)
            #assign ownership of quadrature points to elements
            self.gq_element_track[ci][:] = self.gq_element_depart[ci]


        #todo make sure activeComponents set explicitly?
        self.particle_tracker.forwardTrack(self.gq_t_depart,
                                           self.gq_t_track,
                                           nPoints_track,
                                           self.gq_x_depart,
                                           self.gq_element_track,
                                           self.gq_x_track,
                                           self.gq_flag_track)


    def generateSSIPs(self):
        """
        The general idea is to track interpolation points back to old time level and create quadrature rules
        that include the backtracked images of the interpolation points
        ** assumes solution already backtracked **


        After tracking solution backwards, generate lookup table to determine points that are contained in each
        element (inverse of element_track). In this step, we should take care of the duplicate entries in the interpolation points.
        We also need to make sure that

        Then call dimension specific routines to create quadrature points
        on each element that contain the tracked points.

        """
        #have to redimension tracking arrays
        self.gq_x_track_offsets={}; self.gq_x_track={}; self.gq_t_track={}; self.gq_t_depart={}; self.gq_dt_track={}; self.gq_flag_track={}; self.gq_element_track={};
        self.gq_dV={}; self.gq={}; self.gq_last={}; self.gq_x_depart={}; self.gq_element_depart={}; self.gq_flag_depart={};
        #TODO make these user options
        #TODO make sub element quadrature type an option
        boundaryTolerance = 1.0e-6#1.0e-4;
        neighborTolerance = 1.0e-8#1.0e-4
        #mwf debug
        x_ssip = {}; x_ssip_offsets= {}
        useC = True
        if self.transport.nSpace_global == 2:
            #mwf debug
            #import pdb
            #pdb.set_trace()
            for ci in range(self.transport.nc):
                #determine which elements have SSIPs in them (remove duplicates and project to boundaries)
                x_ssip_offsets[ci],x_ssip[ci] = cellam.generateArraysForTrackedSSIPs(boundaryTolerance,
                                                                                     neighborTolerance,
                                                                                     self.transport.mesh.nodeArray,
                                                                                     self.transport.mesh.elementNodesArray,
                                                                                     self.transport.mesh.elementBoundariesArray,
                                                                                     self.elementBoundaryOuterNormalsArray,
                                                                                     self.transport.mesh.elementBoundaryBarycentersArray,
                                                                                     self.element_track_ip[ci],
                                                                                     self.flag_track_ip[ci],
                                                                                     self.x_track_ip[ci])

                #for debugging, loop through elements extract points and get back local quadrature points and weights
                import TriangleTools
                gq_dV_tmp = {}; gq_x_depart_tmp = {}; gq_element_depart = {}
                nPoints_global = 0
                for eN in range(self.transport.mesh.nElements_global):
                    if x_ssip_offsets[0][eN+1] > x_ssip_offsets[0][eN]:
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()

                        points = x_ssip[ci][x_ssip_offsets[0][eN]:x_ssip_offsets[0][eN+1]]
                        #the arrays are returned as nSubElement x nQuadraturePoints_subElement
                        gq_dV_tmp[eN],gq_x_depart_tmp[eN] = TriangleTools.testGenerateSSIPtriangulation(points)
                        nPoints_global += gq_dV_tmp[eN].shape[0]*gq_dV_tmp[eN].shape[1]
                    else:
                        nPoints_global += self.transport.q['dV'][eN].shape[0]
                #build actual arrays
                self.gq_element_depart[ci] = numpy.zeros((nPoints_global,),'i')
                self.gq_dV[ci]             = numpy.zeros((nPoints_global,),'d')
                self.gq_x_depart[ci]       = numpy.zeros((nPoints_global,3),'d')
                nSoFar = 0
                for eN in range(self.transport.mesh.nElements_global):
                    if gq_dV_tmp.has_key(eN):
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()
                        nPoints_eN = gq_dV_tmp[eN].shape[0]*gq_dV_tmp[eN].shape[1]
                        self.gq_dV[ci][nSoFar:nSoFar+nPoints_eN] = gq_dV_tmp[eN].flat[:]
                        self.gq_x_depart[ci][nSoFar:nSoFar+nPoints_eN].flat[:] = gq_x_depart_tmp[eN].flat[:]
                        self.gq_element_depart[ci][nSoFar:nSoFar+nPoints_eN] = eN
                        nSoFar += nPoints_eN
                    else: #copy over default quadrature
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()

                        nPoints_eN = self.transport.q['dV'][eN].shape[0]
                        self.gq_dV[ci][nSoFar:nSoFar+nPoints_eN] = self.transport.q['dV'][eN].flat[:]
                        self.gq_x_depart[ci][nSoFar:nSoFar+nPoints_eN].flat[:] = self.transport.q['x'][eN].flat[:]
                        self.gq_element_depart[ci][nSoFar:nSoFar+nPoints_eN] = eN
                        nSoFar += nPoints_eN

                #
                #generate other arrays that are needed
                #for now have to resize everthing here
                self.gq_x_track[ci]        = numpy.copy(self.gq_x_depart[ci])
                self.gq_t_track[ci]        = numpy.zeros((nPoints_global,),'d')
                self.gq_t_depart[ci]       = numpy.zeros((nPoints_global,),'d')
                self.gq_dt_track[ci]       = numpy.zeros((nPoints_global,),'d')
                self.gq_flag_track[ci]     = numpy.zeros((nPoints_global,),'i')
                self.gq_flag_depart[ci]     = numpy.zeros((nPoints_global,),'i')
                self.gq_element_track[ci]  = numpy.zeros((nPoints_global,),'i')
                self.gq[('u',ci)]          = numpy.zeros((nPoints_global,),'d')
                self.gq[('m',ci)]          = numpy.zeros((nPoints_global,),'d')
                self.gq_last[('u',ci)]     = numpy.zeros((nPoints_global,),'d')
                self.gq_last[('m',ci)]     = numpy.zeros((nPoints_global,),'d')
                for cj in self.transport.coefficients.mass[ci].keys():
                    self.gq[('dm',ci,cj)]      = numpy.zeros((nPoints_global,),'d')
                    self.gq_last[('dm',ci,cj)] = numpy.zeros((nPoints_global,),'d')

                self.gq[('x',ci)]          = self.gq_x_depart[ci] #simple alias for coeffficient evaluation
                self.gq_last[('x',ci)]     = self.gq_x_depart[ci] #simple alias for coeffficient evaluation
            #ci
        elif self.transport.nSpace_global == 1:
            if useC:
                for ci in range(self.transport.nc):
                    #mwf debug
                    #import pdb
                    #pdb.set_trace()
                    self.gq_element_depart[ci],self.gq_dV[ci],self.gq_x_depart[ci] = cellam.generateQuadratureArraysForSSIPs(boundaryTolerance,
                                                                                                                             neighborTolerance,
                                                                                                                   self.transport.mesh.nodeArray,
                                                                                                                   self.transport.mesh.elementNodesArray,
                                                                                                                   self.transport.mesh.elementBoundariesArray,
                                                                                                                   self.elementBoundaryOuterNormalsArray,
                                                                                                                   self.transport.mesh.elementBoundaryBarycentersArray,
                                                                                                                   self.element_track_ip[ci],
                                                                                                                   self.flag_track_ip[ci],
                                                                                                                   self.x_track_ip[ci],
                                                                                                                   self.transport.q['x'],
                                                                                                                   self.transport.q['dV'])

                    nPoints_global = self.gq_element_depart[ci].shape[0]

                    #for now have to resize everthing here
                    self.gq_x_track[ci]        = numpy.copy(self.gq_x_depart[ci])
                    self.gq_t_track[ci]        = numpy.zeros((nPoints_global,),'d')
                    self.gq_t_depart[ci]       = numpy.zeros((nPoints_global,),'d')
                    self.gq_dt_track[ci]       = numpy.zeros((nPoints_global,),'d')
                    self.gq_flag_track[ci]     = numpy.zeros((nPoints_global,),'i')
                    self.gq_flag_depart[ci]     = numpy.zeros((nPoints_global,),'i')
                    self.gq_element_track[ci]  = numpy.zeros((nPoints_global,),'i')
                    self.gq[('u',ci)]          = numpy.zeros((nPoints_global,),'d')
                    self.gq[('m',ci)]          = numpy.zeros((nPoints_global,),'d')
                    self.gq_last[('u',ci)]     = numpy.zeros((nPoints_global,),'d')
                    self.gq_last[('m',ci)]     = numpy.zeros((nPoints_global,),'d')
                    for cj in self.transport.coefficients.mass[ci].keys():
                        self.gq[('dm',ci,cj)]      = numpy.zeros((nPoints_global,),'d')
                        self.gq_last[('dm',ci,cj)] = numpy.zeros((nPoints_global,),'d')

                    self.gq[('x',ci)]          = self.gq_x_depart[ci] #simple alias for coeffficient evaluation
                    self.gq_last[('x',ci)]     = self.gq_x_depart[ci] #simple alias for coeffficient evaluation
                #ci
            else:
                #start by allocating memory on the fly and then make smarter
                #temporaries
                elementsToTrackedPoints = {}
                x_track_gq_offsets = {}
                x_track_gq         = {}
                dV_track_gq        = {}
                #todo allow for using only 1 component to determine SSIPs
                for ci in range(self.transport.nc):
                    elementsToTrackedPoints[ci] = {}
                    #mwf debug
                    #import pdb
                    #pdb.set_trace()
                    for k in range(len(self.element_track_ip[ci].flat)):
                        eN = self.element_track_ip[ci].flat[k]
                        if eN >= 0 and self.flag_track_ip[ci].flat[k] >= -1:
                            if elementsToTrackedPoints[ci].has_key(eN):
                                #todo: make sure only add points that are far enough away from existing points using a tolerance
                                elementsToTrackedPoints[ci][eN].add((self.x_track_ip[ci].flat[k*3+0],self.x_track_ip[ci].flat[k*3+1],self.x_track_ip[ci].flat[k*3+2]))
                            else:
                                #start with nodal points then add those that are tracked
                                elementsToTrackedPoints[ci][eN] = set([(self.transport.mesh.nodeArray[nN,0],self.transport.mesh.nodeArray[nN,1],self.transport.mesh.nodeArray[nN,2]) for nN in self.transport.mesh.elementNodesArray[eN]])
                                #todo: make sure only add points that are far enough away from existing points using a tolerance and
                                #      if the point is too close to a boundary, the project to the boundary and check that point is not too close
                                #      to an existing point
                                elementsToTrackedPoints[ci][eN] |= set([(self.x_track_ip[ci].flat[3*k+0],self.x_track_ip[ci].flat[3*k+1],self.x_track_ip[ci].flat[k*3+2])])
                    #
                    x_track_gq_offsets[ci] = numpy.zeros((self.transport.mesh.nElements_global+1,),'i')
                    #these will have to be converted to arrays
                    x_track_gq_tmp = {}; dV_track_gq_tmp = {}
                    if self.transport.nSpace_global == 1:
                        subQuadratureOrder = 2
                        subQuadratureType  = Quadrature.GaussEdge#Quadrature.CompositeTrapezoidalEdge#Quadrature.GaussEdge
                        #count number of points
                        for eN in range(self.transport.mesh.nElements_global):
                            if not elementsToTrackedPoints[ci].has_key(eN):
                                x_track_gq_offsets[ci][eN+1] = x_track_gq_offsets[ci][eN]+len(self.transport.q['dV'][eN])
                                #copy over q's integration points and weights to temporary data structures
                                dV_track_gq_tmp[eN] = numpy.copy(self.transport.q['dV'][eN])
                                x_track_gq_tmp[eN]  = numpy.copy(self.transport.q['x'][eN])
                            else:
                                #options are to generate quadrature physical directly or map back to reference
                                #mwf debug
                                #import pdb
                                #pdb.set_trace()
                                #subdivide element according to SSIPs then generate
                                #Gaussian quadrature on each sub-interval
                                #do manipulations in physical space first since that's
                                #how triangle would handle it I believe
                                #manually grab the points, sort, and subdivide
                                #generate a triangulation of element
                                tmpEdgeMesh = sorted(elementsToTrackedPoints[ci][eN])
                                #number of elements in sub-triangulation
                                nElements_base= len(tmpEdgeMesh)-1
                                subElementQuadrature = subQuadratureType()
                                subElementQuadrature.setOrder(subQuadratureOrder)
                                nSubElementPoints = len(subElementQuadrature.points)
                                nQuadraturePointsNew = nElements_base*nSubElementPoints
                                x_track_gq_offsets[ci][eN+1] = x_track_gq_offsets[ci][eN]+nQuadraturePointsNew
                                dV_track_gq_tmp[eN] = numpy.zeros((nQuadraturePointsNew,),'d')
                                x_track_gq_tmp[eN]  = numpy.zeros((nQuadraturePointsNew,3),'d')
                                #loop through each 'base' element in sub element triangulation and
                                #allocate the quadrature points and weights from the quadrature rule
                                #short-cut that may or may not be ok is to generate affine mapping on the fly
                                np_last = 0
                                for eN_local in range(nElements_base):
                                    d = numpy.zeros((3,),'d')
                                    for I in range(3):
                                        d[I]=tmpEdgeMesh[eN_local+1][I]-tmpEdgeMesh[eN_local][I]
                                    volume = numpy.sqrt(numpy.dot(d,d))
                                    for p,w in zip(subElementQuadrature.points,subElementQuadrature.weights):
                                        for I in range(3):
                                            x_track_gq_tmp[eN][np_last,I] = tmpEdgeMesh[eN_local][I]*(1.0-p[0]) + tmpEdgeMesh[eN_local+1][I]*p[0]
                                        dV_track_gq_tmp[eN][np_last]  = w*volume
                                        np_last += 1
                            #else has tracked points
                        #eN
                        nPoints_global = x_track_gq_offsets[ci][-1]
                        self.gq_x_track[ci] = numpy.zeros((nPoints_global,3),'d')
                        self.gq_dV[ci]= numpy.zeros((nPoints_global,),'d')
                        for eN in range(self.transport.mesh.nElements_global):
                            self.gq_x_track[ci][x_track_gq_offsets[ci][eN]:x_track_gq_offsets[ci][eN+1],:] =x_track_gq_tmp[eN][:,:]
                            self.gq_dV[ci][x_track_gq_offsets[ci][eN]:x_track_gq_offsets[ci][eN+1]]=dV_track_gq_tmp[eN][:]
                        #
                        self.gq_x_track_offsets[ci]= numpy.copy(x_track_gq_offsets[ci])
                        self.gq_x_depart[ci]       = numpy.copy(self.gq_x_track[ci])
                        #for now have to resize everthing here
                        self.gq_t_track[ci]        = numpy.zeros((nPoints_global,),'d')
                        self.gq_t_depart[ci]       = numpy.zeros((nPoints_global,),'d')
                        self.gq_dt_track[ci]       = numpy.zeros((nPoints_global,),'d')
                        self.gq_flag_track[ci]     = numpy.zeros((nPoints_global,),'i')
                        self.gq_flag_depart[ci]     = numpy.zeros((nPoints_global,),'i')
                        self.gq_element_track[ci]  = numpy.zeros((nPoints_global,),'i')
                        self.gq_element_depart[ci]  = numpy.zeros((nPoints_global,),'i')
                        self.gq[('u',ci)]          = numpy.zeros((nPoints_global,),'d')
                        self.gq[('m',ci)]          = numpy.zeros((nPoints_global,),'d')
                        self.gq_last[('u',ci)]     = numpy.zeros((nPoints_global,),'d')
                        self.gq_last[('m',ci)]     = numpy.zeros((nPoints_global,),'d')
                        for cj in self.transport.coefficients.mass[ci].keys():
                            self.gq[('dm',ci,cj)]      = numpy.zeros((nPoints_global,),'d')
                            self.gq_last[('dm',ci,cj)] = numpy.zeros((nPoints_global,),'d')

                        self.gq[('x',ci)]          = self.gq_x_depart[ci] #simple alias for coeffficient evaluation
                        self.gq_last[('x',ci)]     = self.gq_x_depart[ci] #simple alias for coeffficient evaluation
                        #go ahead and assign element_depart
                        for eN in range(self.transport.mesh.nElements_global):
                            start = self.gq_x_track_offsets[ci][eN]; finish = self.gq_x_track_offsets[ci][eN+1]
                            self.gq_element_depart[ci][start:finish] = eN

                    #1d
                    #mwf debug
                    #import pdb
                    #pdb.set_trace()
                #ci loop for generating SSIPs
            #not useC
        #1d
        #todo what about allowing x to be consistent with usual approach
        self.gq['x']          = self.gq_x_depart[0] #simple alias for coeffficient evaluation
        self.gq_last['x']     = self.gq_x_depart[0] #simple alias for coeffficient evaluation

        #mwf debug
        #print "generateSSIPs t= %g useC= %g sum(self.gq_dV[0].flat)= %g " % (self.transport.timeIntegration.t,useC,sum(self.gq_dV[0]))
        #print "eN el_track_ip[0] flag track[0] x_track_ip[0]"
        #for eN in range(self.x_track_ip[0].shape[0]):
        #    print "%d %s %s %s " % (eN,self.element_track_ip[0][eN],self.flag_track_ip[0][eN],self.x_track_ip[0][eN])
        #print "i x dV ele"
        #for i in range(self.gq_x_depart[0].shape[0]):
        #    print "%g %g %g %g" % (i,self.gq_x_depart[0][i,0],self.gq_dV[0][i],self.gq_element_depart[0][i])
        #

        for ci in range(self.transport.nc):
            self.gq_flag_depart[ci].fill(-1); self.gq_t_depart[ci].fill(self.transport.timeIntegration.tLast)

            cellam.evaluateSolutionAtTrackedPoints(self.transport.nSpace_global,
                                                   self.transport.nDOF_trial_element[ci],
                                                   self.gq_x_depart[ci].shape[0],
                                                   self.transport.mesh.nElements_global,
                                                   self.transport.mesh.nNodes_global,
                                                   self.transport.mesh.nNodes_element,
                                                   self.transport.mesh.nElementBoundaries_element,
                                                   self.transport.mesh.nodeArray,
                                                   self.transport.mesh.elementNodesArray,
                                                   self.transport.mesh.elementNeighborsArray,
                                                   self.elementBoundaryOuterNormalsArray,
                                                   self.gq_x_depart[ci],
                                                   self.gq_t_depart[ci],
                                                   self.gq_element_depart[ci],
                                                   self.gq_flag_depart[ci],
                                                   self.transport.u[ci].femSpace.dofMap.l2g,
                                                   self.u_dof_last[ci],#todo put this in time integration?
                                                   self.gq_last[('u',ci)])

        self.transport.coefficients.evaluateMassOnly(self.transport.timeIntegration.tLast,self.gq_last)
        #mwf debug
        #import pdb
        #pdb.set_trace()


    #
    def approximateNewMassIntegralUsingSSIPs(self,elementRes):
        """
        approximate weak integral
        \int_{\Omega} m^{n+1} w^{n+1} \dV using variable quadrature based on SSIPs
        """
        log("LADRellam evaluating new mass integral with SSIPs and forwardtracking",level=2)

        assert not self.useBackwardTrackingForOldMass, "no need to use SSIPs and backward tracking for mass"

        #evaluate time solution at SSIPs
        #mwf debug
        #import pdb
        #pdb.set_trace()
        #def setupInitialElementLocations(ci,q_e):
        #    for eN in range(self.transport.mesh.nElements_global):
        #        start = self.gq_x_track_offsets[ci][eN]; finish = self.gq_x_track_offsets[ci][eN+1]
        #        q_e[ci][start:finish] = eN

        #have to evaluate new time solution and mass at SSIPs
        for ci in range(self.transport.nc):
            self.gq_t_depart[ci].fill(self.transport.timeIntegration.t)
            self.gq_flag_depart[ci].fill(-1)
            #should be already set setupInitialElementLocations(ci,self.gq_element_depart)

            cellam.evaluateSolutionAtTrackedPoints(self.transport.nSpace_global,
                                                   self.transport.nDOF_trial_element[ci],
                                                   self.gq_x_depart[ci].shape[0],
                                                   self.transport.mesh.nElements_global,
                                                   self.transport.mesh.nNodes_global,
                                                   self.transport.mesh.nNodes_element,
                                                   self.transport.mesh.nElementBoundaries_element,
                                                   self.transport.mesh.nodeArray,
                                                   self.transport.mesh.elementNodesArray,
                                                   self.transport.mesh.elementNeighborsArray,
                                                   self.elementBoundaryOuterNormalsArray,
                                                   self.gq_x_depart[ci],
                                                   self.gq_t_depart[ci],
                                                   self.gq_element_depart[ci],
                                                   self.gq_flag_depart[ci],
                                                   self.transport.u[ci].femSpace.dofMap.l2g,
                                                   self.transport.u[ci].dof,
                                                   self.gq[('u',ci)])

        self.transport.coefficients.evaluateMassOnly(self.transport.timeIntegration.t,self.gq)

        for ci in range(self.transport.nc):
            #todo do away with scaling
            #for now just use old mass routine
            self.gq[('m',ci)] *= -1.0
            cellam.updateOldMass_weak_arbitraryQuadrature(self.transport.nSpace_global,
                                                          self.transport.nDOF_test_element[ci],
                                                          self.transport.mesh.nElements_global,
                                                          self.transport.mesh.nNodes_global,
                                                          self.transport.mesh.nNodes_element,
                                                          self.transport.mesh.nElementBoundaries_element,
                                                          self.gq_x_depart[ci].shape[0],
                                                          self.transport.mesh.nodeArray,
                                                          self.transport.mesh.elementNodesArray,
                                                          self.transport.mesh.elementNeighborsArray,
                                                          self.elementBoundaryOuterNormalsArray,
                                                          self.gq_dV[ci],
                                                          self.gq_x_depart[ci],
                                                          self.gq_t_depart[ci],
                                                          self.gq_element_depart[ci],
                                                          self.gq_flag_depart[ci],
                                                          self.transport.u[ci].femSpace.dofMap.l2g,
                                                          self.gq[('m',ci)],
                                                          elementRes[ci])
            self.gq[('m',ci)] *= -1.0

    def approximateOldMassIntegralUsingSSIPs(self,elementRes):
        """
        approximate weak integral
        \int_{\Omega} m^{n} w^{n+1} \dV using variable quadrature based on SSIPs
        """
        #todo need to figure out how to handle case like initial step where points
        #may not be tracked backwards yet

        #by default, using just element quadrature array points (g)
        log("LADRellam evaluating old mass integral with SSIPs and forwardtracking",level=2)

        assert not self.useBackwardTrackingForOldMass, "no need to use SSIPs and backward tracking for mass"

        #assume old time solution already evaluated at SSIPs
        #mwf debug
        #import pdb
        #pdb.set_trace()
        for ci in range(self.transport.nc):
            cellam.updateOldMass_weak_arbitraryQuadrature(self.transport.nSpace_global,
                                                          self.transport.nDOF_test_element[ci],
                                                          self.transport.mesh.nElements_global,
                                                          self.transport.mesh.nNodes_global,
                                                          self.transport.mesh.nNodes_element,
                                                          self.transport.mesh.nElementBoundaries_element,
                                                          self.gq_x_track[ci].shape[0],
                                                          self.transport.mesh.nodeArray,
                                                          self.transport.mesh.elementNodesArray,
                                                          self.transport.mesh.elementNeighborsArray,
                                                          self.elementBoundaryOuterNormalsArray,
                                                          self.gq_dV[ci],
                                                          self.gq_x_track[ci],
                                                          self.gq_t_track[ci],
                                                          self.gq_element_track[ci],
                                                          self.gq_flag_track[ci],
                                                          self.transport.u[ci].femSpace.dofMap.l2g,
                                                          self.gq_last[('m',ci)],
                                                          elementRes[ci])
