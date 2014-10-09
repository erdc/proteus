#!/usr/bin/env python

from proteus import MeshTools,FemTools,Comm,Quadrature,cfemIntegrals,Archiver
import numpy,math
"""
particle tracking tools

TODO
   Clean up TrackingChooser and debug
   Decide where things like skipping points with zero solution should be set (outside of track routine,
    but in the Tracker class or done by the calling routine when setting up flag_track?)

  Decide what the best approach is for copying velocity information. Should these be shallow copies (current)
    or deep copies. Should we add a flag to provide the deep copy option?

"""
import ctracking,ftracking
def TrackingChooser(mesh,nd,problemCoefficients):
    assert False
    tracker = None
    if nd == 1 and problemCoefficients.velocitySpaceFlag == 'c0p1':
        tracker = SteadyState_LinearAdvection_C0P1Velocity_AnalyticalTracking_1d(mesh,nd,
                                                                                 coefficients.velocity_l2g,
                                                                                 coefficients.velocity_dof)
    elif nd == 2 and problemCoefficients.velocitySpaceFlag == 'rt0':
        tracker = SteadyState_LinearAdvection_RT0Velocity_AnalyticalTracking_2d(mesh,nd,
                                                                                coefficients.velocity_l2g,
                                                                                coefficients.velocity_dof)
    return tracker

class Tracking_Base:
    """
    crude interface for particle tracking algorithms
    """
    def __init__(self,mesh,nd,uDict=None,activeComponentList=[0],params={}):
        self.mesh    = mesh    #spatial mesh
        self.nd      = nd      #space dimension
        self.params  = params  #algorithm parameters (e.g., zero tolerances etc)
        self.uDict   = uDict
        self.activeComponentList = activeComponentList #which components of solution are relevant for the tracking
    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        raise NotImplementedError
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        raise NotImplementedError
    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel
        """
        pass

    def setFromOptions(self,options):
        """
        chance to set parameter flags etc
        """
        if 'particleTracking_params' in dir(options):
            for key in options.particleTracking_params.keys():
                self.params[key]=options.particleTracking_params[key]
    def updateTransportInformation(self,transport):
        """
        hack for now to allow tracking to grab information it needs from the transport
        and its coefficients until I figure out a better interface
        """
        pass

class SteadyState_LinearAdvection_C0P1Velocity_AnalyticalTracking_1d(Tracking_Base):
    """
    exact element by element tracking for a steady piecewise linear velocity field in 1D
    TODO
      setFromOptions
    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs={0:None}, #dictionary of component velocity dofs
                 activeComponentList=[0],
                 params={('zeroTol',0):1.0e-5}):
        Tracking_Base.__init__(self,mesh,nd,activeComponentList=activeComponentList,params=params)
        self.component_velocity_l2g=component_velocity_l2g
        self.component_velocity_dofs=component_velocity_dofs
    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1.0
        for ci in self.activeComponentList:
            ctracking.trackPointsC0P1Velocity1d(self.mesh.nElements_global,
                                                self.mesh.nNodes_global,
                                                self.mesh.nNodes_element,
                                                self.mesh.nElementBoundaries_element,
                                                self.mesh.nodeArray,
                                                self.mesh.elementNodesArray,
                                                self.mesh.elementNeighborsArray,
                                                self.component_velocity_l2g[ci],
                                                self.component_velocity_dofs[ci],
                                                direction,
                                                t_depart[ci],
                                                t_track[ci],
                                                nPoints[ci],
                                                self.params[('zeroTol',ci)],
                                                x_depart[ci],
                                                element_track[ci],
                                                x_track[ci],
                                                flag_track[ci])
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = -1.0
        for ci in self.activeComponentList:
            ctracking.trackPointsC0P1Velocity1d(self.mesh.nElements_global,
                                                self.mesh.nNodes_global,
                                                self.mesh.nNodes_element,
                                                self.mesh.nElementBoundaries_element,
                                                self.mesh.nodeArray,
                                                self.mesh.elementNodesArray,
                                                self.mesh.elementNeighborsArray,
                                                self.component_velocity_l2g[ci],
                                                self.component_velocity_dofs[ci],
                                                direction,
                                                t_depart[ci],
                                                t_track[ci],
                                                nPoints[ci],
                                                self.params[('zeroTol',ci)],
                                                x_depart[ci],
                                                element_track[ci],
                                                x_track[ci],
                                                flag_track[ci])

    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel
        ignores time information right now since steady state
        """
        self.component_velocity_dofs[component] = trackingVelocity_dofs
        if trackingVelocity_l2g != None:
            self.component_velocity_l2g[component] = trackingVelocity_l2g


class SteadyState_LinearAdvection_RT0Velocity_AnalyticalTracking_2d(Tracking_Base):
    """
    exact element by element tracking for a steady velocity field with local representation in RT0 in 2D
    TODO
      turn local velocity reprensentation flag into an enum for this and postprocessing
    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs={0:None}, #dictionary of component velocity dofs
                 activeComponentList=[0],
                 params={('zeroTol',0):1.0e-5,
                         ('localVelocityRepresentationFlag',0):2}):
        Tracking_Base.__init__(self,mesh,nd,activeComponentList=activeComponentList,params=params)
        self.component_velocity_l2g=component_velocity_l2g
        self.component_velocity_dofs=component_velocity_dofs

    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1.0
        for ci in self.activeComponentList:
            ctracking.trackPointsRT0Velocity2d(self.params[('localVelocityRepresentationFlag',ci)],
                                               self.mesh.nElements_global,
                                               self.mesh.nNodes_global,
                                               self.mesh.nNodes_element,
                                               self.mesh.nElementBoundaries_element,
                                               self.mesh.nodeArray,
                                               self.mesh.elementNodesArray,
                                               self.mesh.elementNeighborsArray,
                                               self.mesh.elementBoundariesArray,
                                               self.mesh.elementBoundaryBarycentersArray,
                                               self.elementBoundaryOuterNormalsArray,
                                               self.component_velocity_l2g[ci],
                                               self.component_velocity_dofs[ci],
                                               direction,
                                               t_depart[ci],
                                               t_track[ci],
                                               nPoints[ci],
                                               self.params[('zeroTol',ci)],
                                               x_depart[ci],
                                               element_track[ci],
                                               x_track[ci],
                                               flag_track[ci],
                                               0)
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = -1.0
        for ci in self.activeComponentList:
            ctracking.trackPointsRT0Velocity2d(self.params[('localVelocityRepresentationFlag',ci)],
                                               self.mesh.nElements_global,
                                               self.mesh.nNodes_global,
                                               self.mesh.nNodes_element,
                                               self.mesh.nElementBoundaries_element,
                                               self.mesh.nodeArray,
                                               self.mesh.elementNodesArray,
                                               self.mesh.elementNeighborsArray,
                                               self.mesh.elementBoundariesArray,
                                               self.mesh.elementBoundaryBarycentersArray,
                                               self.elementBoundaryOuterNormalsArray,
                                               self.component_velocity_l2g[ci],
                                               self.component_velocity_dofs[ci],
                                               direction,
                                               t_depart[ci],
                                               t_track[ci],
                                               nPoints[ci],
                                               self.params[('zeroTol',ci)],
                                               x_depart[ci],
                                               element_track[ci],
                                               x_track[ci],
                                               flag_track[ci],
                                               0)

    def updateTransportInformation(self,transport):
        """
        hack for now to allow tracking to grab information it needs from the transport
        and its coefficients until I figure out a better interface
        """
        assert 'elementBoundaryOuterNormalsArray' in dir(transport), "RT0 tracking needs elementBoundaryOuterNormalsArray for now"
        self.elementBoundaryOuterNormalsArray = transport.elementBoundaryOuterNormalsArray
    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel
        ignores time information right now since steady state
        """
        self.component_velocity_dofs[component] = trackingVelocity_dofs
        if trackingVelocity_l2g != None:
            self.component_velocity_l2g[component] = trackingVelocity_l2g

class SteadyState_LinearAdvection_C0P1Velocity_PT123(Tracking_Base):
    """
    element by element tracking for a steady continuous, piecewise linear velocity field
      using Pearce's PT123 code
    TODO
      setup input tracking flags to determine if a point is a node or not


    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs={0:None}, #dictionary of component velocity dofs
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('dt_init',0):0.001, #initial time step estimate if using constant value
                         ('dt_init_flag',0):1, #how to pick initial time step (0 --> const, 1 --> CFL=1)
                         ('rk_flag',0):45}):#RK type
        Tracking_Base.__init__(self,mesh,nd,activeComponentList=activeComponentList,params=params)
        #build lookup array for indicating which nodes are on the boundary
        self.nodeOnBoundaryArray = numpy.zeros((self.mesh.nNodes_global,),'i')
        ctracking.setNodeOnBoundaryArray(self.mesh.exteriorElementBoundariesArray,
                                         self.mesh.elementBoundaryNodesArray,
                                         self.nodeOnBoundaryArray)

        #Right now, pt123  uses velocity local to global map that's logically nElements x nNodes x dim
        #standard C0P1 FemSpace uses a l2g that's consistent with a scalar value and assumes extra dofs
        #are just tacked on
        self.component_velocity_dofs=component_velocity_dofs
        #can't have shallow copy in 1d because increment to account for fortran base 1
        #and l2g is itself a shallow copy of elementNodesArray
        #if self.nd == 1: self.component_velocity_l2g=component_velocity_l2g
        #else:
        self.component_velocity_l2g={}
        for ci in component_velocity_l2g.keys():
            if component_velocity_l2g[ci] == None:
                self.component_velocity_l2g[ci]=None
            else:
                nVDOF_element=self.nd*self.mesh.nNodes_element
                self.component_velocity_l2g[ci] = numpy.zeros((self.mesh.nElements_global,self.mesh.nNodes_element,self.nd),'i')
                for I in range(self.nd):
                    self.component_velocity_l2g[ci][:,:,I]=component_velocity_l2g[ci]*self.nd + I


    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            #
            maxeq = 3;
            nVDOF_element=self.nd*self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            #this is supposed to be for steady-state problems so time levels in velocity field shouldn't matter
            #as long as they are not the same
            t_velocity_0 = t_depart[ci].min()
            t_velocity_1 = t_track[ci].max()
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"
            #NOTE current version of pt123 uses velocity local to global map that's logically nElements x nNodes x dim
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase


            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=2,
                            maxeq=3,
                            iverbose=0)

            self.mesh.elementNodesArray -= fbase
            self.mesh.nodeElementsArray -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        no need to set direction = -1 if time values are decreasing
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            maxeq = 3
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #TODO get different velocity degrees of freedom layouts matched
            nVDOF_element=self.nd*self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            #this is supposed to be for steady-state problems so time levels in velocity field shouldn't matter
            t_velocity_0 = t_depart[ci].max()
            t_velocity_1 = t_track[ci].min()
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"

            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase
            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=2,#element based C0-P1
                            maxeq=3,
                            iverbose=0)

            self.mesh.elementNodesArray     -= fbase
            self.mesh.nodeElementsArray     -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase

    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel
        ignores time information right now since steady state
        """
        self.component_velocity_dofs[component] = trackingVelocity_dofs
        if trackingVelocity_l2g != None:
            nVDOF_element=self.nd*self.mesh.nNodes_element
            if trackingVelocity_l2g.shape[-1] == self.mesh.nNodes_element:
                    #standard C0P1 FemSpace uses a l2g that's consistent with a scalar value and assumes extra dofs
                    #are just tacked on
                self.component_velocity_l2g[component] = numpy.zeros((self.mesh.nElements_global,self.mesh.nNodes_element,self.nd),'i')
                for I in range(self.nd):
                    self.component_velocity_l2g[component][:,:,I]=trackingVelocity_l2g*self.nd + I
            else:
                assert ((len(trackingVelocity_l2g.shape) == 2 and trackingVelocity_l2g.shape[-1] == nVDOF_element) or
                        (len(trackingVelocity_l2g.shape) == 3 and trackingVelocity_l2g.shape[-1] == self.nd))
                self.component_velocity_l2g[component] = trackingVelocity_l2g

class SteadyState_LinearAdvection_BDM1Velocity_PT123(SteadyState_LinearAdvection_C0P1Velocity_PT123):
    """
    element by element tracking for a steady BDM1 velocity field
      using Pearce's PT123 code
    TODO
      setup input tracking flags to determine if a point is a node or not


    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs={0:None}, #dictionary of component velocity dofs
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('dt_init',0):0.001, #initial time step estimate if using constant value
                         ('dt_init_flag',0):1, #how to pick initial time step (0 --> const, 1 --> CFL=1)
                         ('rk_flag',0):45}):
        Tracking_Base.__init__(self,mesh,nd,activeComponentList=activeComponentList,params=params)
        #build lookup array for indicating which nodes are on the boundary
        self.nodeOnBoundaryArray = numpy.zeros((self.mesh.nNodes_global,),'i')
        ctracking.setNodeOnBoundaryArray(self.mesh.exteriorElementBoundariesArray,
                                         self.mesh.elementBoundaryNodesArray,
                                         self.nodeOnBoundaryArray)

        #Right now, pt123  uses velocity local to global map that's logically nElements x nNodes x dim
        #The PostProcessingTools BDF representation uses a l2g that uses a basis
        # use standard basis
        #   \vec N_i = \lambda_{i/d} \vec e_{i%d}
        #That is the dofs are locally (say in 2d) [v^x_0,v^y_0,v^x_1,v^y_1,v^x_2,v^y_2]
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.component_velocity_dofs=component_velocity_dofs
        self.component_velocity_l2g=component_velocity_l2g

class SteadyState_LinearAdvection_RT0Velocity_PT123(Tracking_Base):
    """
    element by element tracking for a steady continuous, RT0 velocity field
      on simplices using Pearce's PT123 code
    TODO
      setup input tracking flags to determine if a point is a node or not


    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs={0:None}, #dictionary of component velocity dofs
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('localVelocityRepresentationFlag',0):2,#RT0 velocity representation, 2 -- flux based rep, 1 \vec a + b\vec x
                         ('dt_init',0):0.001, #initial time step estimate if using constant value
                         ('dt_init_flag',0):1, #how to pick initial time step (0 --> const, 1 --> CFL=1)
                         ('rk_flag',0):45}):
        Tracking_Base.__init__(self,mesh,nd,activeComponentList=activeComponentList,params=params)
        #build lookup array for indicating which nodes are on the boundary
        self.nodeOnBoundaryArray = numpy.zeros((self.mesh.nNodes_global,),'i')
        ctracking.setNodeOnBoundaryArray(self.mesh.exteriorElementBoundariesArray,
                                         self.mesh.elementBoundaryNodesArray,
                                         self.nodeOnBoundaryArray)

        #Right now, pt123  uses velocity local to global map that's logically nElements x nNodes x dim
        #standard C0P1 FemSpace uses a l2g that's consistent with a scalar value and assumes extra dofs
        #are just tacked on
        self.component_velocity_dofs=component_velocity_dofs
        self.component_velocity_l2g = component_velocity_l2g

    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            #
            maxeq = 3;
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            #this is supposed to be for steady-state problems so time levels in velocity field shouldn't matter
            #as long as they are not the same
            t_velocity_0 = t_depart[ci].min()
            t_velocity_1 = t_track[ci].max()
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"
            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase


            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=idve,
                            maxeq=3,
                            iverbose=0)

            self.mesh.elementNodesArray -= fbase
            self.mesh.nodeElementsArray -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        no need to set direction = -1 if time values are decreasing
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            maxeq = 3
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #TODO get different velocity degrees of freedom layouts matched
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            #this is supposed to be for steady-state problems so time levels in velocity field shouldn't matter
            t_velocity_0 = t_depart[ci].max()
            t_velocity_1 = t_track[ci].min()

            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase
            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=idve,#element based C0-P1
                            maxeq=3,
                            iverbose=0)

            self.mesh.elementNodesArray     -= fbase
            self.mesh.nodeElementsArray     -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase

    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel
        ignores time information right now since steady state
        """
        self.component_velocity_dofs[component] = trackingVelocity_dofs
        if trackingVelocity_l2g != None:
            self.component_velocity_l2g[component] = trackingVelocity_l2g

class SteadyState_LinearAdvection_RT0Velocity_PT123A(Tracking_Base):
    """
    element by element tracking for a steady continuous, RT0 velocity field
      on simplices using Pearce's PT123 code
    TODO
      setup input tracking flags to determine if a point is a node or not


    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs={0:None}, #dictionary of component velocity dofs
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('localVelocityRepresentationFlag',0):2,#RT0 velocity representation, 2 -- flux based rep, 1 \vec a + b\vec x
                         ('dt_init',0):0.001, #initial time step estimate if using constant value
                         ('dt_init_flag',0):1, #how to pick initial time step (0 --> const, 1 --> CFL=1)
                         ('rk_flag',0):45}):
        Tracking_Base.__init__(self,mesh,nd,activeComponentList=activeComponentList,params=params)
        #build lookup array for indicating which nodes are on the boundary
        self.nodeOnBoundaryArray = numpy.zeros((self.mesh.nNodes_global,),'i')
        ctracking.setNodeOnBoundaryArray(self.mesh.exteriorElementBoundariesArray,
                                         self.mesh.elementBoundaryNodesArray,
                                         self.nodeOnBoundaryArray)

        #Right now, pt123  uses velocity local to global map that's logically nElements x nNodes x dim
        #standard C0P1 FemSpace uses a l2g that's consistent with a scalar value and assumes extra dofs
        #are just tacked on
        self.component_velocity_dofs=component_velocity_dofs
        self.component_velocity_l2g = component_velocity_l2g

    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            #
            maxeq = 3;
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            #this is supposed to be for steady-state problems so time levels in velocity field shouldn't matter
            #as long as they are not the same
            t_velocity_0 = t_depart[ci].min()
            t_velocity_1 = t_track[ci].max()
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"
            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase

            #mwf debug
            #import pdb
            #pdb.set_trace()
            ftracking.pt123a(self.mesh.nNodes_global,
                             self.mesh.nElements_global,
                             self.mesh.nNodes_element,
                             self.nd,
                             nPoints[ci],
                             direction,
                             output_id,
                             self.params[('atol_tracking',ci)],
                             self.params[('rtol_tracking',ci)],
                             self.params[('sf_tracking',ci)],
                             self.params[('dn_safe_tracking',ci)],
                             self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                             self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                             self.mesh.nodeElementOffsets,
                             self.mesh.nodeElementsArray,
                             self.nodeOnBoundaryArray,
                             self.elementBoundaryOuterNormalsArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*self.nd),
                             self.elementBoundaryBarycentersArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*3),
                             self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                             self.component_velocity_dofs[ci].reshape(nVDOF_total),
                             self.component_velocity_dofs[ci].reshape(nVDOF_total),
                             t_velocity_0,
                             t_velocity_1,
                             x_depart[ci].reshape(nPoints[ci]*maxeq),
                             t_depart[ci].reshape(nPoints[ci]),
                             t_track[ci].reshape(nPoints[ci]),
                             element_track[ci].reshape(nPoints[ci]),
                             flag_track[ci].reshape(nPoints[ci]),
                             x_track[ci].reshape(nPoints[ci]*maxeq),
                             idve=idve,
                             maxeq=3,
                             iverbose=0)

            self.mesh.elementNodesArray -= fbase
            self.mesh.nodeElementsArray -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        no need to set direction = -1 if time values are decreasing
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            maxeq = 3
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #TODO get different velocity degrees of freedom layouts matched
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            #this is supposed to be for steady-state problems so time levels in velocity field shouldn't matter
            t_velocity_0 = t_depart[ci].max()
            t_velocity_1 = t_track[ci].min()

            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase
            #mwf debug
            #import pdb
            #pdb.set_trace()
            ftracking.pt123a(self.mesh.nNodes_global,
                             self.mesh.nElements_global,
                             self.mesh.nNodes_element,
                             self.nd,
                             nPoints[ci],
                             direction,
                             output_id,
                             self.params[('atol_tracking',ci)],
                             self.params[('rtol_tracking',ci)],
                             self.params[('sf_tracking',ci)],
                             self.params[('dn_safe_tracking',ci)],
                             self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                             self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                             self.mesh.nodeElementOffsets,
                             self.mesh.nodeElementsArray,
                             self.nodeOnBoundaryArray,
                             self.elementBoundaryOuterNormalsArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*self.nd),
                             self.elementBoundaryBarycentersArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*3),
                             self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                             self.component_velocity_dofs[ci].reshape(nVDOF_total),
                             self.component_velocity_dofs[ci].reshape(nVDOF_total),
                             t_velocity_0,
                             t_velocity_1,
                             x_depart[ci].reshape(nPoints[ci]*maxeq),
                             t_depart[ci].reshape(nPoints[ci]),
                             t_track[ci].reshape(nPoints[ci]),
                             element_track[ci].reshape(nPoints[ci]),
                             flag_track[ci].reshape(nPoints[ci]),
                             x_track[ci].reshape(nPoints[ci]*maxeq),
                             idve=idve,#element based C0-P1
                             maxeq=3,
                             iverbose=0)

            self.mesh.elementNodesArray     -= fbase
            self.mesh.nodeElementsArray     -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase

    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel
        ignores time information right now since steady state
        """
        self.component_velocity_dofs[component] = trackingVelocity_dofs
        if trackingVelocity_l2g != None:
            self.component_velocity_l2g[component] = trackingVelocity_l2g

    def updateTransportInformation(self,transport):
        """
        hack for now to allow tracking to grab information it needs from the transport
        and its coefficients until I figure out a better interface
        """
        assert 'elementBoundaryOuterNormalsArray' in dir(transport), "RT0 tracking needs elementBoundaryOuterNormalsArray for now"
        self.elementBoundaryOuterNormalsArray = transport.elementBoundaryOuterNormalsArray
        self.elementBoundaryBarycentersArray = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,3),'d')
        #todo put this loop in c
        for eN in range(self.mesh.nElements_global):
            for ebN in range(self.mesh.nElementBoundaries_element):
                ebN_global = self.mesh.elementBoundariesArray[eN,ebN]
                self.elementBoundaryBarycentersArray[eN,ebN,:] = self.mesh.elementBoundaryBarycentersArray[ebN_global,:]
class LinearAdvection_C0P1Velocity_PT123(SteadyState_LinearAdvection_C0P1Velocity_PT123):
    """
    element by element tracking for a linear in time, continous piecewise linear velocity field
      using Pearce's PT123 code
    Note that the velocity and times dictionaries are shallow copies
    TODO

    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs_0={0:None}, #dictionary of component velocity dofs at time level 1
                 component_velocity_dofs_1={0:None},
                 component_velocity_times_0={0:None},#time levels for velocities
                 component_velocity_times_1={0:None},
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('dt_init',0):0.001, #initial time step estimate if using constant value
                         ('dt_init_flag',0):1, #how to pick initial time step (0 --> const, 1 --> CFL=1)
                         ('rk_flag',0):45}):
        SteadyState_LinearAdvection_C0P1Velocity_PT123.__init__(self,mesh,nd,
                                                                component_velocity_l2g,
                                                                component_velocity_dofs_0,
                                                                activeComponentList=activeComponentList,
                                                                params=params)
        self.component_velocity_times_0 = component_velocity_times_0
        self.component_velocity_times_1 = component_velocity_times_1
        #
        self.component_velocity_dofs_0  = self.component_velocity_dofs
        self.component_velocity_dofs_1  = component_velocity_dofs_1

    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            #
            maxeq = 3;
            nVDOF_element=self.nd*self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"
            #NOTE current version of pt123 uses velocity local to global map that's logically nElements x nNodes x dim
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase


            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=2,#element based C0P1 velocity field
                            maxeq=3,
                            iverbose=1)

            self.mesh.elementNodesArray -= fbase
            self.mesh.nodeElementsArray -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        no need to set direction = -1 if time values are decreasing
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            maxeq = 3
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #TODO get different velocity degrees of freedom layouts matched
            nVDOF_element=self.nd*self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"

            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase
            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=2,#element based C0P1 velocity field
                            maxeq=3,
                            iverbose=0)

            self.mesh.elementNodesArray     -= fbase
            self.mesh.nodeElementsArray     -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel

        """
        if trackingVelocity_l2g != None:
            nVDOF_element=self.nd*self.mesh.nNodes_element
            if trackingVelocity_l2g.shape[-1] == self.mesh.nNodes_element:
                #standard C0P1 FemSpace uses a l2g that's consistent with a scalar value and assumes extra dofs
                #are just tacked on
                self.component_velocity_l2g[component] = numpy.zeros((self.mesh.nElements_global,self.mesh.nNodes_element,self.nd),'i')
                for I in range(self.nd):
                    self.component_velocity_l2g[component][:,:,I]=trackingVelocity_l2g*self.nd + I
            else:
                assert ((len(trackingVelocity_l2g.shape) == 2 and trackingVelocity_l2g.shape[-1] == nVDOF_element) or
                        (len(trackingVelocity_l2g.shape) == 3 and trackingVelocity_l2g.shape[-1] == self.nd))
                self.component_velocity_l2g[component] = trackingVelocity_l2g

        if timeLevel == 1:
            self.component_velocity_dofs_1[component] = trackingVelocity_dofs
            self.component_velocity_times_1[component]= time
        else:
            self.component_velocity_dofs_0[component] = trackingVelocity_dofs
            self.component_velocity_times_0[component]= time


class LinearAdvection_BDM1Velocity_PT123(LinearAdvection_C0P1Velocity_PT123):
    """
    element by element tracking for a linear in time,  velocity field with BDM1 representation in space
      using Pearce's PT123 code
    TODO
      setup input tracking flags to determine if a point is a node or not


    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs_0={0:None}, #dictionary of component velocity dofs at time level 1
                 component_velocity_dofs_1={0:None},
                 component_velocity_times_0={0:None},#time levels for velocities
                 component_velocity_times_1={0:None},
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('dt_init',0):0.001, #initial time step estimate if using constant value
                         ('dt_init_flag',0):1, #how to pick initial time step (0 --> const, 1 --> CFL=1)
                         ('rk_flag',0):45}):#RK type
        Tracking_Base.__init__(self,mesh,nd,activeComponentList=activeComponentList,params=params)
        #build lookup array for indicating which nodes are on the boundary
        self.nodeOnBoundaryArray = numpy.zeros((self.mesh.nNodes_global,),'i')
        ctracking.setNodeOnBoundaryArray(self.mesh.exteriorElementBoundariesArray,
                                         self.mesh.elementBoundaryNodesArray,
                                         self.nodeOnBoundaryArray)

        self.component_velocity_l2g=component_velocity_l2g
        #just to be consistent with C0P1 for now
        self.component_velocity_dofs=component_velocity_dofs_0
        #
        self.component_velocity_times_0 = component_velocity_times_0
        self.component_velocity_times_1 = component_velocity_times_1
        #
        self.component_velocity_dofs_0  = self.component_velocity_dofs
        self.component_velocity_dofs_1  = component_velocity_dofs_1

class LinearAdvection_RT0Velocity_PT123(SteadyState_LinearAdvection_RT0Velocity_PT123):
    """
    element by element tracking for a linear in time, RT0 velocity field (on simplices)
      using Pearce's PT123 code
    Note that the velocity and times dictionaries are shallow copies
    TODO

    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs_0={0:None}, #dictionary of component velocity dofs at time level 1
                 component_velocity_dofs_1={0:None},
                 component_velocity_times_0={0:None},#time levels for velocities
                 component_velocity_times_1={0:None},
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('localVelocityRepresentationFlag',0):2,#RT0 velocity representation, 2 -- flux based rep, 1 \vec a + b\vec x
                         ('dt_init',0):0.001, #initial time step estimate if using constant value
                         ('dt_init_flag',0):1, #how to pick initial time step (0 --> const, 1 --> CFL=1)
                         ('rk_flag',0):45}):#RK type
        SteadyState_LinearAdvection_RT0Velocity_PT123.__init__(self,mesh,nd,
                                                               component_velocity_l2g,
                                                               component_velocity_dofs_0,
                                                               activeComponentList=activeComponentList,
                                                               params=params)
        self.component_velocity_times_0 = component_velocity_times_0
        self.component_velocity_times_1 = component_velocity_times_1
        #
        self.component_velocity_dofs_0  = self.component_velocity_dofs
        self.component_velocity_dofs_1  = component_velocity_dofs_1

    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            #
            maxeq = 3;
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase


            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=idve,
                            maxeq=3,
                            iverbose=0)

            self.mesh.elementNodesArray -= fbase
            self.mesh.nodeElementsArray -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        no need to set direction = -1 if time values are decreasing
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            maxeq = 3
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #TODO get different velocity degrees of freedom layouts matched
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]

            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase
            ftracking.pt123(self.mesh.nNodes_global,
                            self.mesh.nElements_global,
                            self.mesh.nNodes_element,
                            self.nd,
                            nPoints[ci],
                            direction,
                            output_id,
                            self.params[('atol_tracking',ci)],
                            self.params[('rtol_tracking',ci)],
                            self.params[('sf_tracking',ci)],
                            self.params[('dn_safe_tracking',ci)],
                            self.params[('dt_init',ci)],
                            self.params[('dt_init_flag',ci)],self.params[('rk_flag',ci)],
                            self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                            self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                            self.mesh.elementDiametersArray,
                            self.mesh.nodeElementOffsets,
                            self.mesh.nodeElementsArray,
                            self.nodeOnBoundaryArray,
                            self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                            self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                            self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                            t_velocity_0,
                            t_velocity_1,
                            x_depart[ci].reshape(nPoints[ci]*maxeq),
                            t_depart[ci].reshape(nPoints[ci]),
                            t_track[ci].reshape(nPoints[ci]),
                            element_track[ci].reshape(nPoints[ci]),
                            flag_track[ci].reshape(nPoints[ci]),
                            x_track[ci].reshape(nPoints[ci]*maxeq),
                            idve=idve,
                            maxeq=3,
                            iverbose=0)

            self.mesh.elementNodesArray     -= fbase
            self.mesh.nodeElementsArray     -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel

        """
        if trackingVelocity_l2g != None:
            self.component_velocity_l2g[component]    = trackingVelocity_l2g
        if timeLevel == 1:
            self.component_velocity_dofs_1[component] = trackingVelocity_dofs
            self.component_velocity_times_1[component]= time
        else:
            self.component_velocity_dofs_0[component] = trackingVelocity_dofs
            self.component_velocity_times_0[component]= time


class LinearAdvection_RT0Velocity_PT123A(SteadyState_LinearAdvection_RT0Velocity_PT123A):
    """
    analytical element by element tracking for a RT0 velocity field (on simplices)
      V = V_0 + v_x (X-X_0) + V_t (t-t_0) + v_xt (X-X_0)(t-t_0)

    using Pearce's PT123A framework and Russell, Healy 00 approach (sort of)
    idvt = 1 --> v_xt = 0

    Note that the velocity and times dictionaries are shallow copies
    TODO

    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs_0={0:None}, #dictionary of component velocity dofs at time level 1
                 component_velocity_dofs_1={0:None},
                 component_velocity_times_0={0:None},#time levels for velocities
                 component_velocity_times_1={0:None},
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7,#tolerance for traking in-element tests
                         ('localVelocityRepresentationFlag',0):2,#RT0 velocity representation, 2 -- flux based rep, 1 \vec a + b\vec x
                         ('temporalVariationFlag',0):1}):
        SteadyState_LinearAdvection_RT0Velocity_PT123A.__init__(self,mesh,nd,
                                                               component_velocity_l2g,
                                                               component_velocity_dofs_0,
                                                               activeComponentList=activeComponentList,
                                                               params=params)
        self.component_velocity_times_0 = component_velocity_times_0
        self.component_velocity_times_1 = component_velocity_times_1
        #
        self.component_velocity_dofs_0  = self.component_velocity_dofs
        self.component_velocity_dofs_1  = component_velocity_dofs_1

    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            #
            maxeq = 3;
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            #how is local velocity represented?
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            #what type of temporal variation is assumed, 0 -- steady state, 1 -- linear in x-t, 2 -- bilinear in x-t
            assert self.params[('temporalVariationFlag',ci)] in [0,1,2], "self.params[('temporalVariationFlag',ci)]= %s must be in [0,1,2] " % (self.params[('temporalVariationFlag',ci)])
            idvt = self.params[('temporalVariationFlag',ci)]

            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase


            ftracking.pt123a(self.mesh.nNodes_global,
                             self.mesh.nElements_global,
                             self.mesh.nNodes_element,
                             self.nd,
                             nPoints[ci],
                             direction,
                             output_id,
                             self.params[('atol_tracking',ci)],
                             self.params[('rtol_tracking',ci)],
                             self.params[('sf_tracking',ci)],
                             self.params[('dn_safe_tracking',ci)],
                             self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                             self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                             self.mesh.nodeElementOffsets,
                             self.mesh.nodeElementsArray,
                             self.nodeOnBoundaryArray,
                             self.elementBoundaryOuterNormalsArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*self.nd),
                             self.elementBoundaryBarycentersArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*3),
                             self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                             self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                             self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                             t_velocity_0,
                             t_velocity_1,
                             x_depart[ci].reshape(nPoints[ci]*maxeq),
                             t_depart[ci].reshape(nPoints[ci]),
                             t_track[ci].reshape(nPoints[ci]),
                             element_track[ci].reshape(nPoints[ci]),
                             flag_track[ci].reshape(nPoints[ci]),
                             x_track[ci].reshape(nPoints[ci]*maxeq),
                             idve=idve,
                             idvt=idvt,
                             maxeq=3,
                             iverbose=0)

            self.mesh.elementNodesArray -= fbase
            self.mesh.nodeElementsArray -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        no need to set direction = -1 if time values are decreasing
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            maxeq = 3
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #TODO get different velocity degrees of freedom layouts matched
            nVDOF_element=self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            assert self.params[('localVelocityRepresentationFlag',ci)] in [1,2], "self.params[('localVelocityRepresentationFlag',ci)]= %s not supported " % (self.params[('localVelocityRepresentationFlag',ci)])
            idve = 3
            if self.params[('localVelocityRepresentationFlag',ci)] == 1:
                idve = 4
            #what type of temporal variation is assumed, 0 -- steady state, 1 -- linear in x-t, 2 -- bilinear in x-t
            assert self.params[('temporalVariationFlag',ci)] in [0,1,2], "self.params[('temporalVariationFlag',ci)]= %s must be in [0,1,2] " % (self.params[('temporalVariationFlag',ci)])
            idvt = self.params[('temporalVariationFlag',ci)]
            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]

            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase
            ftracking.pt123a(self.mesh.nNodes_global,
                             self.mesh.nElements_global,
                             self.mesh.nNodes_element,
                             self.nd,
                             nPoints[ci],
                             direction,
                             output_id,
                             self.params[('atol_tracking',ci)],
                             self.params[('rtol_tracking',ci)],
                             self.params[('sf_tracking',ci)],
                             self.params[('dn_safe_tracking',ci)],
                             self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                             self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                             self.mesh.nodeElementOffsets,
                             self.mesh.nodeElementsArray,
                             self.nodeOnBoundaryArray,
                             self.elementBoundaryOuterNormalsArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*self.nd),
                             self.elementBoundaryBarycentersArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*3),
                             self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                             self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                             self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                             t_velocity_0,
                             t_velocity_1,
                             x_depart[ci].reshape(nPoints[ci]*maxeq),
                             t_depart[ci].reshape(nPoints[ci]),
                             t_track[ci].reshape(nPoints[ci]),
                             element_track[ci].reshape(nPoints[ci]),
                             flag_track[ci].reshape(nPoints[ci]),
                             x_track[ci].reshape(nPoints[ci]*maxeq),
                             idve=idve,
                             idvt=idvt,
                             maxeq=3,
                             iverbose=0)

            self.mesh.elementNodesArray     -= fbase
            self.mesh.nodeElementsArray     -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel

        """
        if trackingVelocity_l2g != None:
            self.component_velocity_l2g[component]    = trackingVelocity_l2g
        if timeLevel == 1:
            self.component_velocity_dofs_1[component] = trackingVelocity_dofs
            self.component_velocity_times_1[component]= time
        else:
            self.component_velocity_dofs_0[component] = trackingVelocity_dofs
            self.component_velocity_times_0[component]= time


class LinearAdvection_C0P1Velocity_PT123A(SteadyState_LinearAdvection_C0P1Velocity_PT123):
    """
    analytical element by element tracking for a velocity field (on simplices)
     V = V_0 + v_x (X-X_0) + V_t (t-t_0)
      using Pearce's PT123A framework and Russell, Healy 00 approach (sort of)
    Apply this to a C0P1 velocity field and allow some error
    Note that the velocity and times dictionaries are shallow copies
    TODO

    """
    def __init__(self,mesh,nd,
                 component_velocity_l2g={0:None}, #dictionary of component velocity l2g mappings
                 component_velocity_dofs_0={0:None}, #dictionary of component velocity dofs at time level 1
                 component_velocity_dofs_1={0:None},
                 component_velocity_times_0={0:None},#time levels for velocities
                 component_velocity_times_1={0:None},
                 activeComponentList=[0],
                 params={('atol_tracking',0):1.0e-7,#RK time integration tolerances
                         ('rtol_tracking',0):0.0,
                         ('sf_tracking',0):0.9,     #safety factor for RK integration
                         ('dn_safe_tracking',0):1.0e-7}):#tolerance for traking in-element tests
        SteadyState_LinearAdvection_C0P1Velocity_PT123.__init__(self,mesh,nd,
                                                                component_velocity_l2g,
                                                                component_velocity_dofs_0,
                                                                activeComponentList=activeComponentList,
                                                                params=params)
        self.component_velocity_times_0 = component_velocity_times_0
        self.component_velocity_times_1 = component_velocity_times_1
        #
        self.component_velocity_dofs_0  = self.component_velocity_dofs
        self.component_velocity_dofs_1  = component_velocity_dofs_1

    def forwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            #
            maxeq = 3;
            nVDOF_element=self.nd*self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"
            #NOTE current version of pt123 uses velocity local to global map that's logically nElements x nNodes x dim
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase

            ftracking.pt123a(self.mesh.nNodes_global,
                             self.mesh.nElements_global,
                             self.mesh.nNodes_element,
                             self.nd,
                             nPoints[ci],
                             direction,
                             output_id,
                             self.params[('atol_tracking',ci)],
                             self.params[('rtol_tracking',ci)],
                             self.params[('sf_tracking',ci)],
                             self.params[('dn_safe_tracking',ci)],
                             self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                             self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                             self.mesh.nodeElementOffsets,
                             self.mesh.nodeElementsArray,
                             self.nodeOnBoundaryArray,
                             self.elementBoundaryOuterNormalsArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*self.nd),
                             self.elementBoundaryBarycentersArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*3),
                             self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                             self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                             self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                             t_velocity_0,
                             t_velocity_1,
                             x_depart[ci].reshape(nPoints[ci]*maxeq),
                             t_depart[ci].reshape(nPoints[ci]),
                             t_track[ci].reshape(nPoints[ci]),
                             element_track[ci].reshape(nPoints[ci]),
                             flag_track[ci].reshape(nPoints[ci]),
                             x_track[ci].reshape(nPoints[ci]*maxeq),
                             idve=2,
                             idvt=1,
                             maxeq=3,
                             iverbose=0)

            self.mesh.elementNodesArray -= fbase
            self.mesh.nodeElementsArray -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase
    def backwardTrack(self,
                     t_depart,           #point departure times
                     t_track,            #target end time
                     nPoints,            #total number of points
                     x_depart,           #departure points
                     element_track,      #in/out element locations
                     x_track,            #arrival points
                     flag_track):
        """
        track quadrature points in x_depart forward from t_depart --> t_track (input is target)
        loads
           x_track   : location of point at end of tracking
           t_track   : time tracking ended
           flag_track: -3  did not track (e.g., v = 0 or u = 0 and algorithm has flag to skip
                           zero solution values)
                       -2  point exited the domain somewhere in (t_depart,t_target)
                       -1  point in interior at t_target
                      >=0  reserved to indicate mesh vertex id (e.g., for Pearce's tracking)
          element_track : element containing point at end of tracking
        no need to set direction = -1 if time values are decreasing
        """
        direction = 1; output_id=6;
        for ci in self.activeComponentList:
            maxeq = 3
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #TODO get different velocity degrees of freedom layouts matched
            nVDOF_element=self.nd*self.mesh.nNodes_element
            nVDOF_total  =len(self.component_velocity_dofs[ci].flat)
            t_velocity_0 = self.component_velocity_times_0[ci]
            t_velocity_1 = self.component_velocity_times_1[ci]
            #assert abs(t_velocity_1-t_velocity_0) > 0.0, "pt123 requires velocity tracking times different"

            #rather than modify fortran code that assume's base 1 node, element id's etc just update the
            #c data structures here and then convert back to base zero after call
            fbase = 1
            self.mesh.elementNodesArray += fbase
            self.mesh.nodeElementsArray += fbase
            self.component_velocity_l2g[ci] += fbase
            element_track[ci] += fbase
            flag_track[ci]    += fbase
            ftracking.pt123a(self.mesh.nNodes_global,
                             self.mesh.nElements_global,
                             self.mesh.nNodes_element,
                             self.nd,
                             nPoints[ci],
                             direction,
                             output_id,
                             self.params[('atol_tracking',ci)],
                             self.params[('rtol_tracking',ci)],
                             self.params[('sf_tracking',ci)],
                             self.params[('dn_safe_tracking',ci)],
                             self.mesh.nodeArray.reshape(self.mesh.nNodes_global*maxeq),
                             self.mesh.elementNodesArray.reshape(self.mesh.nElements_global*self.mesh.nNodes_element),
                             self.mesh.nodeElementOffsets,
                             self.mesh.nodeElementsArray,
                             self.nodeOnBoundaryArray,
                             self.elementBoundaryOuterNormalsArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*self.nd),
                             self.elementBoundaryBarycentersArray.reshape(self.mesh.nElements_global*self.mesh.nElementBoundaries_element*3),
                             self.component_velocity_l2g[ci].reshape(self.mesh.nElements_global*nVDOF_element),
                             self.component_velocity_dofs_0[ci].reshape(nVDOF_total),
                             self.component_velocity_dofs_1[ci].reshape(nVDOF_total),
                             t_velocity_0,
                             t_velocity_1,
                             x_depart[ci].reshape(nPoints[ci]*maxeq),
                             t_depart[ci].reshape(nPoints[ci]),
                             t_track[ci].reshape(nPoints[ci]),
                             element_track[ci].reshape(nPoints[ci]),
                             flag_track[ci].reshape(nPoints[ci]),
                             x_track[ci].reshape(nPoints[ci]*maxeq),
                             idve=2,
                             idvt=1,
                             maxeq=3,
                             iverbose=0)

            self.mesh.elementNodesArray     -= fbase
            self.mesh.nodeElementsArray     -= fbase
            self.component_velocity_l2g[ci] -= fbase
            element_track[ci] -= fbase
            flag_track[ci]    -= fbase

    def setTrackingVelocity(self,trackingVelocity_dofs,component,time,timeLevel=0,trackingVelocity_l2g=None):
        """
        load tracking velocity degrees of freedom for component at time and call it timeLevel

        """
        if trackingVelocity_l2g != None:
            nVDOF_element=self.nd*self.mesh.nNodes_element
            if trackingVelocity_l2g.shape[-1] == self.mesh.nNodes_element:
                #standard C0P1 FemSpace uses a l2g that's consistent with a scalar value and assumes extra dofs
                #are just tacked on
                self.component_velocity_l2g[component] = numpy.zeros((self.mesh.nElements_global,self.mesh.nNodes_element,self.nd),'i')
                for I in range(self.nd):
                    self.component_velocity_l2g[component][:,:,I]=trackingVelocity_l2g*self.nd + I
            else:
                assert ((len(trackingVelocity_l2g.shape) == 2 and trackingVelocity_l2g.shape[-1] == nVDOF_element) or
                        (len(trackingVelocity_l2g.shape) == 3 and trackingVelocity_l2g.shape[-1] == self.nd))
                self.component_velocity_l2g[component] = trackingVelocity_l2g

        if timeLevel == 1:
            self.component_velocity_dofs_1[component] = trackingVelocity_dofs
            self.component_velocity_times_1[component]= time
        else:
            self.component_velocity_dofs_0[component] = trackingVelocity_dofs
            self.component_velocity_times_0[component]= time


    def updateTransportInformation(self,transport):
        """
        hack for now to allow tracking to grab information it needs from the transport
        and its coefficients until I figure out a better interface
        """
        assert 'elementBoundaryOuterNormalsArray' in dir(transport), "Analytical tracking needs elementBoundaryOuterNormalsArray for now"
        self.elementBoundaryOuterNormalsArray = transport.elementBoundaryOuterNormalsArray
        self.elementBoundaryBarycentersArray = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,3),'d')
        #todo put this loop in c
        for eN in range(self.mesh.nElements_global):
            for ebN in range(self.mesh.nElementBoundaries_element):
                ebN_global = self.mesh.elementBoundariesArray[eN,ebN]
                self.elementBoundaryBarycentersArray[eN,ebN,:] = self.mesh.elementBoundaryBarycentersArray[ebN_global,:]



######################################################################
#routines for testing tracking
def setupMesh_1d(opts,p):
    assert p.nd == 1, "1d only for now"
    #spatial mesh,
    mlMesh = MeshTools.MultilevelEdgeMesh(opts.nnx,1,1,p.L[0],refinementLevels=1)
    #track on finest level only
    mesh=mlMesh.meshList[-1]
    return mesh
def setupMesh_2d(opts,p):
    assert p.nd == 2, "2d only for now"
    #spatial mesh,
    mlMesh = MeshTools.MultilevelTriangularMesh(opts.nnx,opts.nny,1,p.L[0],p.L[1],refinementLevels=1)
    #track on finest level only
    mesh=mlMesh.meshList[-1]
    return mesh
def setupMesh_3d(opts,p):
    assert p.nd == 3, "2d only for now"
    #spatial mesh,
    mlMesh = MeshTools.MultilevelTetrahedralMesh(opts.nnx,opts.nny,opts.nnz,p.L[0],p.L[1],p.L[2],refinementLevels=1)
    #track on finest level only
    mesh=mlMesh.meshList[-1]
    return mesh
def setupSolution_C0P1(opts,p,mesh,t):
    """
    returns FEM space, FEM function, and interpolation values
    """
    #solution representation, not used for linear problems
    trialSpace= FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh,p.nd)
    u = FemTools.FiniteElementFunction(trialSpace,name="u")
    #evaluate solution dofs from interpolation conditions using problem's initial condition
    #need to put t0 in p?
    t=opts.t_start
    u_interpolation_values = numpy.zeros((mesh.nElements_global,trialSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints),'d')
    for eN in range(mesh.nElements_global):
        for k in range(trialSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
            u_interpolation_values[eN,k] = p.initialConditions[0].uOfXT(trialSpace.interpolationPoints[eN,k],t)
    #evaluate solution dofs from interpolation conditions
    u.projectFromInterpolationConditions(u_interpolation_values)
    return trialSpace,u,u_interpolation_values
def setupVelocity_C0P1(opts,p,mesh,t):
    """
    returns FEM space and FEM function, and interpolation values
    """
    velocitySpace= FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh,p.nd)
    velocity = FemTools.FiniteElementFunction(velocitySpace,dim_dof=p.nd,isVector=True,name="velocity")
    velocity_interpolation_values = numpy.zeros((mesh.nElements_global,velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints,p.nd),'d')
    for eN in range(mesh.nElements_global):
        for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
            velocity_interpolation_values[eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],t)
    #evaluate velocity degrees of freedom from its interpolation conditions
    velocity.projectFromInterpolationConditions(velocity_interpolation_values)

    return velocitySpace,velocity,velocity_interpolation_values
def setupVelocity_RT0(opts,p,mesh,geometricSpace,t,nd=2):
    """
    initially returns just a l2g mapping and dof map for an RT0 velocity representation on
    simplexes because we don't have a fem function or space for this yet

    On element e:
      \vec v_h = \sum^d_{i=0}V^i\vec N_{e,i}
    for
      \vec N_{e,i} = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}), i=0,...,d

   The degrees of freedom are flux integrals through faces
      V^i = \int_{\gamma_i}\vec v\dot n_{i}\ds

    velocity dofs is logically nElements_global x nElementBoundaries_element

    geometricSpace defines geometry for mesh (affine simplicial for now)
    """
    velocity_dofs = numpy.zeros((mesh.nElements_global*mesh.nElementBoundaries_element),'d')
    velocity_l2g  = numpy.arange((mesh.nElements_global*mesh.nElementBoundaries_element),dtype='i').reshape((mesh.nElements_global,mesh.nElementBoundaries_element))
    velocity_interpolation_values = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element),'d')

    #setup information for boundary integrals
    if nd == 1:
        boundaryQuadrature = Quadrature.GaussPoint()
    elif nd == 3:
        boundaryQuadrature = Quadrature.GaussTriangle(4)
    else:
        boundaryQuadrature = Quadrature.GaussEdge(4)
    elementBoundaryQuadraturePoints = numpy.array([pt for pt in boundaryQuadrature.points],dtype='d')
    elementBoundaryQuadratureWeights= numpy.array([wt for wt in boundaryQuadrature.weights],dtype='d')

    nElementBoundaryQuadraturePoints_elementBoundary = elementBoundaryQuadraturePoints.shape[0]
    ebq = {}
    ebq['x'] = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
    ebq['inverse(J)'] = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,p.nd,p.nd),'d')
    ebq['g'] = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,max(1,p.nd-1),max(1,p.nd-1)),'d')
    ebq['sqrt(det(g))'] = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary),'d')
    ebq['n'] = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,p.nd),'d')
    ebq['dS'] = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary),'d')

    ebq['velocity'] = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,p.nd),'d')

    geometricSpace.elementMaps.getValuesTrace(elementBoundaryQuadraturePoints,
                                              ebq['x'])
    geometricSpace.elementMaps.getJacobianValuesTrace(elementBoundaryQuadraturePoints,
                                                      ebq['inverse(J)'],
                                                      ebq['g'],
                                                      ebq['sqrt(det(g))'],
                                                      ebq['n'])
    cfemIntegrals.calculateElementBoundaryIntegrationWeights(ebq['sqrt(det(g))'],
                                                             elementBoundaryQuadratureWeights,
                                                             ebq['dS'])
    #manually evaluate velocity degrees of freedom from its interpolation conditions
    for eN in range(mesh.nElements_global):
        for ebN in range(mesh.nElementBoundaries_element):
            integral = 0.0
            for kb in range(nElementBoundaryQuadraturePoints_elementBoundary):
                v = p.analyticalSolutionParticleVelocity[0].uOfXT(ebq['x'][eN,ebN,kb],t)
                for I in range(p.nd):
                    integral += v[I]*ebq['n'][eN,ebN,kb,I]*ebq['dS'][eN,ebN,kb]
                #mwf debug
                #print "setup RT0 eN= %s ebN= %s kb=%s v=[%s,%s] n=[%s,%s] integral=%s " % (eN,ebN,kb,v[0],v[1],ebq['n'][eN,ebN,kb,0],ebq['n'][eN,ebN,kb,1],integral)
            velocity_interpolation_values[eN,ebN] = integral
            velocity_dofs[velocity_l2g[eN,ebN]] = velocity_interpolation_values[eN,ebN]

    return velocity_dofs,velocity_l2g,velocity_interpolation_values,ebq,elementBoundaryQuadraturePoints,elementBoundaryQuadratureWeights
def evaluateVelocity_RT0(opts,p,mesh,t,ebq,velocity_dofs,velocity_l2g,velocity_interpolation_values):
    """
    evaluate RT0 velocity

    On element e:
      \vec v_h = \sum^d_{i=0}V^i\vec N_{e,i}
    for
      \vec N_{e,i} = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}), i=0,...,d

   The degrees of freedom are flux integrals through faces
      V^i = \int_{\gamma_i}\vec v\dot n_{i}\ds

    velocity dofs is logically nElements_global x nElementBoundaries_element

    geometricSpace defines geometry for mesh (affine simplicial for now)
    """

    #manually evaluate velocity degrees of freedom from its interpolation conditions
    #todo move this from python
    from subsurfaceTransportFunctions import calculateNormalFlux
    try:
        p.analyticalSolutionParticleVelocity[0].uOfXTv(ebq['x'],t,ebq['velocity'])
        calculateNormalFlux(ebq['velocity'],ebq['n'],ebq['dS'],velocity_interpolation_values)
        velocity_dofs[velocity_l2g] = velocity_interpolation_values
    except:
        print "WARNING Tracking.evaluateVelocity_RT0 try uOfXTv failed"
        for eN in range(mesh.nElements_global):
            for ebN in range(mesh.nElementBoundaries_element):
                integral = 0.0
                for kb in range(ebq['x'].shape[-2]):
                    v = p.analyticalSolutionParticleVelocity[0].uOfXT(ebq['x'][eN,ebN,kb],t)
                    for I in range(p.nd):
                        integral += v[I]*ebq['n'][eN,ebN,kb,I]*ebq['dS'][eN,ebN,kb]
                    #mwf debug
                    #print "setup RT0 eN= %s ebN= %s kb=%s v=[%s,%s] n=[%s,%s] integral=%s " % (eN,ebN,kb,v[0],v[1],ebq['n'][eN,ebN,kb,0],ebq['n'][eN,ebN,kb,1],integral)
                velocity_interpolation_values[eN,ebN] = integral
                velocity_dofs[velocity_l2g[eN,ebN]] = velocity_interpolation_values[eN,ebN]

def setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False):
    """
    return points to track in reference space
    """
    #assumes have the same number of integration points per element to start tracking
    # generate points on reference element, map to domain
    # these could be arbitrary but test first for standard quadrature points
    if p.nd == 1:
        if useGaussianQuadrature:
            quadrature = Quadrature.GaussianEdge(order=opts.nnq)
        else:
            quadrature = Quadrature.CompositeTrapezoidalEdge(order=opts.nnq)
        #(q['x_hat'],q['w_hat'],q['x_hat_indeces']) = Quadrature.buildUnion({0:quadrature})
        x_hat = numpy.array([x for x in quadrature.points],dtype='d')
        return x_hat
    elif p.nd == 2:
        #need CompositeTrapezoidalTriangle
        if useGaussianQuadrature:
            quadrature = Quadrature.GaussTriangle(order=opts.nnq)
        else:
            quadrature = Quadrature.CompositeTrapezoidalTriangle(order=opts.nnq)

        x_hat = numpy.array([x for x in quadrature.points],dtype='d')
        return x_hat
    elif p.nd == 3:
        #need CompositeTrapezoidalTriangle
        quadrature = Quadrature.GaussTetrahedron(order=opts.nnq)

        x_hat = numpy.array([x for x in quadrature.points],dtype='d')
        return x_hat
    else:
        return None
def setupTrackingDataArrays(t_start,t_target,mesh,q,nq_per_element):
    #input/arrays
    #in/out -- element where point is located
    q['element_track'] = numpy.zeros((mesh.nElements_global,nq_per_element),'i')
    for eN in range(mesh.nElements_global):
        q['element_track'][eN,:] = eN
    #in/out -- departure, arrival times for points
    q['t_depart'] = numpy.zeros((mesh.nElements_global,nq_per_element),'d')
    q['t_track'] = numpy.zeros((mesh.nElements_global,nq_per_element),'d')
    q['t_depart'].fill(t_start)
    q['t_track'].fill(t_target)
    #out -- exit points
    q['x_track'] = numpy.zeros((mesh.nElements_global,nq_per_element,3),'d')
    #in/out --- status of point
    #in
    #0 : track point, later can make >= 0 track, 1 --> vertex
    #out
    #0 : point in interior
    #1 : point exited domain
    #-1: did not track
    q['flag_track']= numpy.zeros((mesh.nElements_global,nq_per_element),'i')
    #on input -1 means try, it's an interior point
    q['flag_track'].fill(-1)
def setupArbitraryTrackingDataArrays(t_start,t_target,mesh,elementBoundaryNormals,q):
    """
    setup arrays assuming user has given specific points to track
    """
    #input/arrays
    #in/out -- element where point is located
    q['element_track'] = numpy.zeros(q['x_track'].shape[0],'i')
    nd = elementBoundaryNormals.shape[-1]
    for i in range(q['x_track'].shape[0]):#find where the points are with brute force
        eN = 0; foundElement = False
        while eN < mesh.nElements_global and not foundElement:
            outside = False
            for ebN in range(mesh.nElementBoundaries_element):
                ebN_global = mesh.elementBoundariesArray[eN,ebN]
                dxf = numpy.dot((q['x_track'][i,:nd]-mesh.elementBoundaryBarycentersArray[ebN_global,:nd]),elementBoundaryNormals[eN,ebN])
                outside = outside or dxf > 1.0e-8
            if not outside: foundElement = True
            eN += 1
        assert foundElement
        q['element_track'][i]=eN-1
    #mwf debug
    #import pdb
    #pdb.set_trace()
    #in/out -- departure, arrival times for points
    q['t_depart'] = numpy.zeros(q['x_track'].shape[0],'d')
    q['t_track'] = numpy.zeros(q['x_track'].shape[0],'d')
    q['t_depart'].fill(t_start)
    q['t_track'].fill(t_target)
    #in/out --- status of point
    #in
    #0 : track point, later can make >= 0 track, 1 --> vertex
    #out
    #0 : point in interior
    #1 : point exited domain
    #-1: did not track
    q['flag_track']= numpy.zeros(q['x_track'].shape[0],'i')
    #on input -1 means try, it's an interior point
    q['flag_track'].fill(-1)


def writePT123inputMesh(opts,p,mesh,filebase="test_pt123"):
    """
    create input mesh file for running actual PT123 code
    """
    base = 1 #base 1 numbering
    #write out mesh
    mesh_prefix = '1dm'; ele_label   = 'GE2'; node_label = 'GN'
    if p.nd == 2:
        mesh_prefix = '2dm'; ele_label   = 'GE3';
    elif p.nd == 3:
        mesh_prefix = '3dm'; ele_label   = 'GE4';

    fmesh = open(filebase+'.'+mesh_prefix,'w')
    fmesh.write('MESH \n')
    for eN in range(mesh.nElements_global):
        fmesh.write("%s   %d " % (ele_label,eN+base))
        for nN in range(mesh.nNodes_element):
            fmesh.write("   %d " % (mesh.elementNodesArray[eN,nN]+base))
        fmesh.write("\n")
    #
    for nN in range(mesh.nNodes_global):
        fmesh.write("%s   %d " % (node_label,nN+base))
        for I in range(p.nd):
            fmesh.write("   %g " % (mesh.nodeArray[nN,I]))
        fmesh.write("\n")
    #
    fmesh.write("ENDR\n")
    fmesh.close()
#
def writePT123nodalVelocity(opts,p,mesh,tnList,q,velocitySpace,velocityFEFunction,
                            filebase="test_pt123"):
    #write out velocity file (have to loop through time steps)
    vel_prefix = 'vn1';
    if p.nd == 2:
        vel_prefix = 'vn2';
    elif p.nd == 3:
        vel_prefix = 'vn3';

    fvel = open(filebase+'.'+vel_prefix,'w')
    fvel.write("     %d     %d    %d \n" % (mesh.nNodes_global,p.nd,len(tnList)))

    for i,t in enumerate(tnList):
        fvel.write("TS    %f \n" % t)

        try:
            p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
        except TypeError:
            for eN in range(mesh.nElements_global):
                for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                    q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],t)
        #evaluate velocity degrees of freedom from its interpolation conditions
        velocityFEFunction.projectFromInterpolationConditions(q['velocity_interpolation_values'])
        for i in range(mesh.nNodes_global):
            for j in range(p.nd):
                fvel.write("  %f  " % velocityFEFunction.dof[i*p.nd+j])
            #mwf PT123 doc incorrect, only reads up to neq
            #for j in range(p.nd,3):
            #    fvel.write("  %g  " % 0.0)
            fvel.write("\n")
    #
    fvel.write("ENDR")
    fvel.close()
def writePT123elementVelocity(opts,p,mesh,params,tnList,q,velocitySpace,velocityFEFunction,
                              filebase="test_pt123"):
    #write out velocity file (have to loop through time steps)
    vel_prefix = 've1';
    if p.nd == 2:
        vel_prefix = 've2';
    elif p.nd == 3:
        vel_prefix = 've3';

    fvel = open(filebase+'.'+vel_prefix,'w')
    fvel.write("     %d     %d    %d \n" % (mesh.nElements_global,p.nd,len(tnList)))

    for i,t in enumerate(tnList):
        fvel.write("TS    %f \n" % t)

        try:
            p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
        except TypeError:
            for eN in range(mesh.nElements_global):
                for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                    q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],t)
        #evaluate velocity degrees of freedom from its interpolation conditions
        velocityFEFunction.projectFromInterpolationConditions(q['velocity_interpolation_values'])
        for eN in range(mesh.nElements_global):
            for nN in range(mesh.nNodes_element):
                I = velocitySpace.dofMap.l2g[eN,nN]
                for j in range(p.nd):
                    fvel.write("  %f  " % velocityFEFunction.dof[I*p.nd+j])
                #mwf PT123 doc incorrect, only reads up to neq
                #for j in range(p.nd,3):
                #    fvel.write("  %g " % 0.0)
            fvel.write("\n")
    #
    fvel.write("ENDR")
    fvel.close()
def writePT123RT0velocityAsElementVelocity(opts,p,mesh,params,tnList,q,ebq,
                                           filebase="test_pt123"):
    #write out velocity file (have to loop through time steps)
    vel_prefix = 've1';
    if p.nd == 2:
        vel_prefix = 've2';
    elif p.nd == 3:
        vel_prefix = 've3';

    fvel = open(filebase+'.'+vel_prefix,'w')
    fvel.write("     %d     %d    %d \n" % (mesh.nElements_global,p.nd,len(tnList)))

    local_nodes_array = numpy.zeros((mesh.nNodes_element,3),'d')
    local_velocity_array = numpy.zeros((mesh.nNodes_element,p.nd),'d')
    local_element_array  = numpy.zeros((mesh.nNodes_element,),'i')
    #mwf debug
    #import pdb
    #pdb.set_trace()
    for i,t in enumerate(tnList):
        fvel.write("TS    %f \n" % t)

        evaluateVelocity_RT0(opts,p,mesh,t,ebq,q['velocity_new_dof'],q['velocity_l2g'],q['velocity_interpolation_values'])

        #evaluate RT0 velocity at nodes
        from proteus import cpostprocessing
        for eN in range(mesh.nElements_global):
            local_element_array.fill(eN)
            local_velocity_array.fill(0.0)
            local_nodes_array.fill(0.0)
            for nN in range(mesh.nNodes_element):
                local_nodes_array[nN,:] = mesh.nodeArray[mesh.elementNodesArray[eN,nN],:]
            #
            cpostprocessing.getRT0velocityValuesFluxRep_arbitraryElementMembership(mesh.nodeArray,
                                                                                   mesh.elementNodesArray,
                                                                                   q['abs(det(J))'],
                                                                                   local_nodes_array,
                                                                                   local_element_array,
                                                                                   q['velocity_new_dof'],
                                                                                   local_velocity_array)
            for nN in range(mesh.nNodes_element):
                for j in range(p.nd):
                    fvel.write("  %f  " % local_velocity_array[nN,j])
                #mwf PT123 doc incorrect, only reads up to neq
                #for j in range(p.nd,3):
                #    fvel.write("  %g " % 0.0)
            fvel.write("\n")
    #
    fvel.write("ENDR")
    fvel.close()

def writePT123nodalVolumeFraction(opts,p,mesh,params,tnList,q,
                                  filebase="test_pt123"):
    #write out element volume fraction file (have to loop through time steps)
    emc_prefix = 'nemc%d' % p.nd;

    femc = open(filebase+'.'+emc_prefix,'w')
    femc.write("     %d     %d  \n" % (mesh.nNodes_global,len(tnList))) #constant for now

    for i,t in enumerate(tnList):
        femc.write("TS    %f \n" % t)

        for eN in range(mesh.nNodes_global):
            femc.write("%f \n" % 1.0)
    #
    femc.write("ENDR")
    femc.close()

def writePT123elementVolumeFraction(opts,p,mesh,params,tnList,q,
                              filebase="test_pt123"):
    #write out element volume fraction file (have to loop through time steps)
    emc_prefix = 'nemc%d' % p.nd;

    femc = open(filebase+'.'+emc_prefix,'w')
    femc.write("     %d     %d  \n" % (mesh.nElements_global,len(tnList))) #constant for now

    for i,t in enumerate(tnList):
        femc.write("TS    %f \n" % t)

        for eN in range(mesh.nElements_global):
            femc.write("%f \n" % 1.0)
    #
    femc.write("ENDR")
    femc.close()

def writePT123particleFile(opts,p,t,mesh,params,tnList,q,analyticalTracking,filebase="test_pt123",base=1):

    part_suffix = "pt%d" % p.nd
    fpart = open(filebase+"."+part_suffix,'w')
    if analyticalTracking:
        fpart.write("-1    ID_RK\n")
    else:
        fpart.write("45    ID_RK\n")

    fpart.write("%d  NP\n" % q['x_depart_start'].shape[0])
    for i in range(q['x_depart_start'].shape[0]):
        fpart.write("%d  %g  %g  %g \n" % (q['element_depart_start'][i]+base,
                                          q['x_depart_start'][i,0],
                                          q['x_depart_start'][i,1],
                                          q['x_depart_start'][i,2]))

    #
    ibf = 1 #forward
    if tnList[-1] < tnList[0]:
        ibf = -1
    fpart.write("%d   IBF (1 = forward 2 = backward) \n" % ibf)
    fpart.write("%f   T_START \n" % t)
    fpart.write("%f  %f  %d     DT_PT, DT_INPUT, ID_DT\n" % (tnList[-1],
                                                             params[('dt_init',0)],
                                                             params[('dt_init_flag',0)]))
    fpart.write("0   NT_PT_OUTPUT \n") #output every step for now
    fpart.write("%g %g %g     ATOL,RTOL,SF \n" % (params[('atol_tracking',0)],params[('rtol_tracking',0)],params[('sf_tracking',0)]))

    fpart.write("ENDR\n")
    fpart.close()

    #write sup file for post processing
    fsup = open(filebase+"_pv.sup",'w')
    fsup.write("%d  %d  %d    NPT  NEQ  NTSTEP\n" % (q['x_depart'].shape[0],p.nd,len(tnList)))
    for t in tnList:
        fsup.write("%f \n" % t)
    #
    fsup.close()
########################################################################
#begin actual tests
########################################################################

def test0(opts):
    """
    1d, constant velocity in space and time
    """
    p = __import__(opts.problem[:-3])

    assert p.nd == 1, "1d only for now"
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    assert opts.nnq > 0, "nnq = %s not ok must be > 0" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    mesh = setupMesh_1d(opts,p)
    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,opts.t_start)

    #velocity representation
    #P^1 in 1d
    velocitySpace,velocity,q['velocity_interpolation_values'] = setupVelocity_C0P1(opts,p,sdmesh,opts.t_start)

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    setupTrackingDataArrays(opts.t_start,opts.t_target,sdmesh,q,nq_per_element)

    #tracking-routine specific options
    #set options controlling tracking strategy
    skipPointsWithZeroSolution = False

    #epsilon
    zeroTol = 1.0e-5
    #which way to track (1. forward, -1. backward)
    direction = 1.0
    if opts.t_target < opts.t_start:
        direction = -1
    ctracking.trackPointsConstVelocity1d(sdmesh.nElements_global,                #sdmesh representation
                                         sdmesh.nNodes_global,
                                         sdmesh.nNodes_element,
                                         sdmesh.nElementBoundaries_element,
                                         sdmesh.nodeArray,
                                         sdmesh.elementNodesArray,
                                         sdmesh.elementNeighborsArray,
                                         velocity.dof[0],                     #constant velocity tracking
                                         direction,
                                         q['t_depart'],                       #point departure times
                                         q['t_track'],                        #target end time
                                         sdmesh.nElements_global*nq_per_element,#total number of points
                                         zeroTol,                             #for zero tests
                                         q['x'],                              #departure points
                                         q['element_track'],                  #in/out element locations
                                         q['x_track'],                        #arrival points
                                         q['flag_track'])                     #tracking status

    #mwf debug
    #pdb.set_trace()

    #exact or reference solution
    q['x_track_exact']= numpy.copy(q['x'])
    q['t_track_exact']= numpy.copy(q['t_depart'])
    velocity_exact = velocity.dof[0]
    #constant velocity
    for eN in range(sdmesh.nElements_global):
        for k in range(nq_per_element):
            dt = opts.t_target - q['t_track_exact'][eN,k]
            dx = velocity_exact*dt
            if q['x_track_exact'][eN,k,0] + dx > p.L[0]:
                dx = p.L[0]-q['x_track_exact'][eN,k,0]
                dt = dx/velocity_exact
            elif q['x_track_exact'][eN,k,0] + dx < 0.0:
                dx = 0.0 - q['x_track_exact'][eN,k,0]
                dt = dx/velocity_exact
            q['t_track_exact'][eN,k]   += dt
            q['x_track_exact'][eN,k,0] += dx
    #crude output
    fdep = open('x_depart.grf','w')
    farr = open('x_arrive.grf','w')
    fflg = open('x_flag.grf','w')
    ferr = open('x_err.grf','w')
    fexa = open('x_exa.grf','w')
    l2errX= 0.0; l2errT= 0.0;
    for eN in range(sdmesh.nElements_global):
        for k in range(nq_per_element):
            fdep.write("%12.5e %12.5e \n " % (q['x'][eN,k,0],q['t_depart'][eN,k]))
            farr.write("%12.5e %12.5e \n " % (q['x_track'][eN,k,0],q['t_track'][eN,k]))
            fflg.write("%12.5e %d     \n " % (q['x_track'][eN,k,0],q['flag_track'][eN,k]))
            ferr.write("%12.5e %12.5e \n"  % (q['x_track'][eN,k,0],q['x_track'][eN,k,0]-q['x_track_exact'][eN,k,0]))
            fexa.write("%12.5e %12.5e \n"  % (q['x_track_exact'][eN,k,0],q['t_track_exact'][eN,k]))
            l2errX += (q['x_track'][eN,k,0]-q['x_track_exact'][eN,k,0])**2
            l2errT += (q['t_track'][eN,k]-q['t_track_exact'][eN,k])**2
    #
    fdep.close()
    farr.close()
    fflg.close()
    ferr.close()
    #
    l2errX = math.sqrt(l2errX)/(sdmesh.nElements_global*nq_per_element)
    l2errT = math.sqrt(l2errT)/(sdmesh.nElements_global*nq_per_element)
    print "T= %s ell_2 error in spatial locations = %s ell_2 error in temporal locations = %s " % (opts.t_target,l2errX,l2errT)

def test1(opts):
    """
    1d, steady, linear velocity in space
    """
    p = __import__(opts.problem[:-3])

    assert p.nd == 1, "1d only for now"
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    assert opts.nnq > 0, "nnq = %s not ok must be > 1" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    mesh = setupMesh_1d(opts,p)
    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,opts.t_start)

    #velocity representation
    #P^1 in 1d
    velocitySpace,velocity,q['velocity_interpolation_values'] = setupVelocity_C0P1(opts,p,sdmesh,opts.t_start)

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    setupTrackingDataArrays(opts.t_start,opts.t_target,sdmesh,q,nq_per_element)

    #tracking-routine specific options
    #set options controlling tracking strategy
    skipPointsWithZeroSolution = False

    #epsilon
    zeroTol = 1.0e-5
    #which way to track (1. forward, -1. backward)
    direction = 1.0
    if opts.t_target < opts.t_start:
        direction = -1
    ctracking.trackPointsC0P1Velocity1d(sdmesh.nElements_global,                #mesh representation
                                        sdmesh.nNodes_global,
                                        sdmesh.nNodes_element,
                                        sdmesh.nElementBoundaries_element,
                                        sdmesh.nodeArray,
                                        sdmesh.elementNodesArray,
                                        sdmesh.elementNeighborsArray,
                                        velocitySpace.dofMap.l2g,
                                        velocity.dof,                       #constant velocity tracking
                                        direction,
                                        q['t_depart'],                       #point departure times
                                        q['t_track'],                        #target end time
                                        sdmesh.nElements_global*nq_per_element,#total number of points
                                        zeroTol,                             #for zero tests
                                        q['x'],                              #departure points
                                        q['element_track'],                  #in/out element locations
                                        q['x_track'],                        #arrival points
                                        q['flag_track'])                     #tracking status

    #mwf debug
    #pdb.set_trace()

    #exact or reference solution
    q['x_track_exact']= numpy.copy(q['x'])
    q['t_track_exact']= numpy.copy(q['t_depart'])
    velocity_exact = velocity.dof[0]
    if "analyticalSolutionParticleTrajectory" in dir(p):
        for eN in range(mesh.nElements_global):
            for k in range(nq_per_element):
                dt = opts.t_target - q['t_track_exact'][eN,k]
                x = p.analyticalSolutionParticleTrajectory[0].uOfXT(q['x_track_exact'][eN,k],dt)
                if x > p.L[0]:
                    dx = p.L[0]-q['x_track_exact'][eN,k,0]
                    dt = p.analyticalSolutionParticleTrajectory[0].dtOfdX(q['x_track_exact'][eN,k],
                                                                  numpy.array([dx,0.0,0.0],dtype='d'))
                    x = p.L[0]
                elif x < 0.0:
                    dx = 0.0 - q['x_track_exact'][eN,k,0]
                    dt = p.analyticalSolutionParticleTrajectory[0].dtOfdX(q['x_track_exact'][eN,k],
                                                                  numpy.array([dx,0.0,0.0],dtype='d'))
                    x = 0.0
                #mwf debug
                if dt == None:
                    pdb.set_trace()
                q['t_track_exact'][eN,k]   += dt
                q['x_track_exact'][eN,k,0] = x
    else:
        #constant velocity
        for eN in range(sdmesh.nElements_global):
            for k in range(nq_per_element):
                dt = opts.t_target - q['t_track_exact'][eN,k]
                dx = velocity_exact*dt
                if q['x_track_exact'][eN,k,0] + dx > p.L[0]:
                    dx = p.L[0]-q['x_track_exact'][eN,k,0]
                    dt = dx/velocity_exact
                elif q['x_track_exact'][eN,k,0] + dx < 0.0:
                    dx = 0.0 - q['x_track_exact'][eN,k,0]
                    dt = dx/velocity_exact
                q['t_track_exact'][eN,k]   += dt
                q['x_track_exact'][eN,k,0] += dx
    #crude output
    fdep = open('x_depart.grf','w')
    farr = open('x_arrive.grf','w')
    fflg = open('x_flag.grf','w')
    ferr = open('x_err.grf','w')
    fexa = open('x_exa.grf','w')
    l2errX= 0.0; l2errT= 0.0;
    for eN in range(sdmesh.nElements_global):
        for k in range(nq_per_element):
            fdep.write("%12.5e %12.5e \n " % (q['x'][eN,k,0],q['t_depart'][eN,k]))
            farr.write("%12.5e %12.5e \n " % (q['x_track'][eN,k,0],q['t_track'][eN,k]))
            fflg.write("%12.5e %d     \n " % (q['x_track'][eN,k,0],q['flag_track'][eN,k]))
            ferr.write("%12.5e %12.5e \n"  % (q['x_track'][eN,k,0],q['x_track'][eN,k,0]-q['x_track_exact'][eN,k,0]))
            fexa.write("%12.5e %12.5e \n"  % (q['x_track_exact'][eN,k,0],q['t_track_exact'][eN,k]))
            l2errX += (q['x_track'][eN,k,0]-q['x_track_exact'][eN,k,0])**2
            l2errT += (q['t_track'][eN,k]-q['t_track_exact'][eN,k])**2
    #
    fdep.close()
    farr.close()
    fflg.close()
    ferr.close()
    #
    l2errX = math.sqrt(l2errX)/(sdmesh.nElements_global*nq_per_element)
    l2errT = math.sqrt(l2errT)/(sdmesh.nElements_global*nq_per_element)
    print "T= %s ell_2 error in spatial locations = %s ell_2 error in temporal locations = %s " % (opts.t_target,l2errX,l2errT)

def test2(opts):
    """
    2d, steady, linear velocity in space
    """
    p = __import__(opts.problem[:-3])

    assert p.nd == 2, "2d only for now"
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    assert opts.nnq > 0, "nnq = %s not ok must be > 0" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    mesh = setupMesh_2d(opts,p)
    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,opts.t_start)


    #velocity representation
    #RT0 in 2d
    #since don't have a full FEM representation for vector spaces this requires more work
    q['velocity_dof'],q['velocity_l2g'],q['velocity_interpolation_values'],ebq,q['elementBoundaryQuadraturePoints'],q['elementBoundaryQuadratureWeights'] = setupVelocity_RT0(opts,p,sdmesh,trialSpace,opts.t_start,nd=2)
    #uses flux rep
    localVelocityRepresentationFlag = 2

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    setupTrackingDataArrays(opts.t_start,opts.t_target,sdmesh,q,nq_per_element)

    #save velocity field for visualization
    q['J']            = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['inverse(J)']   = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['det(J)']       = numpy.zeros((q['x'].shape[0],q['x'].shape[1]),'d')
    trialSpace.elementMaps.getJacobianValues(q['x_hat'],q['J'],q['inverse(J)'],q['det(J)'])

    q['abs(det(J))']  = numpy.absolute(q['det(J)'])
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    from proteus import cpostprocessing
    cpostprocessing.getRT0velocityValuesFluxRep_arbitraryElementMembership(sdmesh.nodeArray,
                                                                           sdmesh.elementNodesArray,
                                                                           q['abs(det(J))'],
                                                                           q['x'],
                                                                           q['element_track'],
                                                                           q['velocity_dof'],
                                                                           q['velocity_depart'])


    ########## visualization and output stuff
    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","test2",useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    if opts.view:
        try:  #run time visualization?
            from proteus import Viewers
            from proteusGraphical import vtkViewers
            Viewers.viewerOn("test2",'vtk')
            viewer = Viewers.V_base(p)
        except:
            pass
    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x'],opts.t_start,init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['x_track'].flat[:] = q['x']; q['t_viz'] = numpy.copy(q['t_depart'])
    q['x_depart'] = numpy.copy(q['x'])
    #for exact or reference solution
    q['x_track_exact']= numpy.copy(q['x_depart'])
    q['t_track_exact']= numpy.copy(q['t_depart'])
    if opts.view:
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray,
                                                vectors=q['velocity_depart'],viewTypes=['spheres'])#viewTypes=['spheres','streamlines'])
            if opts.wait:
                raw_input("Hit any key to continue")
        except:
            pass



    #2d need outer normals on each face
    #right now have element boundary quadrature array calculated values already so can grab 1st one
    #since its affine geometry
    #In actual code can skip most of the ebq info by just using
    #  \vec n = (\ten{J}^{-T} \vec \hat{n})/||\ten{J}^{-T} \vec \hat{n})||

    elementBoundaryOuterNormals = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,p.nd),'d')
    elementBoundaryOuterNormalsOrig = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,p.nd),'d')
    #mwf orig
    for eN in range(mesh.nElements_global):
        for ebN in range(mesh.nElementBoundaries_element):
            elementBoundaryOuterNormalsOrig[eN,ebN,:] = ebq['n'][eN,ebN,0,:]
    boundaryNormals = numpy.array(trialSpace.elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
    ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                            q['inverse(J)'],
                                            elementBoundaryOuterNormals)
    assert numpy.max(numpy.absolute(elementBoundaryOuterNormalsOrig.flat[:]-elementBoundaryOuterNormals.flat[:])) < 1.0e-8
    #tracking-routine specific options
    #set options controlling tracking strategy
    skipPointsWithZeroSolution = False

    #epsilon
    zeroTol = 1.0e-5
    #which way to track (1. forward, -1. backward)
    direction = 1.0
    if opts.t_target < opts.t_start:
        direction = -1


    #
    t     = opts.t_start
    dtout = (opts.t_target - opts.t_start)/opts.nnt
    for i in range(int(opts.nnt)):
        #update for next step
        t += dtout

        q['t_track'].fill(t)
        ctracking.trackPointsRT0Velocity2d(localVelocityRepresentationFlag,
                                           sdmesh.nElements_global,                #mesh representation
                                           sdmesh.nNodes_global,
                                           sdmesh.nNodes_element,
                                           sdmesh.nElementBoundaries_element,
                                           sdmesh.nodeArray,
                                           sdmesh.elementNodesArray,
                                           sdmesh.elementNeighborsArray,
                                           sdmesh.elementBoundariesArray,
                                           sdmesh.elementBoundaryBarycentersArray,
                                           elementBoundaryOuterNormals,
                                           q['velocity_l2g'],
                                           q['velocity_dof'],                       #constant velocity tracking
                                           direction,
                                           q['t_depart'],                       #point departure times
                                           q['t_track'],                        #target end time
                                           sdmesh.nElements_global*nq_per_element,#total number of points
                                           zeroTol,                             #for zero tests
                                           q['x'],                              #departure points
                                           q['element_track'],                  #in/out element locations
                                           q['x_track'],                        #arrival points
                                           q['flag_track'])                     #tracking status


        #evaluate velocity at tracked points too
        cpostprocessing.getRT0velocityValuesFluxRep_arbitraryElementMembership(sdmesh.nodeArray,
                                                                               sdmesh.elementNodesArray,
                                                                               q['abs(det(J))'],
                                                                               q['x_track'],
                                                                               q['element_track'],
                                                                               q['velocity_dof'],
                                                                               q['velocity_track'])






        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)

        if opts.view:
            q['t_viz'].flat[:] = q['t_track']
            try: #run time visualization?
                vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
                if opts.wait:
                    raw_input("Hit any key to continue")
            except:
                pass

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x'].flat[:] = q['x_track'].flat
        #can also set flag so that don't track points that exited the domain
        #q['flag_track'][numpy.where(q['flag_track'] == 1)] = -2
        #mwf debug
        #pdb.set_trace()
    #output step loop
    #exact or reference solution
    if "analyticalSolutionParticleTrajectory" in dir(p):
        for eN in range(mesh.nElements_global):
            for k in range(nq_per_element):
                dt = t - q['t_track_exact'][eN,k]
                x0= q['x_track_exact'][eN,k]
                x = p.analyticalSolutionParticleTrajectory[0].uOfXT(x0,dt)
#                 #check distances to the boundaries
#                 if x[0] < 0.0:
#                     n  = numpy.array([-1.0,0.0]); xb = numpy.array([0.0,0.5])
#                     dx = xb-x0[:p.nd]
#                     dtb= p.analyticalSolutionParticleTrajectory[0].dtOfdX(x0,
#                                                                           dx,n)
#                     if abs(dtb) < abs(dt):
#                         dt = dtb
#                 if x[0] > p.L[0]:
#                     n  = numpy.array([1.0,0.0]); xb = numpy.array([p.L[0],0.5])
#                     dx = xb-x0[:p.nd]
#                     dtb= p.analyticalSolutionParticleTrajectory[0].dtOfdX(x0,
#                                                                           dx,n)
#                     if abs(dtb) < abs(dt):
#                         dt = dtb

#                 if x[1] < 0.0:
#                     n  = numpy.array([0.0,-1.0]); xb = numpy.array([0.5,0.0])
#                     dx = xb-x0[:p.nd]
#                     dtb= p.analyticalSolutionParticleTrajectory[0].dtOfdX(x0,
#                                                                           dx,n)
#                     if abs(dtb) < abs(dt):
#                         dt = dtb

#                 if x[1] > p.L[1]:
#                     n  = numpy.array([0.0,1.0]); xb = numpy.array([0.5,p.L[1]])
#                     dx = xb-x0[:p.nd]
#                     dtb= p.analyticalSolutionParticleTrajectory[0].dtOfdX(x0,
#                                                                           dx,n)
#                     if abs(dtb) < abs(dt):
#                         dt = dtb

                x = p.analyticalSolutionParticleTrajectory[0].uOfXT(x0,dt)
                #mwf debug
                if dt == None:
                    pdb.set_trace()
                q['t_track_exact'][eN,k]          += dt
                q['x_track_exact'][eN,k][:p.nd]    = x

    writerExact = Archiver.XdmfWriter()
    arGridExact = writerExact.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track_exact'],t,init=True,tCount=0,spaceSuffix="_particles_exact")
    writerExact.writeScalarXdmf_particles(archive,q['t_track_exact'],"t_particle_exact",tCount=0)

    #crude error calculation
    fdep = open('x_depart.grf','w')
    farr = open('x_arrive.grf','w')
    fflg = open('x_flag.grf','w')
    ferr = open('x_err.grf','w')
    fexa = open('x_exa.grf','w')
    l2errX= 0.0; l2errT= 0.0;
    for eN in range(sdmesh.nElements_global):
        for k in range(nq_per_element):
            fdep.write("%12.5e %12.5e %12.5e \n " % (q['x'][eN,k,0],q['x'][eN,k,1],q['t_depart'][eN,k]))
            farr.write("%12.5e %12.5e %12.5e \n " % (q['x_track'][eN,k,0],q['x_track'][eN,k,1],q['t_track'][eN,k]))
            fflg.write("%12.5e %12.5e %d     \n " % (q['x_track'][eN,k,0],q['x_track'][eN,k,1],q['flag_track'][eN,k]))
            ferr.write("%12.5e %12.5e %12.5e \n"  % (q['x_track'][eN,k,0],q['x_track'][eN,k,1],numpy.sqrt(numpy.dot(q['x_track'][eN,k]-q['x_track_exact'][eN,k],
                                                                                                                    q['x_track'][eN,k]-q['x_track_exact'][eN,k]))))
            fexa.write("%12.5e %12.5e %12.5e \n"  % (q['x_track_exact'][eN,k,0],q['x_track_exact'][eN,k,1],q['t_track_exact'][eN,k]))
            l2errX += (q['x_track'][eN,k,0]-q['x_track_exact'][eN,k,0])**2
            l2errT += (q['t_track'][eN,k]-q['t_track_exact'][eN,k])**2
    #
    fdep.close()
    farr.close()
    fflg.close()
    ferr.close()
    #
    l2errX = math.sqrt(l2errX)/(sdmesh.nElements_global*nq_per_element)
    l2errT = math.sqrt(l2errT)/(sdmesh.nElements_global*nq_per_element)
    print "T= %s ell_2 error in spatial locations = %s ell_2 error in temporal locations = %s " % (opts.t_target,l2errX,l2errT)
    #accommodate vtk updating paradigm
    if opts.view:
        q['t_viz'].flat[:] = q['t_track']
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
        except:
            pass

    #close out archives and visualization
    archive.close()

    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass

def test3(opts):
    """
    2d, steady, linear velocity in space, check error by tracking backwards at end of time
    """
    p = __import__(opts.problem[:-3])

    assert p.nd == 2, "2d only for now"
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    assert opts.nnq > 0, "nnq = %s not ok must be > 0" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    mesh = setupMesh_2d(opts,p)
    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,opts.t_start)


    #velocity representation
    #RT0 in 2d
    #since don't have a full FEM representation for vector spaces this requires more work
    q['velocity_dof'],q['velocity_l2g'],q['velocity_interpolation_values'],ebq,q['elementBoundaryQuadraturePoints'],q['elementBoundaryQuadratureWeights'] = setupVelocity_RT0(opts,p,sdmesh,trialSpace,opts.t_start,nd=2)

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    setupTrackingDataArrays(opts.t_start,opts.t_target,sdmesh,q,nq_per_element)

    #save velocity field for visualization
    q['J']            = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['inverse(J)']   = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['det(J)']       = numpy.zeros((q['x'].shape[0],q['x'].shape[1]),'d')
    trialSpace.elementMaps.getJacobianValues(q['x_hat'],q['J'],q['inverse(J)'],q['det(J)'])

    q['abs(det(J))']  = numpy.absolute(q['det(J)'])
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    from proteus import cpostprocessing
    cpostprocessing.getRT0velocityValuesFluxRep_arbitraryElementMembership(sdmesh.nodeArray,
                                                                           sdmesh.elementNodesArray,
                                                                           q['abs(det(J))'],
                                                                           q['x'],
                                                                           q['element_track'],
                                                                           q['velocity_dof'],
                                                                           q['velocity_depart'])


    ########## visualization and output stuff
    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","test3",useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    if opts.view:
        try:  #run time visualization?
            from proteus import Viewers
            from proteusGraphical import vtkViewers
            Viewers.viewerOn("test2",'vtk')
            viewer = Viewers.V_base(p)
        except:
            pass
    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x'],opts.t_start,init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['x_track'].flat[:] = q['x']; q['t_viz'] = numpy.copy(q['t_depart'])
    q['x_depart'] = numpy.copy(q['x'])
    #for exact or reference solution
    q['x_track_exact']= numpy.copy(q['x_depart'])
    q['t_track_exact']= numpy.copy(q['t_depart'])
    if opts.view:
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray,
                                                vectors=q['velocity_depart'],viewTypes=['spheres'])#viewTypes=['spheres','streamlines'])
            if opts.wait:
                raw_input("Hit any key to continue")
        except:
            pass



    #2d need outer normals on each face
    #right now have element boundary quadrature array calculated values already so can grab 1st one
    #since its affine geometry
    #In actual code can skip most of the ebq info by just using
    #  \vec n = (\ten{J}^{-T} \vec \hat{n})/||\ten{J}^{-T} \vec \hat{n})||

    elementBoundaryOuterNormals = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,p.nd),'d')
    elementBoundaryOuterNormalsOrig = numpy.zeros((mesh.nElements_global,mesh.nElementBoundaries_element,p.nd),'d')
    #mwf orig
    #\todo move out of python
    for eN in range(mesh.nElements_global):
        for ebN in range(mesh.nElementBoundaries_element):
            elementBoundaryOuterNormalsOrig[eN,ebN,:] = ebq['n'][eN,ebN,0,:]
    boundaryNormals = numpy.array(trialSpace.elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
    ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                            q['inverse(J)'],
                                            elementBoundaryOuterNormals)
    assert numpy.max(numpy.absolute(elementBoundaryOuterNormalsOrig.flat[:]-elementBoundaryOuterNormals.flat[:])) < 1.0e-8
    #tracking-routine specific options
    #set options controlling tracking strategy
    skipPointsWithZeroSolution = False

    #epsilon
    zeroTol = 1.0e-5
    #which way to track (1. forward, -1. backward)
    direction = 1.0
    if opts.t_target < opts.t_start:
        direction = -1


    #
    t     = opts.t_start
    dtout = (opts.t_target - opts.t_start)/opts.nnt
    for i in range(int(opts.nnt)):
        #update for next step
        t += dtout

        q['t_track'].fill(t)
        ctracking.trackPointsRT0Velocity2d(sdmesh.nElements_global,                #mesh representation
                                           sdmesh.nNodes_global,
                                           sdmesh.nNodes_element,
                                           sdmesh.nElementBoundaries_element,
                                           sdmesh.nodeArray,
                                           sdmesh.elementNodesArray,
                                           sdmesh.elementNeighborsArray,
                                           sdmesh.elementBoundariesArray,
                                           sdmesh.elementBoundaryBarycentersArray,
                                           elementBoundaryOuterNormals,
                                           q['velocity_l2g'],
                                           q['velocity_dof'],                       #constant velocity tracking
                                           direction,
                                           q['t_depart'],                       #point departure times
                                           q['t_track'],                        #target end time
                                           sdmesh.nElements_global*nq_per_element,#total number of points
                                           zeroTol,                             #for zero tests
                                           q['x'],                              #departure points
                                           q['element_track'],                  #in/out element locations
                                           q['x_track'],                        #arrival points
                                           q['flag_track'])                     #tracking status


        #evaluate velocity at tracked points too
        cpostprocessing.getRT0velocityValuesFluxRep_arbitraryElementMembership(sdmesh.nodeArray,
                                                                               sdmesh.elementNodesArray,
                                                                               q['abs(det(J))'],
                                                                               q['x_track'],
                                                                               q['element_track'],
                                                                               q['velocity_dof'],
                                                                               q['velocity_track'])






        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)

        if opts.view:
            q['t_viz'].flat[:] = q['t_track']
            try: #run time visualization?
                vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
                if opts.wait:
                    raw_input("Hit any key to continue")
            except:
                pass

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x'].flat[:] = q['x_track'].flat
        #can also set flag so that don't track points that exited the domain
        #q['flag_track'][numpy.where(q['flag_track'] == 1)] = -2
        #mwf debug
        #pdb.set_trace()
    #output step loop

    #first just test stepping forward with reversed velocity
    #direction = -1
    q['velocity_dof'] *= -1.0
    toutput = t
    for i in range(int(opts.nnt),2*(int(opts.nnt))):
        #update for next step
        toutput += dtout
        t       += dtout#t       -= dtout
        q['t_track'].fill(t)#q['t_track'] = numpy.minimum(q['t_track'],t)
        ctracking.trackPointsRT0Velocity2d(sdmesh.nElements_global,                #mesh representation
                                           sdmesh.nNodes_global,
                                           sdmesh.nNodes_element,
                                           sdmesh.nElementBoundaries_element,
                                           sdmesh.nodeArray,
                                           sdmesh.elementNodesArray,
                                           sdmesh.elementNeighborsArray,
                                           sdmesh.elementBoundariesArray,
                                           sdmesh.elementBoundaryBarycentersArray,
                                           elementBoundaryOuterNormals,
                                           q['velocity_l2g'],
                                           q['velocity_dof'],                       #constant velocity tracking
                                           direction,
                                           q['t_depart'],                       #point departure times
                                           q['t_track'],                        #target end time
                                           sdmesh.nElements_global*nq_per_element,#total number of points
                                           zeroTol,                             #for zero tests
                                           q['x'],                              #departure points
                                           q['element_track'],                  #in/out element locations
                                           q['x_track'],                        #arrival points
                                           q['flag_track'])                     #tracking status


        #evaluate velocity at tracked points too
        cpostprocessing.getRT0velocityValuesFluxRep_arbitraryElementMembership(sdmesh.nodeArray,
                                                                               sdmesh.elementNodesArray,
                                                                               q['abs(det(J))'],
                                                                               q['x_track'],
                                                                               q['element_track'],
                                                                               q['velocity_dof'],
                                                                               q['velocity_track'])






        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],toutput,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)

        if opts.view:
            q['t_viz'].flat[:] = q['t_track']
            try: #run time visualization?
                vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
                if opts.wait:
                    raw_input("Hit any key to continue")
            except:
                pass

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x'].flat[:] = q['x_track'].flat
        #can also set flag so that don't track points that exited the domain
        #q['flag_track'][numpy.where(q['flag_track'] == 1)] = -2
        #mwf debug
        #pdb.set_trace()
    #output step loop

    writerExact = Archiver.XdmfWriter()
    arGridExact = writerExact.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track_exact'],toutput,init=True,tCount=0,spaceSuffix="_particles_exact")
    writerExact.writeScalarXdmf_particles(archive,q['t_track_exact'],"t_particle_exact",tCount=0)

    #\todo get rid of?
    #crude error calculation
    fdep = open('x_depart.grf','w')
    farr = open('x_arrive.grf','w')
    fflg = open('x_flag.grf','w')
    ferr = open('x_err.grf','w')
    fexa = open('x_exa.grf','w')
    l2errX= 0.0; l2errT= 0.0;
    for eN in range(sdmesh.nElements_global):
        for k in range(nq_per_element):
            fdep.write("%12.5e %12.5e %12.5e \n " % (q['x'][eN,k,0],q['x'][eN,k,1],q['t_depart'][eN,k]))
            farr.write("%12.5e %12.5e %12.5e \n " % (q['x_track'][eN,k,0],q['x_track'][eN,k,1],q['t_track'][eN,k]))
            fflg.write("%12.5e %12.5e %d     \n " % (q['x_track'][eN,k,0],q['x_track'][eN,k,1],q['flag_track'][eN,k]))
            ferr.write("%12.5e %12.5e %12.5e \n"  % (q['x_track'][eN,k,0],q['x_track'][eN,k,1],numpy.sqrt(numpy.dot(q['x_track'][eN,k]-q['x_track_exact'][eN,k],
                                                                                                                    q['x_track'][eN,k]-q['x_track_exact'][eN,k]))))
            fexa.write("%12.5e %12.5e %12.5e \n"  % (q['x_track_exact'][eN,k,0],q['x_track_exact'][eN,k,1],q['t_track_exact'][eN,k]))
            l2errX += (q['x_track'][eN,k,0]-q['x_track_exact'][eN,k,0])**2
            l2errT += (q['t_track'][eN,k]-q['t_track_exact'][eN,k])**2
    #
    fdep.close()
    farr.close()
    fflg.close()
    ferr.close()
    #
    l2errX = math.sqrt(l2errX)/(sdmesh.nElements_global*nq_per_element)
    l2errT = math.sqrt(l2errT)/(sdmesh.nElements_global*nq_per_element)
    print "T= %s ell_2 error in spatial locations = %s ell_2 error in temporal locations = %s " % (opts.t_target,l2errX,l2errT)
    #accommodate vtk updating paradigm
    if opts.view:
        q['t_viz'].flat[:] = q['t_track']
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
        except:
            pass

    #close out archives and visualization
    archive.close()

    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass

def test4(opts):
    """
    2d, steady, linear velocity in space
    """
    p = __import__(opts.problem[:-3])

    assert p.nd == 2, "2d only for now"
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    assert opts.nnq > 0, "nnq = %s not ok must be > 1" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    mesh = setupMesh_2d(opts,p)
    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,opts.t_start)

    #velocity representation
    #P^1 in 2d
    velocitySpace,velocity,q['velocity_interpolation_values'] = setupVelocity_C0P1(opts,p,sdmesh,opts.t_start)

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    setupTrackingDataArrays(opts.t_start,opts.t_target,sdmesh,q,nq_per_element)

    #save velocity field for visualization
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    velocity.getValues(q[('u_shape_functions')],q['velocity_depart'])
    q['velocity_track'] = q['velocity_depart'][:]
    particle_tracker = SteadyState_LinearAdvection_C0P1Velocity_PT123(sdmesh,p.nd,
                                                                      {0:velocitySpace.dofMap.l2g},
                                                                      {0:velocity.dof},
                                                                      activeComponentList=[0])

    particle_tracker.setFromOptions(opts)

    ########## visualization and output stuff
    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","test4",useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    if opts.view:
        try:  #run time visualization?
            from proteus import Viewers
            from proteusGraphical import vtkViewers
            Viewers.viewerOn("test4",'vtk')
            viewer = Viewers.V_base(p)
        except:
            pass
    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x'],opts.t_start,init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['x_track'].flat[:] = q['x']; q['t_viz'] = numpy.copy(q['t_depart'])
    q['x_depart'] = numpy.copy(q['x'])
    if opts.view:
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray,
                                                vectors=q['velocity_depart'],viewTypes=['spheres'])#viewTypes=['spheres','streamlines'])
            if opts.wait:
                raw_input("Hit any key to continue")
        except:
            pass

    #which way to track (1. forward, -1. backward)
    direction = 1.0
    if opts.t_target < opts.t_start:
        direction = -1


    #
    t     = opts.t_start
    dtout = (opts.t_target - opts.t_start)/opts.nnt
    for i in range(int(opts.nnt)):
        #update for next step
        t += dtout

        q['t_track'].fill(t)
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if direction > 0.0:
            particle_tracker.forwardTrack({0:q['t_depart']},
                                          {0:q['t_track']},                          #target end time
                                          {0:sdmesh.nElements_global*nq_per_element},#total number of points
                                          {0:q['x']},                                #departure points
                                          {0:q['element_track']},                    #in/out element locations
                                          {0:q['x_track']},                          #arrival points
                                          {0:q['flag_track']})
        else:
            particle_tracker.backwardTrack({0:q['t_depart']},
                                           {0:q['t_track']},                          #target end time
                                           {0:sdmesh.nElements_global*nq_per_element},#total number of points
                                           {0:q['x']},                                #departure points
                                           {0:q['element_track']},                    #in/out element locations
                                           {0:q['x_track']},                          #arrival points
                                           {0:q['flag_track']})

        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)

        if opts.view:
            q['t_viz'].flat[:] = q['t_track']
            try: #run time visualization?
                vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
                if opts.wait:
                    raw_input("Hit any key to continue")
            except:
                pass

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x'].flat[:] = q['x_track'].flat
        #can also set flag so that don't track points that exited the domain
        #q['flag_track'][numpy.where(q['flag_track'] == 1)] = -2
        #mwf debug
        #pdb.set_trace()
    #output step loop
    #mwf debug
    #pdb.set_trace()

    #accommodate vtk updating paradigm
    if opts.view:
        q['t_viz'].flat[:] = q['t_track']
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
        except:
            pass

    #close out archives and visualization
    archive.close()

    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass


def test5(opts):
    """
    transient, linear velocity in space. Computes error by comparing starting and ending location for points in a specified
      interval
    """
    p = __import__(opts.problem[:-3])

    assert p.nd in [1,2], "1d or 2d only for now"
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    assert opts.nnq > 0, "nnq = %s not ok must be > 1" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    if p.nd == 2:
        mesh = setupMesh_2d(opts,p)
    else:
        mesh = setupMesh_1d(opts,p)

    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,opts.t_start)

    #velocity representation
    #P^1 in 2d
    velocitySpace,velocity,q['velocity_interpolation_values'] = setupVelocity_C0P1(opts,p,sdmesh,opts.t_start)
    component_velocity_times={0:opts.t_start}

    velocity_new = FemTools.FiniteElementFunction(velocitySpace,dim_dof=p.nd,isVector=True,name="velocity_new")
    velocity_new.dof[:] = velocity.dof
    component_velocity_times_new={0:opts.t_start}

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    setupTrackingDataArrays(opts.t_start,opts.t_target,sdmesh,q,nq_per_element)

    #save velocity field for visualization
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    velocity.getValues(q[('u_shape_functions')],q['velocity_depart'])
    q['velocity_track'] = q['velocity_depart'][:]
    particle_tracker = LinearAdvection_C0P1Velocity_PT123(sdmesh,p.nd,
                                                          {0:velocitySpace.dofMap.l2g},
                                                          {0:velocity.dof},
                                                          activeComponentList=[0])

    particle_tracker.setFromOptions(opts)

    ########## visualization and output stuff
    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","test5",useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    if opts.view:
        try:  #run time visualization?
            from proteus import Viewers
            from proteusGraphical import vtkViewers
            Viewers.viewerOn("test5",'vtk')
            viewer = Viewers.V_base(p)
        except:
            pass
    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x'],opts.t_start,init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['x_track'].flat[:] = q['x']; q['t_viz'] = numpy.copy(q['t_depart'])
    q['x_depart'] = numpy.copy(q['x'])
    if opts.view:
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray,
                                                vectors=q['velocity_depart'],viewTypes=['spheres'])#viewTypes=['spheres','streamlines'])
            if opts.wait:
                raw_input("Hit any key to continue")
        except:
            pass

    #which way to track (1. forward, -1. backward)
    direction = 1.0
    if opts.t_target < opts.t_start:
        direction = -1

    #
    t     = opts.t_start
    dtout = (opts.t_target - opts.t_start)/opts.nnt
    for i in range(int(opts.nnt)):
        #update for next step
        velocity.dof[:]=velocity_new.dof[:]
        #don't need this actually because setup as a shallow copy?
        particle_tracker.setTrackingVelocity(velocity.dof,0,t,timeLevel=0)

        t += dtout

        #evaluate velocity
        #\todo move out of python
        try:
            p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
        except TypeError:
            for eN in range(mesh.nElements_global):
                for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                    q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],t)
        #mwf debug
        #import pdb
        #pdb.set_trace()
        #evaluate velocity degrees of freedom from its interpolation conditions
        velocity_new.projectFromInterpolationConditions(q['velocity_interpolation_values'])
        particle_tracker.setTrackingVelocity(velocity_new.dof,0,t,timeLevel=1)

        q['t_track'].fill(t)
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if direction > 0.0:
            particle_tracker.forwardTrack({0:q['t_depart']},
                                          {0:q['t_track']},                          #target end time
                                          {0:sdmesh.nElements_global*nq_per_element},#total number of points
                                          {0:q['x']},                                #departure points
                                          {0:q['element_track']},                    #in/out element locations
                                          {0:q['x_track']},                          #arrival points
                                          {0:q['flag_track']})
        else:
            particle_tracker.backwardTrack({0:q['t_depart']},
                                           {0:q['t_track']},                          #target end time
                                           {0:sdmesh.nElements_global*nq_per_element},#total number of points
                                           {0:q['x']},                                #departure points
                                           {0:q['element_track']},                    #in/out element locations
                                           {0:q['x_track']},                          #arrival points
                                           {0:q['flag_track']})

        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)

        if opts.view:
            q['t_viz'].flat[:] = q['t_track']
            try: #run time visualization?
                if p.nd == 2:
                    vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
                if opts.wait:
                    raw_input("Hit any key to continue")
            except:
                pass

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x'].flat[:] = q['x_track'].flat
        #can also set flag so that don't track points that exited the domain
        #q['flag_track'][numpy.where(q['flag_track'] == 1)] = -2
        #mwf debug
        #pdb.set_trace()
    #output step loop
    #mwf debug
    #pdb.set_trace()

    #accommodate vtk updating paradigm
    if opts.view:
        q['t_viz'].flat[:] = q['t_track']
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_2D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
        except:
            pass

    #close out archives and visualization
    archive.close()

    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass


def test6(opts):
    """
    transient, linear velocity in space. Computes error by comparing starting and ending location for points in a specified
      interval
    """
    p = __import__(opts.problem[:-3])

    assert p.nd in [1,2,3]
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    assert opts.nnq > 0, "nnq = %s not ok must be > 1" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    if p.nd == 3:
        mesh = setupMesh_3d(opts,p)
    elif p.nd == 2:
        mesh = setupMesh_2d(opts,p)
    else:
        mesh = setupMesh_1d(opts,p)

    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,opts.t_start)

    #velocity representation
    #
    velocitySpace,velocity,q['velocity_interpolation_values'] = setupVelocity_C0P1(opts,p,sdmesh,opts.t_start)
    component_velocity_times={0:opts.t_start}

    velocity_new = FemTools.FiniteElementFunction(velocitySpace,dim_dof=p.nd,isVector=True,name="velocity_new")
    velocity_new.dof[:] = velocity.dof
    component_velocity_times_new={0:opts.t_start}

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    setupTrackingDataArrays(opts.t_start,opts.t_target,sdmesh,q,nq_per_element)

    #save velocity field for visualization
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    velocity.getValues(q[('u_shape_functions')],q['velocity_depart'])
    q['velocity_track'] = q['velocity_depart'][:]
    particle_tracker = LinearAdvection_C0P1Velocity_PT123(sdmesh,p.nd,
                                                          {0:velocitySpace.dofMap.l2g},
                                                          {0:velocity.dof},
                                                          activeComponentList=[0])

    particle_tracker.setFromOptions(opts)

    ########## visualization and output stuff
    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","test6",useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    if opts.view:
        try:  #run time visualization?
            from proteus import Viewers
            from proteusGraphical import vtkViewers
            Viewers.viewerOn("test6",'vtk')
            viewer = Viewers.V_base(p)
        except:
            pass
    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x'],opts.t_start,init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['x_track'].flat[:] = q['x']; q['t_viz'] = numpy.copy(q['t_depart'])
    q['x_depart'] = numpy.copy(q['x'])
    if opts.view:
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_3D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray,
                                                vectors=q['velocity_depart'],viewTypes=['spheres'])#viewTypes=['spheres','streamlines'])
            if opts.wait:
                raw_input("Hit any key to continue")
        except:
            pass

    #which way to track (1. forward, -1. backward)
    direction = 1.0
    if opts.t_target < opts.t_start:
        direction = -1

    #
    t     = opts.t_start
    dtout = (opts.t_target - opts.t_start)/opts.nnt
    for i in range(int(opts.nnt)):
        #update for next step
        velocity.dof[:]=velocity_new.dof[:]
        #don't need this actually because setup as a shallow copy?
        particle_tracker.setTrackingVelocity(velocity.dof,0,t,timeLevel=0)

        t += dtout

        #evaluate velocity
        #\todo move out of python
        try:
            p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
        except TypeError:
            for eN in range(mesh.nElements_global):
                for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                    q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],t)
        #mwf debug
        #import pdb
        #pdb.set_trace()
        #evaluate velocity degrees of freedom from its interpolation conditions
        velocity_new.projectFromInterpolationConditions(q['velocity_interpolation_values'])
        particle_tracker.setTrackingVelocity(velocity_new.dof,0,t,timeLevel=1)

        q['t_track'].fill(t)
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if direction > 0.0:
            particle_tracker.forwardTrack({0:q['t_depart']},
                                          {0:q['t_track']},                          #target end time
                                          {0:sdmesh.nElements_global*nq_per_element},#total number of points
                                          {0:q['x']},                                #departure points
                                          {0:q['element_track']},                    #in/out element locations
                                          {0:q['x_track']},                          #arrival points
                                          {0:q['flag_track']})
        else:
            particle_tracker.backwardTrack({0:q['t_depart']},
                                           {0:q['t_track']},                          #target end time
                                           {0:sdmesh.nElements_global*nq_per_element},#total number of points
                                           {0:q['x']},                                #departure points
                                           {0:q['element_track']},                    #in/out element locations
                                           {0:q['x_track']},                          #arrival points
                                           {0:q['flag_track']})

        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)

        if opts.view:
            q['t_viz'].flat[:] = q['t_track']
            try: #run time visualization?
                if p.nd == 3:
                    vtkViewers.viewScalar_pointCloud_3D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
                if opts.wait:
                    raw_input("Hit any key to continue")
            except:
                pass

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x'].flat[:] = q['x_track'].flat
        #can also set flag so that don't track points that exited the domain
        #q['flag_track'][numpy.where(q['flag_track'] == 1)] = -2
        #mwf debug
        #pdb.set_trace()
    #output step loop
    #mwf debug
    #pdb.set_trace()

    #accommodate vtk updating paradigm
    if opts.view:
        q['t_viz'].flat[:] = q['t_track']
        try: #run time visualization?
            vtkViewers.viewScalar_pointCloud_3D(q['x_track'],q['t_viz'],"t",Viewers.windowNumber,sdmesh.nodeArray,sdmesh.elementNodesArray)
        except:
            pass

    #close out archives and visualization
    archive.close()

    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass


def Lu1Dex(opts,tryPT123A=True,writePT123files=True):
    """
    1D example from Lu 94. Single-phase flow in a semi-infinite 1d
    domain. Aquifer is b=10 m thick with Transmissivity T=20 cm^2/s,
    storage coefficient S=0.002, and porosity \theta=0.5. The aquifer
    conductivity and specific storage are

     K = T/b 2.0e-2 cm/s = 1728.0 cm/d, S_s = S/b = 2.0e-6 1/cm

     S_s \pd{h}{t} = K\pd2{h}{x}

     for
       x \in [0,\infty]

     Initial and boundary conditions are

       h(x,0) = H= 10 m

       h(0,t) = 0 m, t > 0

       \pd{h}{x} = 0 x \rightarrow \infty  (double check?)

     Analytical solution is

       h = H \erf(\frac{x}{2\sqrt{tK/S_s}}

      \pd{h}{x} = \frac{H}{\sqrt{\pi t K / S_s}}\exp\left(-\frac{S_s}{4 K t}x^2\right)   , note Lu has extra K in numerator here

      v_x = -\frac{K}{\theta}\pd{h}{x} = -\frac{HK}{\theta\sqrt{\pi t K / S_s}}\exp\left(-\frac{S_s}{4 K t}x^2\right)

     Dense grid solution for travel time to boundary of point placed at x=500.0, t=1 min is 14.35 [d], the Lu and Russell Healy answer is 13.28 [d]
        For a point placed at x=500, t=1000.0 m, the dense grid solution is 21.15 [d], and the Lu answer = 20.77 [d], Russell and Healy answer is 20.78 [d]

     The C0P1 answer on 07/13/10 is 13.2814 for the t=1 min particle and 20.7506153 [d] for the t=1000 [m] particle using atol=dn_safe = 1.0e-7, rtol=0
    """
    from math import pi,sqrt,exp
    #dummy p module
    class LuProblem:
        class LuIC:
            def __init__(self,H=1000.0):
                self.H=H
            def uOfXT(self,x,t):
                return self.H
        class LuParticleVelocity1d:
            def __init__(self,H=1000.0,K=2.0e-2,S=2.0e-6,theta=0.5):
                self.H=H; self.K=K; self.S=S; self.theta=theta
            def uOfXT(self,x,t):
                if t < 1.0e-12:
                    return 0.0
                term1 = -self.K*self.H/(self.theta*sqrt(pi*self.K*t/self.S))
                term2 = exp(-self.S*x[0]*x[0]/(4.0*self.K*t))
                return term1*term2
            def uOfXTv(self,x,t,v):
                v.fill(0.0)
                if t < 1.0e-12:
                    return
                for i in range(len(v.flat)):#1d
                    term1 = -self.K*self.H/(self.theta*sqrt(pi*self.K*t/self.S))
                    term2 = exp(-self.S*x.flat[i*3+0]*x.flat[i*3+0]/(4.0*self.K*t))
                    v.flat[i] = term1*term2

        def __init__(self):
            self.nd = 1
            #time intervals for velocity evaluation from Lu in days
            self.tnList_days=[0.00035,0.001,0.01,0.05,0.2,0.7,1.2,2.0,3.0,5.0,9.0,13.0,17.0,21.0,30.0]
            self.tnList = [t*60.0*60.0*24.0 for t in self.tnList_days]
            self.L=(500.0,1.0,1.0) #cm
            self.H= 1000.0 #cm
            self.K= 2.0e-2#1728.0 cm/d = 2.0e-2 #cm/s
            self.S = 2.0e-6
            self.theta=0.5
            self.initialConditions = {0:LuProblem.LuIC(H=self.H)}
            self.analyticalSolutionParticleVelocity={0:LuProblem.LuParticleVelocity1d(H=self.H,K=self.K,S=self.S,theta=self.theta)}
        def getInitialTrackingLocations(self,q):
            """
            return points that need to be tracked
            """
            nPoints = 1
            x = numpy.zeros((1,3),'d')
            x[0,0] = 500.0
            return x
    p = LuProblem()

    assert opts.nnq > 0, "nnq = %s not ok must be > 1" % opts.nnq

    tnList = p.tnList

    #to collect integration points etc
    q = {}
    mesh = setupMesh_1d(opts,p)

    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,tnList[0])

    #velocity representation
    #
    velocitySpace,velocity,q['velocity_interpolation_values'] = setupVelocity_C0P1(opts,p,sdmesh,tnList[0])

    component_velocity_times={0:tnList[0]}

    velocity_new = FemTools.FiniteElementFunction(velocitySpace,dim_dof=p.nd,isVector=True,name="velocity_new")
    velocity_new.dof[:] = velocity.dof
    component_velocity_times_new={0:tnList[0]}

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    #options for tracking
    particleTracking_params = {}
    particleTracking_params[('dn_safe_tracking',0)]=1.0e-7
    particleTracking_params[('temporalVariationFlag',0)]=2
    particleTracking_params[('atol_tracking',0)]=1.0e-7
    particleTracking_params[('rtol_tracking',0)]=0.0
    particleTracking_params[('sf_tracking',0)]=0.9
    if not tryPT123A or writePT123files:
        particleTracking_params[('dt_init',0)]= 0.001
        particleTracking_params[('dt_init_flag',0)] = 1 #how to pick initial time step (0 --> const, 1 --> CFL=1)
        particleTracking_params[('rk_flag',0)] = 45 #RK type

    opts.particleTracking_params = particleTracking_params

    #find the interval where start time lies
    t = opts.t_start
    istart = 0
    while istart < len(tnList)-2 and tnList[istart+1] <= t:
        istart+=1
    assert tnList[istart] <= t and t < tnList[istart+1]
    #this will get cycled through and set as value at tnList[istart] in loop
    #evaluate velocity
    #\todo move out of python
    try:
        p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
    except TypeError:
        for eN in range(mesh.nElements_global):
            for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],tnList[istart])
    #evaluate velocity degrees of freedom from its interpolation conditions
    velocity_new.projectFromInterpolationConditions(q['velocity_interpolation_values'])

    #get points actually going to track
    q['J']            = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['inverse(J)']   = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['det(J)']       = numpy.zeros((q['x'].shape[0],q['x'].shape[1]),'d')
    trialSpace.elementMaps.getJacobianValues(q['x_hat'],q['J'],q['inverse(J)'],q['det(J)'])

    elementBoundaryOuterNormalsArray = numpy.zeros((sdmesh.nElements_global,sdmesh.nElementBoundaries_element,p.nd),'d')
    boundaryNormals = numpy.array(trialSpace.elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
    ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                            q['inverse(J)'],
                                            elementBoundaryOuterNormalsArray)
    q['x_track']  = p.getInitialTrackingLocations(q)#q['x'];
    q['x_depart'] = numpy.copy(q['x_track'])

    setupArbitraryTrackingDataArrays(opts.t_start,tnList[istart+1],sdmesh,elementBoundaryOuterNormalsArray,q)

    if writePT123files:
        q['element_depart_start'] = numpy.copy(q['element_track'])
        q['x_depart_start'] = numpy.copy(q['x_depart'])

    if tryPT123A:

        class Transport_dummy:
            def __init__(self,elementBoundaryOuterNormalsArray):
                self.elementBoundaryOuterNormalsArray = elementBoundaryOuterNormalsArray
        #
        particle_tracker = LinearAdvection_C0P1Velocity_PT123A(sdmesh,p.nd,
                                                               {0:velocitySpace.dofMap.l2g},
                                                               {0:velocity.dof},
                                                               activeComponentList=[0])

        particle_tracker.updateTransportInformation(Transport_dummy(elementBoundaryOuterNormalsArray))
    else:
        particle_tracker = LinearAdvection_C0P1Velocity_PT123(sdmesh,p.nd,
                                                              {0:velocitySpace.dofMap.l2g},
                                                              {0:velocity.dof},
                                                              activeComponentList=[0])
    particle_tracker.setFromOptions(opts)


    ########## visualization and output stuff
    #save velocity field for visualization
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    velocity_new.getValues(q[('u_shape_functions')],q['velocity_depart'])
    q['velocity_track'] = q['velocity_depart'][:]

    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","Lu1D",useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_depart'],tnList[istart],init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    #writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['t_viz'] = numpy.copy(q['t_depart'])

    #which way to track (1. forward, -1. backward)
    direction = 1.0

    #print out point we want to track
    print """t= %s x_depart[0][0]= %s """ % (t,q['x_depart'][0][0])
    #time value for 'old' velocity
    tvelOld = tnList[istart]

    for i,tout in enumerate(tnList[istart+1:]):
        dtout = tout-t
        #update for next step
        velocity.dof[:]=velocity_new.dof[:]

        particle_tracker.setTrackingVelocity(velocity.dof,0,tvelOld,timeLevel=0)
        #time tracking to
        t += dtout

        #evaluate velocity
        #\todo move out of python
        try:
            p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
        except TypeError:
            for eN in range(mesh.nElements_global):
                for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                    q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],t)
        #evaluate velocity degrees of freedom from its interpolation conditions
        velocity_new.projectFromInterpolationConditions(q['velocity_interpolation_values'])
        particle_tracker.setTrackingVelocity(velocity_new.dof,0,t,timeLevel=1)

        #only update times for points that are going to be tracked
        q['t_track'][q['flag_track'] >= -1] = t

        if direction > 0.0:
            particle_tracker.forwardTrack({0:q['t_depart']},
                                          {0:q['t_track']},                          #target end time
                                          {0:q['x_depart'].shape[0]},                #total number of points
                                          {0:q['x_depart']},                         #departure points
                                          {0:q['element_track']},                    #in/out element locations
                                          {0:q['x_track']},                          #arrival points
                                          {0:q['flag_track']})
        else:
            particle_tracker.backwardTrack({0:q['t_depart']},
                                           {0:q['t_track']},                          #target end time
                                           {0:q['x_depart'].shape[0]},                #total number of points
                                           {0:q['x_depart']},                         #departure points
                                           {0:q['element_track']},                    #in/out element locations
                                           {0:q['x_track']},                          #arrival points
                                           {0:q['flag_track']})

        velocity_new.getValues(q[('u_shape_functions')],q['velocity_track'])
        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        #writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)


        q['t_depart'].flat[:] = q['t_track'].flat
        q['x_depart'].flat[:] = q['x_track'].flat
        tvelOld = t

        print """t= %s t_track[0]= %s x_track[0][0] = %s  """ % (t,q['t_track'][0],q['x_track'][0][0])

    #output step loop

    #
    print "Lu Done, t= %s t_track[0]= %s [s],  %s [d], x_track[0][0] = %s  """ % (t,q['t_track'][0],q['t_track'][0]/(60.0*60.0*24.0),q['x_track'][0][0])

    #close out archives and visualization
    archive.close()

    #try to write out PT123 files
    if writePT123files:
        writePT123inputMesh(opts,p,sdmesh,filebase="testLu1D_pt123")
        writePT123nodalVelocity(opts,p,sdmesh,tnList,q,velocitySpace,velocity_new,
                                filebase="testLu1D_pt123")
        writePT123nodalVolumeFraction(opts,p,mesh,particleTracking_params,tnList,q,
                                      filebase="testLu1D_pt123")
        writePT123particleFile(opts,p,opts.t_start,sdmesh,particleTracking_params,tnList,q,tryPT123A,filebase="testLu1D_pt123")


    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass

def test8(opts,tryPT123A=False,writePT123files=True):
    """
    transient, C0P1 in space. Computes error by comparing starting and ending location for points in a specified
      interval
    """
    p = __import__(opts.problem[:-3])

    assert p.nd in [1,2,3]
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    #need separate nnq for velocity eval and tracking points now
    nPointsToTrack = opts.nnq
    opts.nnq = 1
    assert opts.nnq > 0, "nnq = %s not ok must be > 1" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    dtout = (opts.t_target - opts.t_start)/opts.nnt
    tnList = [opts.t_start + i*dtout for i in range(int(opts.nnt)+1)]
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    if p.nd == 3:
        mesh = setupMesh_3d(opts,p)
    elif p.nd == 2:
        mesh = setupMesh_2d(opts,p)
    else:
        mesh = setupMesh_1d(opts,p)

    sdmesh = mesh.subdomainMesh
    #options for tracking
    particleTracking_params = {}
    particleTracking_params[('dn_safe_tracking',0)]=1.0e-7
    particleTracking_params[('temporalVariationFlag',0)]=2
    particleTracking_params[('atol_tracking',0)]=1.0e-7
    particleTracking_params[('rtol_tracking',0)]=0.0
    particleTracking_params[('sf_tracking',0)]=0.9
    if not tryPT123A or writePT123files:
        particleTracking_params[('dt_init',0)]= 0.001
        particleTracking_params[('dt_init_flag',0)] = 1 #how to pick initial time step (0 --> const, 1 --> CFL=1)
        particleTracking_params[('rk_flag',0)] = 45 #RK type

    opts.particleTracking_params = particleTracking_params
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,tnList[0])

    #velocity representation
    #
    velocitySpace,velocity,q['velocity_interpolation_values'] = setupVelocity_C0P1(opts,p,sdmesh,tnList[0])

    component_velocity_times={0:tnList[0]}

    velocity_new = FemTools.FiniteElementFunction(velocitySpace,dim_dof=p.nd,isVector=True,name="velocity_new")
    velocity_new.dof[:] = velocity.dof
    component_velocity_times_new={0:tnList[0]}

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    #find the interval where start time lies
    t = opts.t_start
    istart = 0
    while istart < len(tnList)-2 and tnList[istart+1] <= t:
        istart+=1
    assert tnList[istart] <= t and t < tnList[istart+1]
    #this will get cycled through and set as value at tnList[istart] in loop
    #evaluate velocity
    #\todo move out of python
    try:
        p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
    except TypeError:
        for eN in range(mesh.nElements_global):
            for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],tnList[istart])
    #evaluate velocity degrees of freedom from its interpolation conditions
    velocity_new.projectFromInterpolationConditions(q['velocity_interpolation_values'])

    #get points actually going to track
    q['J']            = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['inverse(J)']   = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['det(J)']       = numpy.zeros((q['x'].shape[0],q['x'].shape[1]),'d')
    trialSpace.elementMaps.getJacobianValues(q['x_hat'],q['J'],q['inverse(J)'],q['det(J)'])

    elementBoundaryOuterNormalsArray = numpy.zeros((sdmesh.nElements_global,sdmesh.nElementBoundaries_element,p.nd),'d')
    boundaryNormals = numpy.array(trialSpace.elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
    ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                            q['inverse(J)'],
                                            elementBoundaryOuterNormalsArray)
    q['nPointsToTrack'] = nPointsToTrack
    q['x_track']  = p.getInitialTrackingLocations(q)#q['x'];
    q['x_depart'] = numpy.copy(q['x_track'])

    #save initial elements for output pt123 input files

    setupArbitraryTrackingDataArrays(opts.t_start,tnList[istart+1],sdmesh,elementBoundaryOuterNormalsArray,q)

    if writePT123files:
        q['element_depart_start'] = numpy.copy(q['element_track'])
        q['x_depart_start'] = numpy.copy(q['x_depart'])

    label = "PT123"
    if tryPT123A:
        label = "PT123A"

        class Transport_dummy:
            def __init__(self,elementBoundaryOuterNormalsArray):
                self.elementBoundaryOuterNormalsArray = elementBoundaryOuterNormalsArray
        #
        particle_tracker = LinearAdvection_C0P1Velocity_PT123A(sdmesh,p.nd,
                                                               {0:velocitySpace.dofMap.l2g},
                                                               {0:velocity.dof},
                                                               activeComponentList=[0])

        particle_tracker.updateTransportInformation(Transport_dummy(elementBoundaryOuterNormalsArray))
    else:
        particle_tracker = LinearAdvection_C0P1Velocity_PT123(sdmesh,p.nd,
                                                              {0:velocitySpace.dofMap.l2g},
                                                              {0:velocity.dof},
                                                              activeComponentList=[0])
    particle_tracker.setFromOptions(opts)


    ########## visualization and output stuff
    #save velocity field for visualization
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    velocity_new.getValues(q[('u_shape_functions')],q['velocity_depart'])
    q['velocity_track'] = q['velocity_depart'][:]

    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","test8_"+label,useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_depart'],tnList[istart],init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    #writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['t_viz'] = numpy.copy(q['t_depart'])

    #which way to track (1. forward, -1. backward)
    direction = 1.0

    xfile = open('x_track_test8_'+label+'.dat','w')
    tfile = open('t_track_test8_'+label+'.dat','w')
    #print out point we want to track
    #print """t= %s x_depart= %s """ % (t,q['x_depart'])
    numpy.savetxt(tfile,q['t_track'])
    numpy.savetxt(xfile,q['x_depart'])
    #time value for 'old' velocity
    tvelOld = tnList[istart]

    for i,tout in enumerate(tnList[istart+1:]):
        dtout = tout-t
        #update for next step
        velocity.dof[:]=velocity_new.dof[:]

        particle_tracker.setTrackingVelocity(velocity.dof,0,tvelOld,timeLevel=0)
        #time tracking to
        t += dtout

        #evaluate velocity
        #\todo move out of python
        try:
            p.analyticalSolutionParticleVelocity[0].uOfXTv(velocitySpace.interpolationPoints,t,q['velocity_interpolation_values'])
        except TypeError:
            for eN in range(mesh.nElements_global):
                for k in range(velocitySpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                    q['velocity_interpolation_values'][eN,k,:] = p.analyticalSolutionParticleVelocity[0].uOfXT(velocitySpace.interpolationPoints[eN,k],t)
        #evaluate velocity degrees of freedom from its interpolation conditions
        velocity_new.projectFromInterpolationConditions(q['velocity_interpolation_values'])
        particle_tracker.setTrackingVelocity(velocity_new.dof,0,t,timeLevel=1)

        #only update times for points that are going to be tracked
        q['t_track'][q['flag_track'] >= -1] = t

        if direction > 0.0:
            particle_tracker.forwardTrack({0:q['t_depart']},
                                          {0:q['t_track']},                          #target end time
                                          {0:q['x_depart'].shape[0]},                #total number of points
                                          {0:q['x_depart']},                         #departure points
                                          {0:q['element_track']},                    #in/out element locations
                                          {0:q['x_track']},                          #arrival points
                                          {0:q['flag_track']})
        else:
            particle_tracker.backwardTrack({0:q['t_depart']},
                                           {0:q['t_track']},                          #target end time
                                           {0:q['x_depart'].shape[0]},                #total number of points
                                           {0:q['x_depart']},                         #departure points
                                           {0:q['element_track']},                    #in/out element locations
                                           {0:q['x_track']},                          #arrival points
                                           {0:q['flag_track']})

        velocity_new.getValues(q[('u_shape_functions')],q['velocity_track'])
        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        #writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)


        numpy.savetxt(tfile,q['t_track'])
        numpy.savetxt(xfile,q['x_track'])

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x_depart'].flat[:] = q['x_track'].flat
        tvelOld = t

        #print """t= %s t_track= %s x_track = %s  """ % (t,q['t_track'],q['x_track'])

    #output step loop

    xfile.close()
    tfile.close()
    #
    #print "Done, t= %s t_track= %s [s],  x_track = %s  """ % (t,q['t_track'],q['x_track'])

    #close out archives and visualization
    archive.close()

    #try to write out PT123 files
    if writePT123files:
        writePT123inputMesh(opts,p,sdmesh,filebase="test8_pt123")
        writePT123nodalVelocity(opts,p,sdmesh,tnList,q,velocitySpace,velocity_new,
                                filebase="test8_pt123")
        writePT123nodalVolumeFraction(opts,p,mesh,particleTracking_params,tnList,q,
                                      filebase="test8_pt123")
        writePT123particleFile(opts,p,opts.t_start,sdmesh,particleTracking_params,tnList,q,tryPT123A,filebase="test8_pt123")

    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass

def test9(opts,tryPT123A=True,writePT123files=True):
    """
    transient, RT0 in space. Computes error by comparing starting and ending location for points in a specified
      interval
    """
    p = __import__(opts.problem[:-3])

    assert p.nd in [1,2,3]
    assert "analyticalSolutionParticleVelocity" in dir(p), "need p.analyticalSolutionParticleVelocity to test tracking"
    #need separate nnq for velocity eval and tracking points now
    nPointsToTrack = opts.nnq
    opts.nnq = 1
    assert opts.nnq > 0, "nnq = %s not ok must be > 1" % opts.nnq
    if opts.t_target == None:
        opts.t_target = p.T
    dtout = (opts.t_target - opts.t_start)/opts.nnt
    tnList = [opts.t_start + i*dtout for i in range(int(opts.nnt)+1)]
    #to collect integration points etc
    q = {}
    #spatial mesh on 1 level for now,
    if p.nd == 3:
        mesh = setupMesh_3d(opts,p)
    elif p.nd == 2:
        mesh = setupMesh_2d(opts,p)
    else:
        mesh = setupMesh_1d(opts,p)

    sdmesh = mesh.subdomainMesh

    #options for tracking
    particleTracking_params = {}
    particleTracking_params[('dn_safe_tracking',0)]=1.0e-6
    particleTracking_params[('temporalVariationFlag',0)]=2
    particleTracking_params[('atol_tracking',0)]=1.0e-6
    particleTracking_params[('rtol_tracking',0)]=0.0
    particleTracking_params[('sf_tracking',0)]=0.9
    if not tryPT123A or writePT123files:
        particleTracking_params[('dt_init',0)]= 0.001
        particleTracking_params[('dt_init_flag',0)] = 1 #how to pick initial time step (0 --> const, 1 --> CFL=1)
        particleTracking_params[('rk_flag',0)] = 45 #RK type

    opts.particleTracking_params = particleTracking_params
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,tnList[0])

    #velocity representation
    #since don't have a full FEM representation for vector spaces this requires more work
    q['velocity_dof'],q['velocity_l2g'],q['velocity_interpolation_values'],ebq,q['elementBoundaryQuadraturePoints'],q['elementBoundaryQuadratureWeights'] = setupVelocity_RT0(opts,p,sdmesh,trialSpace,opts.t_start,nd=p.nd)
    #uses flux rep
    localVelocityRepresentationFlag = 2
    q['velocity_new_dof'] = numpy.copy(q['velocity_dof'])

    component_velocity_times={0:tnList[0]}

    component_velocity_times_new={0:tnList[0]}

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    #find the interval where start time lies
    t = opts.t_start
    istart = 0
    while istart < len(tnList)-2 and tnList[istart+1] <= t:
        istart+=1
    assert tnList[istart] <= t and t < tnList[istart+1]

    #get points actually going to track
    q['J']            = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['inverse(J)']   = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['det(J)']       = numpy.zeros((q['x'].shape[0],q['x'].shape[1]),'d')
    trialSpace.elementMaps.getJacobianValues(q['x_hat'],q['J'],q['inverse(J)'],q['det(J)'])

    elementBoundaryOuterNormalsArray = numpy.zeros((sdmesh.nElements_global,sdmesh.nElementBoundaries_element,p.nd),'d')
    boundaryNormals = numpy.array(trialSpace.elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
    ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                            q['inverse(J)'],
                                            elementBoundaryOuterNormalsArray)
    q['nPointsToTrack'] = nPointsToTrack
    q['x_track']  = p.getInitialTrackingLocations(q)#q['x'];
    q['x_depart'] = numpy.copy(q['x_track'])

    setupArbitraryTrackingDataArrays(opts.t_start,tnList[istart+1],sdmesh,elementBoundaryOuterNormalsArray,q)

    if writePT123files:
        q['abs(det(J))']  = numpy.absolute(q['det(J)'])
        q['element_depart_start'] = numpy.copy(q['element_track'])
        q['x_depart_start'] = numpy.copy(q['x_depart'])

    label = "PT123"
    if tryPT123A:
        label = "PT123A"

        class Transport_dummy:
            def __init__(self,elementBoundaryOuterNormalsArray):
                self.elementBoundaryOuterNormalsArray = elementBoundaryOuterNormalsArray
        #
        particle_tracker = LinearAdvection_RT0Velocity_PT123A(sdmesh,p.nd,
                                                              {0:q['velocity_l2g']},
                                                              {0:q['velocity_dof']},
                                                              activeComponentList=[0])

        particle_tracker.updateTransportInformation(Transport_dummy(elementBoundaryOuterNormalsArray))
    else:
        particle_tracker = LinearAdvection_RT0Velocity_PT123(sdmesh,p.nd,
                                                             {0:q['velocity_l2g']},
                                                             {0:q['velocity_dof']},
                                                             activeComponentList=[0])
    #setting tolerances in opts here is not very convenient because doesn't have numerics file loaded
    particle_tracker.setFromOptions(opts)

    ########## visualization and output stuff

    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","test9_"+label,useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_depart'],tnList[istart],init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['t_viz'] = numpy.copy(q['t_depart'])

    #which way to track (1. forward, -1. backward)
    direction = 1.0

    xfile = open('x_track_test9_'+label+'.dat','w')
    tfile = open('t_track_test9_'+label+'.dat','w')
    #print out point we want to track
    #print """t= %s x_depart= %s """ % (t,q['x_depart'])
    numpy.savetxt(tfile,q['t_track'])
    numpy.savetxt(xfile,q['x_depart'])

    #time value for 'old' velocity
    tvelOld = tnList[istart]

    for i,tout in enumerate(tnList[istart+1:]):
        dtout = tout-t
        #update for next step
        q['velocity_dof'][:]=q['velocity_new_dof']

        particle_tracker.setTrackingVelocity(q['velocity_dof'],0,tvelOld,timeLevel=0)
        #time tracking to
        t += dtout
        evaluateVelocity_RT0(opts,p,sdmesh,t,ebq,q['velocity_new_dof'],q['velocity_l2g'],q['velocity_interpolation_values'])

        particle_tracker.setTrackingVelocity(q['velocity_new_dof'],0,t,timeLevel=1)

        #only update times for points that are going to be tracked
        q['t_track'][q['flag_track'] >= -1] = t

        if direction > 0.0:
            particle_tracker.forwardTrack({0:q['t_depart']},
                                          {0:q['t_track']},                          #target end time
                                          {0:q['x_depart'].shape[0]},                #total number of points
                                          {0:q['x_depart']},                         #departure points
                                          {0:q['element_track']},                    #in/out element locations
                                          {0:q['x_track']},                          #arrival points
                                          {0:q['flag_track']})
        else:
            particle_tracker.backwardTrack({0:q['t_depart']},
                                           {0:q['t_track']},                          #target end time
                                           {0:q['x_depart'].shape[0]},                #total number of points
                                           {0:q['x_depart']},                         #departure points
                                           {0:q['element_track']},                    #in/out element locations
                                           {0:q['x_track']},                          #arrival points
                                           {0:q['flag_track']})

        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)

        numpy.savetxt(tfile,q['t_track'])
        numpy.savetxt(xfile,q['x_track'])

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x_depart'].flat[:] = q['x_track'].flat
        tvelOld = t

        #print """t= %s t_track= %s x_track = %s  """ % (t,q['t_track'],q['x_track'])

    #output step loop

    #
    #print "Done, t= %s t_track= %s [s],  x_track = %s  """ % (t,q['t_track'],q['x_track'])
    xfile.close()
    tfile.close()
    #close out archives and visualization
    archive.close()

    #try to write out PT123 files
    if writePT123files:
        #mwf debug
        #import pdb
        #pdb.set_trace()
        writePT123inputMesh(opts,p,sdmesh,filebase="test9_pt123")
        writePT123RT0velocityAsElementVelocity(opts,p,sdmesh,particleTracking_params,
                                               tnList,q,ebq,
                                               filebase="test9_pt123")
        writePT123elementVolumeFraction(opts,p,mesh,particleTracking_params,tnList,q,
                                        filebase="test9_pt123")
        writePT123particleFile(opts,p,opts.t_start,sdmesh,particleTracking_params,tnList,q,tryPT123A,filebase="test9_pt123")


    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass

def Pearce_Ex7_non_uniform(opts,tryPT123A=True,writePT123files=False):
    """
    1D example from Pearce with transient velocity

    Domain is [-100,100], domain discretized with 21 nodes

    """
    from math import pi,sqrt,exp
    #dummy p module
    class Problem:
        class IC:
            def __init__(self,val=0.0):
                self.val=val
            def uOfXT(self,x,t):
                return self.val
        def __init__(self,filebase="test_nu"):
            #time intervals for velocity evaluation
            #self.tnList = [0.0,200.0,400.0,600.0,800.0,1000.0,1200.0]
            self.L=(200.0,1.0,1.0)
            self.x0=numpy.array([-100.0,0.0,0.0])
            self.filebase=filebase
            vnfile = open(filebase+".vn1",'r')
            line = vnfile.readline()
            self.nn,self.nd,self.nt=map(int,line.split())
            self.tnList=[]
            self.velocity_dofs=[]
            for it in range(self.nt):
                line = vnfile.readline()
                assert 'TS' in line, "line= %s" % line
                self.tnList.append(float(line.split()[1]))
                self.velocity_dofs.append(numpy.zeros((self.nn,),'d'))
                for jn in range(self.nn):
                    line=vnfile.readline()
                    self.velocity_dofs[-1][jn]=float(line)
            line=vnfile.readline()
            assert 'ENDR' in line, " line= %s " % line
            #create mesh on unit interval then transform output as pt123 does
            self.mlMesh = MeshTools.MultilevelEdgeMesh(self.nn,1,1,self.L[0],refinementLevels=1)
            self.mesh = self.mlMesh.meshList[-1]
            self.initialConditions = {0:Problem.IC()}

        def getInitialTrackingLocations(self,q):
            """
            return points that need to be tracked
            """
            nPoints = self.nn
            x = numpy.copy(self.mesh.nodeArray)
            return x
    p = Problem()

    #mwf debug
    #import pdb
    #pdb.set_trace()

    tnList = p.tnList

    #to collect integration points etc
    q = {}
    mesh = p.mesh

    sdmesh = mesh.subdomainMesh
    #solution representation, not used for linear problems

    trialSpace,u,q['u_interpolation_values'] = setupSolution_C0P1(opts,p,sdmesh,tnList[0])

    #velocity representation
    #
    velocitySpace= FemTools.C0_AffineLinearOnSimplexWithNodalBasis(sdmesh,p.nd)
    velocity = FemTools.FiniteElementFunction(velocitySpace,dim_dof=p.nd,isVector=True,name="velocity")
    velocity.dof[:] = p.velocity_dofs[0]

    component_velocity_times={0:tnList[0]}

    velocity_new = FemTools.FiniteElementFunction(velocitySpace,dim_dof=p.nd,isVector=True,name="velocity_new")
    velocity_new.dof[:] = velocity.dof
    component_velocity_times_new={0:tnList[0]}

    q['x_hat'] = setupTrackingPointsInReferenceSpace(opts,p,useGaussianQuadrature=False)
    nq_per_element = q['x_hat'].shape[0]

    #map points to physical space
    q['x'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,3),'d')
    trialSpace.elementMaps.getValues(q['x_hat'],q['x'])

    #go ahead and get shape functions at initial points to evaluate solution, may not want to use
    #this approach later for nonlinear problems (put shape eval in tracking step if need soln at more points)
    q['u_shape_functions'] = numpy.zeros((sdmesh.nElements_global,nq_per_element,trialSpace.max_nDOF_element),'d')
    trialSpace.getBasisValues(q['x_hat'],q['u_shape_functions'])

    #options for tracking
    particleTracking_params = {}
    particleTracking_params[('dn_safe_tracking',0)]=1.0e-7
    particleTracking_params[('temporalVariationFlag',0)]=2
    particleTracking_params[('atol_tracking',0)]=1.0e-7
    particleTracking_params[('rtol_tracking',0)]=0.0
    particleTracking_params[('sf_tracking',0)]=0.9
    if not tryPT123A or writePT123files:
        particleTracking_params[('dt_init',0)]= 0.001
        particleTracking_params[('dt_init_flag',0)] = 1 #how to pick initial time step (0 --> const, 1 --> CFL=1)
        particleTracking_params[('rk_flag',0)] = 45 #RK type

    opts.particleTracking_params = particleTracking_params

    #find the interval where start time lies
    t = opts.t_start
    istart = 0
    while istart < len(tnList)-2 and tnList[istart+1] <= t:
        istart+=1
    assert tnList[istart] <= t and t < tnList[istart+1]

    #this will get cycled through and set as value at tnList[istart] in loop
    velocity_new.dof[:]=p.velocity_dofs[istart]

    #get points actually going to track
    q['J']            = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['inverse(J)']   = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd,p.nd),'d')
    q['det(J)']       = numpy.zeros((q['x'].shape[0],q['x'].shape[1]),'d')
    trialSpace.elementMaps.getJacobianValues(q['x_hat'],q['J'],q['inverse(J)'],q['det(J)'])

    elementBoundaryOuterNormalsArray = numpy.zeros((sdmesh.nElements_global,sdmesh.nElementBoundaries_element,p.nd),'d')
    boundaryNormals = numpy.array(trialSpace.elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
    ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                            q['inverse(J)'],
                                            elementBoundaryOuterNormalsArray)
    q['x_track']  = p.getInitialTrackingLocations(q)#q['x'];
    q['x_depart'] = numpy.copy(q['x_track'])

    setupArbitraryTrackingDataArrays(opts.t_start,tnList[istart+1],sdmesh,elementBoundaryOuterNormalsArray,q)

    if writePT123files:
        q['element_depart_start'] = numpy.copy(q['element_track'])
        q['x_depart_start'] = numpy.copy(q['x_depart'])

    if tryPT123A:

        class Transport_dummy:
            def __init__(self,elementBoundaryOuterNormalsArray):
                self.elementBoundaryOuterNormalsArray = elementBoundaryOuterNormalsArray
        #
        particle_tracker = LinearAdvection_C0P1Velocity_PT123A(sdmesh,p.nd,
                                                               {0:velocitySpace.dofMap.l2g},
                                                               {0:velocity.dof},
                                                               activeComponentList=[0])

        particle_tracker.updateTransportInformation(Transport_dummy(elementBoundaryOuterNormalsArray))
    else:
        particle_tracker = LinearAdvection_C0P1Velocity_PT123(sdmesh,p.nd,
                                                              {0:velocitySpace.dofMap.l2g},
                                                              {0:velocity.dof},
                                                              activeComponentList=[0])
    particle_tracker.setFromOptions(opts)


    ########## visualization and output stuff
    #save velocity field for visualization
    q['velocity_depart'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')
    q['velocity_track'] = numpy.zeros((q['x'].shape[0],q['x'].shape[1],p.nd),'d')

    velocity_new.getValues(q[('u_shape_functions')],q['velocity_depart'])
    q['velocity_track'] = q['velocity_depart'][:]

    #mwf debug
    #pdb.set_trace()
    archive = Archiver.XdmfArchive(".","testEx7_1d_non-uniform",useTextArchive=False)
    import xml.etree.ElementTree as ElementTree
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    writer   = Archiver.XdmfWriter()

    #
    arGrid = writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_depart'],tnList[istart],init=True,tCount=0)
    writer.writeScalarXdmf_particles(archive,q['t_depart'],"t_particle",tCount=0)
    #writer.writeVectorXdmf_particles(archive,q['velocity_depart'],"velocity",tCount=0)

    #to accomodate updating pardigm for visualization copy _track info to be consistent
    #with start locations
    q['t_viz'] = numpy.copy(q['t_depart'])

    #which way to track (1. forward, -1. backward)
    direction = 1.0

    #print out point we want to track
    print """t= %s x_depart[0][0]= %s """ % (t,q['x_depart'][0][0])
    #time value for 'old' velocity
    tvelOld = tnList[istart]

    for i,tout in enumerate(tnList[istart+1:]):
        lastStep = False
        if tout > opts.t_target:
            lastStep = True
        dtout = tout-t
        #update for next step
        velocity.dof[:]=velocity_new.dof[:]

        particle_tracker.setTrackingVelocity(velocity.dof,0,tvelOld,timeLevel=0)
        #time tracking to
        t += dtout
        #evaluate velocity
        velocity_new.dof[:]=p.velocity_dofs[istart+1+i]
        particle_tracker.setTrackingVelocity(velocity_new.dof,0,t,timeLevel=1)
        #mwf debug
        print "Pearce non-uniform 1d t= %s velocity_new= %s " %(t,velocity_new.dof)
        #only update times for points that are going to be tracked
        q['t_track'][q['flag_track'] >= -1] = min(t,opts.t_target)

        if direction > 0.0:
            particle_tracker.forwardTrack({0:q['t_depart']},
                                          {0:q['t_track']},                          #target end time
                                          {0:q['x_depart'].shape[0]},                #total number of points
                                          {0:q['x_depart']},                         #departure points
                                          {0:q['element_track']},                    #in/out element locations
                                          {0:q['x_track']},                          #arrival points
                                          {0:q['flag_track']})
        else:
            particle_tracker.backwardTrack({0:q['t_depart']},
                                           {0:q['t_track']},                          #target end time
                                           {0:q['x_depart'].shape[0]},                #total number of points
                                           {0:q['x_depart']},                         #departure points
                                           {0:q['element_track']},                    #in/out element locations
                                           {0:q['x_track']},                          #arrival points
                                           {0:q['flag_track']})

        velocity_new.getValues(q[('u_shape_functions')],q['velocity_track'])
        #mwf scale by domain offset for output
        q['x_track'][:,0] += p.x0[0]
        writer.writeMeshXdmf_particles(archive,sdmesh,p.nd,q['x_track'],t,init=False,meshChanged=True,arGrid=arGrid,tCount=i+1)
        writer.writeScalarXdmf_particles(archive,q['t_track'],"t_particle",tCount=i+1)
        #writer.writeVectorXdmf_particles(archive,q['velocity_track'],"velocity",tCount=i+1)
        #convert back to normalized domain
        q['x_track'][:,0] -= p.x0[0]

        q['t_depart'].flat[:] = q['t_track'].flat
        q['x_depart'].flat[:] = q['x_track'].flat
        tvelOld = t

        print """t= %s t_track[0]= %s x_track[0][0] = %s  """ % (t,q['t_track'][0],q['x_track'][0][0])
        if lastStep: break
    #output step loop

    #
    print "Pearce's Ex7 Done, t= %s t_track[0]= %s,  x_track[0][0] = %s  """ % (t,q['t_track'][0],q['x_track'][0][0])

    #close out archives and visualization
    archive.close()

    #try to write out PT123 files
    if writePT123files:
        writePT123inputMesh(opts,p,sdmesh,filebase="testEx7_1d_non-uniform_pt123")
        writePT123nodalVelocity(opts,p,sdmesh,tnList,q,velocitySpace,velocity_new,
                                filebase="testEx7_1d_non-uniform_pt123")
        writePT123nodalVolumeFraction(opts,p,mesh,particleTracking_params,tnList,q,
                                      filebase="testEx7_1d_non-uniform_pt123")
        writePT123particleFile(opts,p,opts.t_start,sdmesh,particleTracking_params,tnList,q,tryPT123A,filebase="testEx7_1d_non-uniform_pt123")


    try:#run time visualization?
        sys.exit(vtkViewers.g.app.exec_())
    except:
        pass

if __name__=='__main__':
    """
    TODO add verbose flag
         how to visualize?
         add error calculation
    """
    import optparse
    import sys
    import pdb
    try:
        import cProfile as profiler
    except:
        import profile as profiler
    import pstats

    usage  = "usage: %prog [options] "
    parser = optparse.OptionParser(usage=usage)

    parser.add_option('--probDir',
                      default='.',
                      action="store",
                      help="""where to find problem descriptions""")
    parser.add_option("-p", "--problem",
                      help="problem definition for tracking",
                      action="store",
                      type="string",
                      dest="problem",
                      default="la_gauss_ellam_1d_p.py")
    parser.add_option("--nnx",
                      help="number of nodes in the x direction",
                      action="store",
                      type="int",
                      dest="nnx",
                      default=21)
    parser.add_option("--nny",
                      help="number of nodes in the y direction",
                      action="store",
                      type="int",
                      dest="nny",
                      default=21)
    parser.add_option("--nnz",
                      help="number of nodes in the z direction",
                      action="store",
                      type="int",
                      dest="nnz",
                      default=21)
    parser.add_option("--nnq",
                      help="number of points to track per element must be >= 1 for now",
                      action="store",
                      type="int",
                      dest="nnq",
                      default=3)
    parser.add_option("--nnt",
                      help="number of time intervals over which to track",
                      action="store",
                      type="float",
                      dest="nnt",
                      default=1)
    parser.add_option("--t_start",
                      action="store",
                      help="initial, uniform time to start tracking",
                      type="float",
                      dest="t_start",
                      default=0.0)
    parser.add_option("-T","--t_target",
                      action="store",
                      help="initial, uniform time to start tracking",
                      type="float",
                      dest="t_target",
                      default=None)
    parser.add_option("--test_id",
                      help="flag for which test to run",
                      action="store",
                      type="int",
                      dest="test_id",
                      default=0)
    parser.add_option("--view",
                      help="try run time visualization",
                      action="store_true",
                      dest="view",
                      default=False)
    parser.add_option("--wait",
                      help="wait during run time visualization",
                      action="store_true",
                      dest="wait",
                      default=False)
    parser.add_option("--profile",
                      help="profile entire tracking example",
                      action="store_true",
                      dest="profile",
                      default=False)


    (opts,args) = parser.parse_args()

    #pdb.set_trace()

    #initialize mpi
    Comm.argv = sys.argv[:1]
    comm = Comm.init()

    probDir = str(opts.probDir)
    if probDir not in sys.path:
        sys.path.insert(0,probDir)

    testDict = {0:test0,1:test1,2:test2,3:test3,4:test4,5:test5,6:test6,7:Lu1Dex,8:test8,9:test9,
                10:Pearce_Ex7_non_uniform}
    test = None
    assert opts.test_id in testDict.keys(), "test_id= %s not supported yet " % opts.test_id

    if opts.profile:
        name = 'Tracking_test_id_%s' % opts.test_id
        profiler.run('testDict[opts.test_id](opts)',name+'_prof')
        stats = pstats.Stats(name+'_prof')
        stats.strip_dirs()
        stats.dump_stats(name+'_prof_c')
        stats.sort_stats('cumulative')
        stats.print_stats(30)
        stats.sort_stats('time')
        stats.print_stats(30)
    else:
        testDict[opts.test_id](opts)
