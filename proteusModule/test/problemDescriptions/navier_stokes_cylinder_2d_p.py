"""
Incompressible Navier-Stokes flow around a square obstacle in 2D.
"""
from pyadh import *
from pyadh.default_p import *
import sys
from math import *
## \page Tests Test Problems
#\ref navier_stokes_squareobstacle_2d_p.py "Incompressible Navier-Stokes flow around a square obstacle"
# 

##\ingroup test
# \file navier_stokes_squareobstacle_2d_p.py
# @{
#
# \brief Incompressible flow around a square with no-slip boundary conditions on the square.
#
# The governing equations are described by the pyadh::TransportCoefficients::NavierStokes class. The domain and initial/boundary conditions are
#
# \f{eqnarray*}
# \Omega &=& \Omega_E \backslash \Omega_O  \\
# \Omega_E &=& [0,20] \times [0,10]   \\
# \Omega_O &=& [2,3] \times [5,6]   \\
# u&=&0 \mbox{ on } \partial \Omega_O   \\
# v&=&0 \mbox{ on } \partial \Omega_O   \\
# u &=& 100 \mbox{ on } \{0\} \times [0,10]    \\
# v &=& 0 \mbox{ on } [0,20] \times \{0,10\}   \\
# \frac{\partial u}{\partial x} &=& \{0\} \mbox{ on } 20 \times [0,10]   \\
# \frac{\partial v}{\partial x} &=& \{0\} \mbox{ on } 20 \times [0,10]   \\
# \f}
#
# The flow generates vortices behind the obstacle, which break off periodically and are carried down stream (a Karman vortex street).
#
#\image html squareobstacle_p.jpg
#\image html squareobstacle_v.jpg
#
# <A href="https://adh.usace.army.mil/pyadh-images/squareobstacle.avi"> AVI Animation of velocity relative to free stream velocity</A>
#
nd = 2


rho=998.2
nu=1.004e-6
#rho=1.0
#nu=1.0e-3
inflow_height=0.41
bottom_length=2.2
cylinder_radius=0.1/2.0
cylinder_center = (0.15+0.1/2.0,0.15+0.1/2.0)
cylinder_center = (1.1,0.15+0.1/2.0)
Um = 0.01
#Um = 0.3
Um = 1.5
#REs = raw_input("Enter Reynolds Number\n")
REs = "1.0"#"1.0"#"50.0"#"120.0"
name = "cylinder2d_"+REs+"_"+"sgs_"
#name = "cylinder2d_"+REs+"_test"
RE = float(REs)
Ubar = nu*RE/(2.0*cylinder_radius)
Um = 3.0*Ubar/2.0
#T=8.0
#T = 5.0e-2*bottom_length/Um
#Ubar = 2.0*Um/3.0
print "REYNOLDS NUMBER = ",Ubar*2.0*cylinder_radius/nu
#old domain
#import cylinder2d
#polyfile = "cylinder"
#
#boundaryTags = cylinder2d.genPoly(polyfile,
#                                  cross_section=cylinder2d.circular_cross_section,
#                                  height=inflow_height,
#                                  length=bottom_length,
#                                  radius=cylinder_radius,
#                                  center = cylinder_center,
#                                  n_points_on_obstacle=2*41-2)

#new domain
from cylinder2dDomain import *
domain = cylinder2d(cross_section=circular_cross_section,
                    height=inflow_height,
                    length=bottom_length,
                    radius=cylinder_radius,
                    center = cylinder_center,
                    n_points_on_obstacle=2*41-2)
domain.writePoly("cylinder")
boundaryTags = domain.boundaryTags
#
    

initialConditions = None

analyticalSolution = {}
from pyadh import RANS2P
useOpt=False#2D opt stuff is not finished
if useOpt:
    LevelModelType = RANS2P.OneLevelRANS2P
bcsTimeDependent=True
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=rho,nu_0=nu,
                                             rho_1=rho,nu_1=nu,
                                             g=[0.0,-9.8],
                                             nd=2,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)
coefficients.waterLevel = 0.5*inflow_height
mu = coefficients.rho_0*coefficients.nu_0

#now define the Dirichlet boundary conditions

residence_time = bottom_length/Um

# import ode
# movingDomain=True
# class RigidCylinder(AuxiliaryVariables.AV_base):
#     def __init__(self,D=1.0,Ubar=1.0,rho=1.0,center=(0,0,0),radius=1.0):
#         self.mass = math.pi*radius**2 * rho
#         self.C_fact = 2.0/(rho*D*Ubar**2)
#         self.world = ode.World()
#         self.world.setGravity( (0,0,0) )
#         # Create a body inside the world
#         self.body = ode.Body(self.world)
#         self.M = ode.Mass()
#         self.M.setSphere(rho,radius)
#         self.M.mass = 1.0
#         self.body.setMass(self.M)
#         self.body.setPosition(center)
#         self.last_position=center
#         self.position=center
#         self.last_velocity=(0.0,0.0,0.0)
#         self.velocity=(0.0,0.0,0.0)
#     def attachModel(self,model,ar):
#         import copy
#         self.model=model
#         self.ar=ar
#         self.writer = Archiver.XdmfWriter()
#         self.nd = model.levelModelList[-1].nSpace_global
#         m = self.model.levelModelList[-1]
#         flagMax = max(m.mesh.elementBoundaryMaterialTypes)
#         flagMin = min(m.mesh.elementBoundaryMaterialTypes)
#         assert(flagMin == 0)
#         assert(flagMax >= 0)
#         self.nForces=flagMax+1
#         self.levelFlist=[]
#         for m in self.model.levelModelList:
#             if self.nd == 2:
#                 F = numpy.zeros((self.nForces,2),'d')
#             elif self.nd == 3:
#                 F = numpy.zeros((self.nForces,3),'d')
#             else:
#                 logEvent("Can't use stress computation for nd = "+`self.nd`)
#                 F=None
#             self.levelFlist.append(F)
#         self.historyF=[]
#         self.historyF.append(copy.deepcopy(self.levelFlist))
# #         try:
# #             self.dragWindow=Viewers.windowNumber
# #             self.viewForce()
# #         except:
# #             pass
#         return self
#     def get_u(self):
#         #print "velocity0",self.last_velocity[0]
#         return self.last_velocity[0]
#     def get_v(self):
#         #print "velocity1",self.last_velocity[1]
#         return self.last_velocity[1]
#     def calculate(self):
#         import pdb
#         for m,F in zip(self.model.levelModelList,self.levelFlist):
#             F.flat[:]=0.0
#             if self.nd ==3:
#                 cfemIntegrals.calculateExteriorElementBoundaryStress3D(m.mesh.elementBoundaryMaterialTypes,
#                                                                        m.mesh.exteriorElementBoundariesArray,
#                                                                        m.mesh.elementBoundaryElementsArray,
#                                                                        m.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                        m.ebqe[('u',0)],#pressure
#                                                                        m.ebqe[('velocity',1)],#mom_flux_vec_u
#                                                                        m.ebqe[('velocity',2)],#mom_flux_vec_v
#                                                                        m.ebqe[('velocity',3)],#mom_flux_vec_w
#                                                                        m.ebqe[('dS_u',0)],#dS
#                                                                        m.ebqe[('n')],
#                                                                        F)
#             if self.nd == 2:
#                 cfemIntegrals.calculateExteriorElementBoundaryStress2D(m.mesh.elementBoundaryMaterialTypes,
#                                                                        m.mesh.exteriorElementBoundariesArray,
#                                                                        m.mesh.elementBoundaryElementsArray,
#                                                                        m.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                        m.ebqe[('u',0)],#pressure
#                                                                        m.ebqe[('velocity',1)],#mom_flux_vec_u
#                                                                        m.ebqe[('velocity',2)],#mom_flux_vec_v
#                                                                        m.ebqe[('dS_u',0)],#dS
#                                                                        m.ebqe[('n')],
#                                                                        F)
#             #pdb.set_trace()
#             logEvent("Force")
#             logEvent(`F`)
#             Ftot=F[0,:]
#             for ib in range(1,self.nForces):
#                 Ftot+=F[ib,:]
#             logEvent("Total force on all boundaries")
#             logEvent(`Ftot`)
#         logEvent("Drag Force " +`self.model.stepController.t_model`+" "+`F[-1,0]`)
#         logEvent("Lift Force " +`self.model.stepController.t_model`+" "+`F[-1,1]`)
#         self.body.addForce((F[-1,0],F[-1,1],0.0))
#         #assume moving in the x direction
#         #self.body.addForce((0.0,F[-1,1],0.0))
#         self.world.step(self.model.stepController.dt_model)
#         #f = m a = m (v_new - v_old)/dt
#         #f dt/m = v_new - v_old
#         #v_new = v_old + f dt/m 
#         self.velocity = (self.last_velocity[0]+F[-1,0]*self.model.stepController.dt_model/self.mass,
#                          self.last_velocity[1]+F[-1,1]*self.model.stepController.dt_model/self.mass,
#                          0.0)
#         self.position = (self.last_position[0]+self.velocity[0]*self.model.stepController.dt_model,
#                          self.last_position[1]+self.velocity[1]*self.model.stepController.dt_model,
#                          0)
#         self.h = (self.velocity[0]*self.model.stepController.dt_model,
#                   self.velocity[1]*self.model.stepController.dt_model,
#                   0)
#         x,y,z = self.body.getPosition()
#         u,v,w = self.body.getLinearVel()
#         self.position=(x,y,z)
#         self.velocity=(u,v,w)
#         self.h = (self.velocity[0]*self.model.stepController.dt_model,
#                   self.velocity[1]*self.model.stepController.dt_model,
#                   self.velocity[2]*self.model.stepController.dt_model)
#         self.last_velocity=self.velocity
#         self.last_position=self.position
#         print "%1.2fsec: pos=(%6.3f, %6.3f, %6.3f)  vel=(%6.3f, %6.3f, %6.3f)" % \
#             (self.model.stepController.t_model, 
#              self.position[0], self.position[1], self.position[2], 
#              self.velocity[0],self.velocity[1],w)
#         self.h = (self.velocity[0]*self.model.stepController.dt_model,
#                   self.velocity[1]*self.model.stepController.dt_model,
#                   0)
#         print "dt model in object ",self.model.stepController.dt_model
#         logEvent("Drag Coefficient " +`self.model.stepController.t_model`+" "+`self.C_fact*F[-1,0]`)
#         logEvent("Lift Coefficient " +`self.model.stepController.t_model`+" "+`self.C_fact*F[-1,1]`)
# #        for ib in range(self.nForces):
# #             self.writeScalarXdmf(self.C_fact*F[ib,0],"Drag Coefficient %i" % (ib,))
# #             self.writeScalarXdmf(self.C_fact*F[ib,1],"Lift Coefficient %i" % (ib,))
# #        self.historyF.append(copy.deepcopy(self.levelFlist))
# #         try:
# #             self.viewForce()
# #         except:
# #             pass
# #     def viewForce(self):
# #         tList=[]
# #         FxList=[]
# #         FyList=[]
# #         for ti,levelFlist in enumerate(self.historyF):
# #             tList.append(ti)
# #             FxList.append(self.C_fact*levelFlist[-1][-1,0])
# #             FyList.append(self.C_fact*levelFlist[-1][-1,1])
# #         print "In view force",FxList
# #         Viewers.windowNumber=self.dragWindow
# #         vtkViewers.viewScalar_1D(numpy.array(tList),numpy.array(FxList),"t","Fx","Drag Coefficients",Viewers.windowNumber,
# #                                  Pause=False,
# #                                  sortPoints=False,newMesh=True)
# #         Viewers.newPlot()
# #         Viewers.newWindow()
# #         vtkViewers.viewScalar_1D(numpy.array(tList),numpy.array(FyList),"t","Fy","Lift Coefficients",Viewers.windowNumber,
# #                                  Pause=False,
# #                                  sortPoints=False,newMesh=True)
# #         Viewers.newPlot()
# #         Viewers.newWindow()

# rc = RigidCylinder(D=2.0*cylinder_radius,Ubar=Ubar,rho=0.5*rho,center=cylinder_center+(0.0,),radius=cylinder_radius)
def velRamp(t):
#    return 1.0
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def getDBC_pressure_tube(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: -(inflow_height-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

noSlip = []#[boundaryTags['bottom'],boundaryTags['top']]
slip = [boundaryTags['bottom'],boundaryTags['top']]
#noSlip = [boundaryTags['obstacle']]#,boundaryTags['bottom'],boundaryTags['top']]
          
def getDBC_u_tube(x,flag):
    if flag == boundaryTags['upstream']:
        #return lambda x,t: velRamp(t)*4.0*Um*x[1]*(inflow_height-x[1])/(inflow_height**2)
        return lambda x,t: velRamp(t)*Um
    elif flag == boundaryTags['obstacle']:
        #return lambda x,t: rc.get_u()
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['obstacle']:
        #return lambda x,t: rc.get_v()
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}


def getAFBC_p_tube(x,flag):
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0
    elif flag in slip:
        return lambda x,t: 0.0

def getDFBC_downstream(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: 0.0
    elif flag in slip:
        return lambda x,t: 0.0        

def getNoBC(x,flag):
    pass

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,1:getNoBC,2:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},1:{1:getDFBC_downstream},2:{2:getDFBC_downstream}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        a = -(inflow_height-x[1])*coefficients.rho*coefficients.g[1]
        return a
class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}

T = 0.1*residence_time#100.0*residence_time

## @}
