from collections import defaultdict

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from . import Comm
from .AuxiliaryVariables import AV_base

class PointGauges(AV_base):
    """ Monitor fields at specific values.

    This class provides a generic point-wise monitor that can be instantiated and attached to Proteus simulations by
    including them in the list of Auxiliary Variables in problem setup.

    Each PointGauges instance may contain one or more fields, which may contain one or more point locations to
    monitor.  The monitoring is defined over a given time and sample rate, and a filename is also supplied.  All
    results are serialized to a CSV file.

    Parallel Implementation Notes:
    After the gauge has been attached, all processes are partitioned into Gauge Owners and non-Gauge Owners.  The
    calculate method is a "no-op" for non-Owners.  For Gauge Owners, all values are computed individually,
    then collectively transmitted to the "root" process, which is the only process responsible for serializing gauge
    results to disk.  This code has not been aggressively vetted for parallel correctness or scalability.

    """
    def __init__(self, gauges,
                 activeTime=None,
                 sampleRate=0,
                 fileName='point_gauges.csv'):
        """
        Create a set of point gauges that will automatically be serialized as CSV data to the requested file.

        :param gauges: An iterable of "gauges".  Each gauge is specified by a 2-tuple, with the first element in the
        tuple a set of fields to be monitored, and the second element the 3-space representation of the gauge
        location.
        :param activeTime: If not None, a 2-tuple of start time and end time for which the point gauge is active.
        :param sampleRate: The intervals at which samples should be measured.  Note that this is a rough lower
        bound, and that the gauge values could be computed less frequently depending on the time integrator.  The
        default value of zero computes the gauge values at every time step.
        :param fileName: The name of the file to serialize results to.

        Example:

        p = PointGauges(self, gauges=((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))),
                                      (('p',), ((0.5, 0.5, 0),))),
                        activeTime=(0, 2.5),
                        sampleRate=0.2,
                        fileName='combined_gauge_0_0.5_sample_all.csv'):

        This creates a PointGauges object that will monitor the u and v fields at the locations [0.5, 0.5,
        0] and [1, 0.5, 0], and the p field at [0.5, 0.5, 0] at simulation time between = 0 and 2.5 with samples
        taken no more frequently than every 0.2 seconds.  Results will be saved to:
        combined_gauge_0_0.5_sample_all.csv.
        """


        AV_base.__init__(self)
        self.gauges = gauges
        self.measuredFields = set()

        # build up dictionary of location information from gauges
        # dictionary of dictionaries, outer dictionary is keyed by location (3-tuple)
        # inner dictionaries contain monitored fields, and closest node
        # closest_node is None if this process does not own the node
        self.locations = {}
        for gauge in gauges:
            fields, locations = gauge
            self.measuredFields.update(fields)
            for location in locations:
                # initialize new dictionary of information at this location
                if location not in self.locations:
                    l_d = {'fields': set()}
                    self.locations[location] = l_d
                # add any currently unmonitored fields
                self.locations[location]['fields'].update(fields)

        self.activeTime = activeTime
        self.sampleRate = sampleRate
        self.fileName = fileName
        self.file = None  # only the root process should have a file open
        self.flags = {}
        self.files = {}
        self.outputWriterReady = False

    def findNearestNode(self, location):
        """Given a gauge location, attempts to locate the most suitable process for monitoring information about
        this location, as well as the node on the process closest to the location.
        """

        # determine local nearest node distance
        node_distances = np.linalg.norm(self.vertices - location, axis=1)
        nearest_node = np.argmin(node_distances)
        nearest_node_distance = node_distances[nearest_node]

        # determine global nearest node
        comm = Comm.get().comm.tompi4py()
        global_min_distance, owning_proc = comm.allreduce(nearest_node_distance, op=MPI.MINLOC)
        if comm.rank != owning_proc:
            nearest_node = None

        assert owning_proc is not None
        return owning_proc, nearest_node

    def buildQuantityRow(self, m, quantity_id, quantity):
        """ Builds up contributions to gauge operator from the underlying element space
        """

        location, node = quantity
        # get the FiniteElementFunction object for this quantity
        femFun = self.model.levelModelList[-1].u[quantity_id]
        # get the mesh for this quantity
        mesh = femFun.femSpace.mesh
        # search elements that contain the nearest node
        # use nearest node if the location is not found on any elements
        for eOffset in range(mesh.nodeElementOffsets[node], mesh.nodeElementOffsets[node + 1]):
            eN = femFun.femSpace.mesh.nodeElementsArray[eOffset]
            # evalute the inverse map for element eN
            xi = femFun.femSpace.elementMaps.getInverseValue(eN, location)
            # query whether xi lies within the reference element
            if femFun.femSpace.elementMaps.referenceElement.onElement(xi):
                for i, psi in enumerate(femFun.femSpace.referenceFiniteElement.localFunctionSpace.basis):
                    # assign quantity weights here
                    m[quantity_id, femFun.femSpace.dofMap.l2g[eN, i]] = psi(xi)
                break
        else:
            # just use nearest node for now if  we're given a point outside the domain.
            # the ideal thing would be to find the element with the nearest face
            m[quantity_id, node] = 1

    def initOutputWriter(self):
        """ Initialize communication strategy for collective output of gauge data.

        On the root process in this communicator, create a map of quantity owners and the corresponding location in
        their arrays.  This process is responsible for collecting gauge data and saving it to disk.

        Gauge data is globally ordered by field, then by location id (as ordered by globalMeasuredQuantities)
        """

        numLocalQuantities = sum([len(self.measuredQuantities[field]) for field in self.measuredFields])
        self.localQuantitiesBuf = np.zeros(numLocalQuantities)

        if self.gaugeComm.rank != 0:
            self.globalQuantitiesBuf = None
        else:
            self.file = open(self.fileName, 'w')

            self.quantityIDs = [0] * self.gaugeComm.size
            for field in self.measuredFields:
                for id in range(len(self.globalMeasuredQuantities[field])):
                    location, owningProc = self.globalMeasuredQuantities[field][id]
                    gaugeProc = self.globalGaugeRanks[owningProc]
                    quantityID = self.quantityIDs[gaugeProc]
                    self.quantityIDs[gaugeProc] += 1
                    assert gaugeProc >= 0
                    self.globalMeasuredQuantities[field][id] = location, gaugeProc, quantityID

            numGlobalQuantities = sum([len(self.globalMeasuredQuantities[field]) for field in self.measuredFields])
            self.globalQuantitiesBuf = np.zeros(numGlobalQuantities)

            # determine mapping from global measured quantities to communication buffer
            self.globalQuantitiesMap = np.zeros(numGlobalQuantities, dtype=np.int)
            i = 0
            for field in self.measuredFields:
                for location, gaugeProc, quantityID in self.globalMeasuredQuantities[field]:
                    self.globalQuantitiesMap[i] = sum(self.quantityIDs[:gaugeProc - 1]) + quantityID
                    assert self.globalQuantitiesMap[i] < numGlobalQuantities
                    i += 1

            # a couple consistency checks
            assert sum(self.quantityIDs) == numGlobalQuantities
            assert all(quantityID > 0 for quantityID in self.quantityIDs)

        self.outputWriterReady = True

    def buildGaugeComm(self):
        """ Create a communicator composed only of processes that own gauge quantities.

        Collective over global communicator.  Builds a local communicator for collecting all gauge data.
        This communicator contains only processes that will contain gauge data.
        """

        comm = Comm.get().comm.tompi4py()

        gaugeOwners = set()

        for field in self.measuredFields:
            for location, owningProc in self.globalMeasuredQuantities[field]:
                gaugeOwners.update((owningProc,))

        self.isGaugeOwner = comm.rank in gaugeOwners
        gaugeComm = comm.Split(color=self.isGaugeOwner)
        if self.isGaugeOwner:
            self.gaugeComm = gaugeComm
            gaugeRank = self.gaugeComm.rank
        else:
            self.gaugeComm = None
            gaugeRank = -1
        self.globalGaugeRanks = comm.allgather(gaugeRank)

    def identifyMeasuredQuantities(self):
        """ build measured quantities, a list of fields

        Each field in turn contains a list of gauge locations and their accompanying nearest node
        only local quantities are saved
        """

        self.measuredQuantities = defaultdict(list)
        self.globalMeasuredQuantities = defaultdict(list)

        for location, l_d in self.locations.iteritems():
            owningProc, nearestNode = self.findNearestNode(location)
            l_d['nearest_node'] = nearestNode
            for field in l_d['fields']:
                self.globalMeasuredQuantities[field].append((location, owningProc))
                if l_d['nearest_node'] is not None:
                    self.measuredQuantities[field].append((location, l_d['nearest_node']))

    def buildGaugeOperators(self):
        """ Build the linear algebra operators needed to compute the gauges.

        The operators are all local since the gauge measurements are calculated locally.
        """

        num_owned_nodes = self.model.levelModelList[-1].mesh.nNodes_global

        self.m = []
        self.field_ids = []
        self.dofsVecs = []
        self.gaugesVecs = []

        for field in self.measuredFields:
            m = PETSc.Mat().create(PETSc.COMM_SELF)
            m.setSizes([len(self.measuredQuantities[field]), num_owned_nodes])
            m.setType('aij')
            m.setUp()
            # matrices are a list in same order as fields
            self.m.append(m)
            field_id = self.fieldNames.index(field)
            self.field_ids.append(field_id)
            # dofs are a list in same order as fields as well
            dofs = self.model.levelModelList[-1].u[field_id].dof
            dofsVec = PETSc.Vec().createWithArray(dofs, comm=PETSc.COMM_SELF)
            self.dofsVecs.append(dofsVec)

        for field, m in zip(self.measuredFields, self.m):
            for quantity_id, quantity in enumerate(self.measuredQuantities[field]):
                self.buildQuantityRow(m, quantity_id, quantity)
            gaugesVec = PETSc.Vec().create(comm=PETSc.COMM_SELF)
            gaugesVec.setSizes(len(self.measuredQuantities[field]))
            gaugesVec.setUp()
            self.gaugesVecs.append(gaugesVec)

        for m in self.m:
            m.assemble()

    def attachModel(self, model, ar):
        """ Attach this gauge to the given simulation model.
        """

        self.model = model
        self.fieldNames = model.levelModelList[-1].coefficients.variableNames
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray

        self.m = {}

        self.identifyMeasuredQuantities()
        self.buildGaugeComm()

        if self.isGaugeOwner:
            self.initOutputWriter()
            self.buildGaugeOperators()
            self.outputHeader()
            # this is currently broken for initial time, need to fix initial model time
            # or enforce that calculate is called as soon as possible
            # after model time is set up
            # time = self.get_time()
            time = 0
            self.outputRow(time)
            self.last_output = time
        return self

    def get_time(self):
        """ Returns the current model time"""
        return self.model.levelModelList[-1].timeIntegration.t

    def outputHeader(self):
        """ Outputs a single header for a CSV style file to self.file"""

        assert self.isGaugeOwner

        if self.gaugeComm.rank == 0:
            self.file.write("%10s" % ('time',))
            for field in self.measuredFields:
                for quantity in self.globalMeasuredQuantities[field]:
                    location, gaugeProc, quantityID = quantity
                    self.file.write(",%s [%9.5g %9.5g %9.5g]" % (field, location[0], location[1], location[2]))
            self.file.write('\n')

    def outputRow(self, time):
        """ Outputs a single row of currently calculated gauge data to self.file"""

        assert self.isGaugeOwner

        self.localQuantitiesBuf = np.concatenate([gaugesVec.getArray() for gaugesVec in self.gaugesVecs])
        self.gaugeComm.Gatherv(self.localQuantitiesBuf, self.globalQuantitiesBuf)

        if self.gaugeComm.rank == 0:
            self.file.write("%10.4e" % time)
            for id in self.globalQuantitiesMap:
                self.file.write(", %36.18e" % (self.globalQuantitiesBuf[id],))
            self.file.write('\n')
            # disable this for better performance, but risk of data loss on crashes
            self.file.flush()
        self.last_output = time

    def calculate(self):
        """ Computes current gauge values, updates open output files
        """

        if not self.isGaugeOwner:
            return

        time = self.get_time()

        if self.activeTime[0] <= time <= self.activeTime[1] and time >= self.last_output + self.sampleRate:
            for field, m, dofsVec, gaugesVec in zip(self.measuredFields, self.m, self.dofsVecs, self.gaugesVecs):
                m.mult(dofsVec, gaugesVec)
            self.outputRow(time)


# this has not been ported to the new-style format
class LineGauges(AV_base):
    def __init__(self, gaugeEndpoints={'pressure_1': ((0.5, 0.5, 0.0), (0.5, 1.8, 0.0))}, linePoints=10):

        AV_base.__init__(self)
        self.endpoints = gaugeEndpoints
        self.flags = {}
        self.linepoints = {}
        self.files = {}  #while open later
        pointFlag = 1000
        for name, (pStart, pEnd) in self.endpoints.iteritems():
            self.flags[name] = pointFlag
            p0 = np.array(pStart)
            direction = np.array(pEnd) - p0
            self.linepoints[name] = []
            for scale in np.linspace(0.0, 1.0, linePoints):
                self.linepoints[name].append(p0 + scale * direction)
            pointFlag += 1

    def attachModel(self, model, ar):
        self.model = model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.tt = model.levelModelList[-1].timeIntegration.t
        self.p = model.levelModelList[-1].u[0].dof
        self.u = model.levelModelList[-1].u[1].dof
        self.v = model.levelModelList[-1].u[2].dof
        return self

    def calculate(self):
        for name, flag in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name + '.txt', 'w')
                for x, y, p, u, v in zip(self.vertices[vnMask, 0], self.vertices[vnMask, 1], self.p[vnMask],
                                         self.u[vnMask], self.v[vnMask]):
                    self.files[name].write(
                        '%22.16e %22.16e %22.16e %22.16e  %22.16e  %22.16e\n' % (self.tt, x, y, p, u, v))


# this has not been ported to the new-style format
class LineGauges_phi(AV_base):
    def __init__(self, gaugeEndpoints={'pressure_1': ((0.5, 0.5, 0.0), (0.5, 1.8, 0.0))}, linePoints=10):

        AV_base.__init__(self)
        self.endpoints = gaugeEndpoints
        self.flags = {}
        self.linepoints = {}
        self.files = {}  #while open later
        pointFlag = 1000
        for name, (pStart, pEnd) in self.endpoints.iteritems():
            self.flags[name] = pointFlag
            p0 = np.array(pStart)
            direction = np.array(pEnd) - p0
            self.linepoints[name] = []
            for scale in np.linspace(0.0, 1.0, linePoints):
                self.linepoints[name].append(p0 + scale * direction)
            pointFlag += 1

    def attachModel(self, model, ar):
        self.model = model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.tt = model.levelModelList[-1].timeIntegration.t
        self.phi = model.levelModelList[-1].u[0].dof
        return self

    def attachAuxiliaryVariables(self, avDict):
        return self

    def calculate(self):

        for name, flag in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name + '_phi.txt', 'w')
                for x, y, phi in zip(self.vertices[vnMask, 0], self.vertices[vnMask, 1], self.phi[vnMask]):
                    self.files[name].write('%22.16e %22.16e %22.16e %22.16e\n' % (self.tt, x, y, phi))
