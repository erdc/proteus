from collections import defaultdict, OrderedDict
from itertools import product

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from . import Comm
from .AuxiliaryVariables import AV_base
from .Profiling import logEvent as log
from . import Geom

def PointGauges(gauges, activeTime=None, sampleRate=0, fileName='point_gauges.csv'):
    """ Create a set of point gauges that will automatically be serialized as CSV data to the requested file.

    :param gauges: An iterable of "gauges".  Each gauge is specified by a 2-tuple, with the first element in the
    tuple a set of fields to be monitored, and the second element a tuple of the 3-space representations of the gauge
    locations.

    See the Gauges class for an explanation of the other parameters.

    Example:

    p = PointGauges(gauges=((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))),
                            (('p',), ((0.5, 0.5, 0),))),
                    activeTime=(0, 2.5),
                    sampleRate=0.2,
                    fileName='combined_gauge_0_0.5_sample_all.csv')

    This creates a PointGauges object that will monitor the u and v fields at the locations [0.5, 0.5,
    0] and [1, 0.5, 0], and the p field at [0.5, 0.5, 0] at simulation time between = 0 and 2.5 with samples
    taken no more frequently than every 0.2 seconds.  Results will be saved to:
    combined_gauge_0_0.5_sample_all.csv.
    """

    # build up dictionary of location information from gauges
    # dictionary of dictionaries, outer dictionary is keyed by location (3-tuple)
    # inner dictionaries contain monitored fields, and closest node
    # closest_node is None if this process does not own the node
    points = OrderedDict()
    fields = list()

    for gauge in gauges:
        gauge_fields, gauge_points = gauge
        for field in gauge_fields:
            if field not in fields:
                fields.append(field)
        for point in gauge_points:
            # initialize new dictionary of information at this location
            if point not in points:
                l_d = {'fields': set()}
                points[point] = l_d
            # add any currently unmonitored fields
            points[point]['fields'].update(gauge_fields)
    return Gauges(fields, activeTime, sampleRate, fileName, points=points)


def LineGauges(gauges, activeTime=None, sampleRate=0, fileName='line_gauges.csv'):
    """ Create a set of line gauges that will automatically be serialized as CSV data to the requested file.

    :param gauges: An iterable of "gauges".  Each gauge is specified by a 2-tuple, with the first element in the
    tuple a set of fields to be monitored, and the second element a list of pairs of endpoints of the gauges in
    3-space representation.

    See the Gauges class for an explanation of the other parameters.
    """

    # expand the product of fields and lines for each gauge

    lines = list()
    fields = list()
    for gauge in gauges:
        gauge_fields, gauge_lines = gauge
        for field in gauge_fields:
            if field not in fields:
                fields.append(field)
        lines.extend(product(gauge_fields, gauge_lines))

    return Gauges(fields, activeTime, sampleRate, fileName, lines=lines)


class Gauges(AV_base):
    """ Monitor fields at specific values.

    This class provides a generic point-wise and line-integral monitor that can be instantiated and attached to
    Proteus simulations by including them in the list of Auxiliary Variables in problem setup.

    Each Gauges instance may contain one or more fields, which may contain one or more locations to
    monitor.  The monitoring is defined over a given time and sample rate, and a filename is also supplied.  All
    results are serialized to a CSV file.

    Parallel Implementation Notes:
    After the gauge has been attached, all processes are partitioned into Gauge Owners and non-Gauge Owners.  The
    calculate method is a "no-op" for non-Owners.  For Gauge Owners, all values are computed individually,
    then collectively transmitted to the "root" process, which is the only process responsible for serializing gauge
    results to disk.  This code has not been aggressively vetted for parallel correctness or scalability.

    """
    def __init__(self, fields, activeTime=None, sampleRate=0, fileName='gauges.csv', points=None, lines=None):
        """
        Create a set of gauges that will automatically be serialized as CSV data to the requested file.


        :param activeTime: If not None, a 2-tuple of start time and end time for which the point gauge is active.
        :param sampleRate: The intervals at which samples should be measured.  Note that this is a rough lower
        bound, and that the gauge values could be computed less frequently depending on the time integrator.  The
        default value of zero computes the gauge values at every time step.
        :param fileName: The name of the file to serialize results to.

        Data is currently column-formatted, with 10 characters allotted to the time field, and 45 characters
        allotted to each point field.
        """

        AV_base.__init__(self)

        self.fields = fields
        self.activeTime = activeTime
        self.sampleRate = sampleRate
        self.fileName = fileName
        self.points = points if points else OrderedDict()
        self.lines = lines if lines else []
        self.file = None  # only the root process should have a file open
        self.flags = {}
        self.files = {}
        self.outputWriterReady = False
        self.last_output = None
        self.pointGaugeMats = []
        self.field_ids = []
        self.dofsVecs = []
        self.pointGaugeVecs = []
        self.segments = []

    def findNearestNode(self, location):
        """Given a gauge location, attempts to locate the most suitable process for monitoring information about
        this location, as well as the node on the process closest to the location.

        Returns a 2-tuple containing an identifier for the closest 'owning' process as well as the local id of the
        node.
        """

        # determine local nearest node distance
        node_distances = np.linalg.norm(self.vertices - location, axis=1)
        nearest_node = np.argmin(node_distances)
        nearest_node_distance = node_distances[nearest_node]

        # determine global nearest node
        comm = Comm.get().comm.tompi4py()
        global_min_distance, owning_proc = comm.allreduce(nearest_node_distance, op=MPI.MINLOC)
        log("Gauges at location: [%g %g %g] assigned to %d" % (location[0], location[1], location[2], owning_proc), 3)
        if comm.rank != owning_proc:
            nearest_node = None

        assert owning_proc is not None
        return owning_proc, nearest_node

    def buildQuantityRow(self, m, femFun, quantity_id, quantity):
        """ Builds up contributions to gauge operator from the underlying element space
        """

        location, node = quantity

        # get the mesh for this quantity
        mesh = femFun.femSpace.mesh
        # search elements that contain the nearest node
        # use nearest node if the location is not found on any elements
        for eOffset in range(mesh.nodeElementOffsets[node], mesh.nodeElementOffsets[node + 1]):
            eN = femFun.femSpace.mesh.nodeElementsArray[eOffset]
            # evaluate the inverse map for element eN
            xi = femFun.femSpace.elementMaps.getInverseValue(eN, location)
            # query whether xi lies within the reference element
            if femFun.femSpace.elementMaps.referenceElement.onElement(xi):
                for i, psi in enumerate(femFun.femSpace.referenceFiniteElement.localFunctionSpace.basis):
                    # assign quantity weights here
                    m[quantity_id, femFun.femSpace.dofMap.l2g[eN, i]] = psi(xi)
                break
        else:
            # just use nearest node for now if we're given a point outside the domain.
            # the ideal thing would be to find the element with the nearest face
            m[quantity_id, node] = 1

    def initOutputWriter(self):
        """ Initialize communication strategy for collective output of gauge data.

        On the root process in this communicator, create a map of quantity owners and the corresponding location in
        their arrays.  This process is responsible for collecting gauge data and saving it to disk.

        Gauge data is globally ordered by field, then by location id (as ordered by globalMeasuredQuantities)
        """

        numLocalQuantities = sum([len(self.measuredQuantities[field]) for field in self.fields])
        self.localQuantitiesBuf = np.zeros(numLocalQuantities)

        if self.gaugeComm.rank != 0:
            self.globalQuantitiesBuf = None
        else:
            self.file = open(self.fileName, 'w')

            self.quantityIDs = [0] * self.gaugeComm.size
            # Assign quantity ids to processors
            for field in self.fields:
                for id in range(len(self.globalMeasuredQuantities[field])):
                    location, owningProc = self.globalMeasuredQuantities[field][id]
                    gaugeProc = self.globalGaugeRanks[owningProc]
                    quantityID = self.quantityIDs[gaugeProc]
                    self.quantityIDs[gaugeProc] += 1
                    assert gaugeProc >= 0
                    self.globalMeasuredQuantities[field][id] = location, gaugeProc, quantityID
                    log("Gauge for %s[%d] at %e %e %e is at P[%d][%d]" % (field, id, location[0], location[1],
                                                                          location[2], gaugeProc, quantityID), 5)

            log("Quantity IDs:\n%s" % str(self.quantityIDs), 5)

            numGlobalQuantities = sum([len(self.globalMeasuredQuantities[field]) for field in self.fields])
            self.globalQuantitiesBuf = np.zeros(numGlobalQuantities)

            # determine mapping from global measured quantities to communication buffers
            self.globalQuantitiesMap = np.zeros(numGlobalQuantities, dtype=np.int)
            i = 0
            for field in self.fields:
                for location, gaugeProc, quantityID in self.globalMeasuredQuantities[field]:
                    self.globalQuantitiesMap[i] = sum(self.quantityIDs[:gaugeProc]) + quantityID
                    assert self.globalQuantitiesMap[i] < numGlobalQuantities
                    i += 1

            log("Global Quantities Map: \n%s" % str(self.globalQuantitiesMap), 5)
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

        for field in self.fields:
            for location, owningProc in self.globalMeasuredQuantities[field]:
                gaugeOwners.update((owningProc,))

        self.isGaugeOwner = comm.rank in gaugeOwners
        gaugeComm = comm.Split(color=self.isGaugeOwner)

        log("Gauge owner: %d" % self.isGaugeOwner, 5)
        if self.isGaugeOwner:
            self.gaugeComm = gaugeComm
            gaugeRank = self.gaugeComm.rank
        else:
            self.gaugeComm = None
            gaugeRank = -1
        self.globalGaugeRanks = comm.allgather(gaugeRank)
        log("Gauge ranks: \n%s" % str(self.globalGaugeRanks), 5)


    def addLineGaugePoints(self, line, line_segments):
        """ Add all gauge points from each line into self.points
        """
        points = self.points

        field, endpoints = line
        for point in line_segments:
            if point in points:
                points[point]['fields'].update(field)
            else:
                points[point] = {'fields':set((field,))}

    def identifyMeasuredQuantities(self):
        """ build measured quantities, a list of fields

        Each field in turn contains a list of gauge locations and their accompanying nearest node
        only local quantities are saved
        """

        self.measuredQuantities = defaultdict(list)
        self.globalMeasuredQuantities = defaultdict(list)

        points = self.points

        for point, l_d in points.iteritems():
            owningProc, nearestNode = self.findNearestNode(point)
            l_d['nearest_node'] = nearestNode
            for field in l_d['fields']:
                self.globalMeasuredQuantities[field].append((point, owningProc))
                if nearestNode is not None:
                    point_id = len(self.measuredQuantities[field])
                    log("Gauge for %s[%d] at %e %e %e is closest to node %d" % (field, point_id, point[0], point[1],
                                                                                point[2], nearestNode), 3)
                    l_d[field] = point_id
                    self.measuredQuantities[field].append((point, nearestNode))

    def buildPointGaugeOperators(self):
        """ Build the linear algebra operators needed to compute the point gauges.

        The operators are all local since the point gauge measurements are calculated locally.
        """

        for field, field_id in zip(self.fields, self.field_ids):
            m = PETSc.Mat().create(PETSc.COMM_SELF)
            m.setSizes([len(self.measuredQuantities[field]), self.num_owned_nodes])
            m.setType('aij')
            m.setUp()
            # matrices are a list in same order as fields
            self.pointGaugeMats.append(m)
            # dofs are a list in same order as fields as well
            dofs = self.u[field_id].dof
            dofsVec = PETSc.Vec().createWithArray(dofs, comm=PETSc.COMM_SELF)
            self.dofsVecs.append(dofsVec)

        for field, field_id, m in zip(self.fields, self.field_ids, self.pointGaugeMats):
            # get the FiniteElementFunction object for this quantity
            femFun = self.u[field_id]
            for quantity_id, quantity in enumerate(self.measuredQuantities[field]):
                location, node = quantity
                log("Gauge for: %s at %e %e %e is on local operator row %d" % (field, location[0], location[1],
                                                                      location[2], quantity_id), 3)
                self.buildQuantityRow(m, femFun, quantity_id, quantity)
            pointGaugesVec = PETSc.Vec().create(comm=PETSc.COMM_SELF)
            pointGaugesVec.setSizes(len(self.measuredQuantities[field]))
            pointGaugesVec.setUp()
            self.pointGaugeVecs.append(pointGaugesVec)

        for m in self.pointGaugeMats:
            m.assemble()

    def getMeshIntersections(self, line):
        field, endpoints = line
        # get Proteus mesh index for this field
        field_id = self.fieldNames.index(field)
        femFun = self.u[field_id]
        mesh = femFun.femSpace.mesh
        referenceElement = femFun.femSpace.elementMaps.referenceElement
        if referenceElement.dim == 2 and referenceElement.nNodes == 3:
            toPolyhedron = Geom.triangleVerticesToNormals
        elif referenceElement.dim == 3 and referenceElement.nNodes == 4:
            toPolyhedron = Geom.tetrahedronVerticesToNormals
        else:
            raise NotImplementedError("Unable to compute mesh intersections for this element type")
        return Geom.getMeshIntersections(mesh, toPolyhedron, endpoints)

    def buildLineGaugeOperators(self, lines, linesSegments):
        """ Build the linear algebra operators needed to compute the line gauges.

        The operators are collective since the line gauge measurements may be across multiple processors.

        Unlike point gauges, each line has its own communicator and field associated with it
        """

        #create lineGaugesVec to store contributions to all lines from this process
        self.lineGaugesVec = PETSc.Vec().create(comm=PETSc.COMM_SELF)
        self.lineGaugesVec.setSizes(len(lines))
        self.lineGaugesVec.setUp()

        # create lineGaugeMats to store coefficients mapping contributions from each field
        # to the line gauges
        self.lineGaugeMats = []
        # size of lineGaugeMats depends on number of local points for each field
        for pointGaugesVec in self.pointGaugeVecs:
            m = PETSc.Mat().create(comm=PETSc.COMM_SELF)
            m.setSizes([len(lines), pointGaugesVec.getSize()])
            m.setType('aij')
            m.setUp()
            self.lineGaugeMats.append(m)

        # Assemble contributions from each point in each line segment

        for lineIndex, (line, segments) in enumerate(zip(self.lines, linesSegments)):
            field, endpoints = line
            fieldIndex = self.fields.index(field)

            # Trapezoid Rule to calculate coefficients here
            for p1, p2 in zip(segments[:-1], segments[1:]):
                segmentLength = np.linalg.norm(np.asarray(p2)-np.asarray(p1))
                for point in p1, p2:
                    point_data = self.points[tuple(point)]
                    # only assign coefficients for locally owned points
                    if field in point_data:
                        pointID = point_data[field]
                        self.lineGaugeMats[fieldIndex].setValue(lineIndex, pointID, segmentLength/2, addv=True)

        for m in self.lineGaugeMats:
            m.assemble()

    def attachModel(self, model, ar):
        """ Attach this gauge to the given simulation model.
        """

        self.model = model
        self.fieldNames = model.levelModelList[-1].coefficients.variableNames
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.num_owned_nodes = model.levelModelList[-1].mesh.nNodes_global
        self.u = model.levelModelList[-1].u
        self.timeIntegration = model.levelModelList[-1].timeIntegration

        for field in self.fields:
            field_id = self.fieldNames.index(field)
            self.field_ids.append(field_id)

        linesSegments = []
        for line in self.lines:
            lineSegments = self.getMeshIntersections(line)
            self.addLineGaugePoints(line, lineSegments)
            linesSegments.append(np.asarray(list(lineSegments)))

        self.identifyMeasuredQuantities()

        self.buildGaugeComm()

        if self.isGaugeOwner:
            self.initOutputWriter()
            self.buildPointGaugeOperators()
            self.buildLineGaugeOperators(self.lines, linesSegments)
            self.outputHeader()
        return self

    def get_time(self):
        """ Returns the current model time"""
        return self.timeIntegration.t

    def outputHeader(self):
        """ Outputs a single header for a CSV style file to self.file"""

        assert self.isGaugeOwner

        if self.gaugeComm.rank == 0:
            self.file.write("%10s" % ('time',))
            for field in self.fields:
                for quantity in self.globalMeasuredQuantities[field]:
                    location, gaugeProc, quantityID = quantity
                    self.file.write(",%12s [%9.5g %9.5g %9.5g]" % (field, location[0], location[1], location[2]))
            self.file.write('\n')

    def outputRow(self, time):
        """ Outputs a single row of currently calculated gauge data to self.file"""

        assert self.isGaugeOwner

        self.localQuantitiesBuf = np.concatenate([gaugesVec.getArray() for gaugesVec in self.pointGaugeVecs])
        self.gaugeComm.Gatherv(self.localQuantitiesBuf, self.globalQuantitiesBuf)

        if self.lines:
            lineGaugeBuf = self.lineGaugesVec.getArray()
            globalGaugeBuf = lineGaugeBuf.copy()
            self.gaugeComm.Reduce(lineGaugeBuf, globalGaugeBuf, op=MPI.SUM)
        else:
            globalGaugeBuf = []

        if self.gaugeComm.rank == 0:
            self.file.write("%10.4e" % time)
            for id in self.globalQuantitiesMap:
                self.file.write(", %43.18e" % (self.globalQuantitiesBuf[id],))
            for lineGauge in globalGaugeBuf:
                self.file.write(", %43.18e" % (lineGauge))
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

        # check that gauge is in its active time region
        if self.activeTime is not None and (self.activeTime[0] > time or self.activeTime[1] < time):
            return

        # check that gauge is ready to be sampled again
        if self.last_output is not None and time < self.last_output + self.sampleRate:
            return

        for m, dofsVec, gaugesVec in zip(self.pointGaugeMats, self.dofsVecs, self.pointGaugeVecs):
            m.mult(dofsVec, gaugesVec)

        # this could be optimized out... but why?
        self.lineGaugesVec.zeroEntries()
        for m, dofsVec in zip(self.lineGaugeMats, self.pointGaugeVecs):
            m.multAdd(dofsVec, self.lineGaugesVec, self.lineGaugesVec)

        self.outputRow(time)