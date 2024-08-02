"""Auxiliary variable classes for storing solutions at points and
along lines to mimic gauges in lab and field experiments.

.. inheritance-diagram:: proteus.Gauges
   :parts: 1

"""
import os
from collections import defaultdict, OrderedDict
from itertools import product

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
from numpy.linalg import norm

from . import Comm
from .AuxiliaryVariables import AV_base
from .Profiling import logEvent
from proteus.MeshTools import triangleVerticesToNormals, tetrahedronVerticesToNormals, getMeshIntersections
from proteus import Profiling


def PointGauges(gauges, activeTime=None, sampleRate=0, fileName='point_gauges.csv'):
    """Create a set of point gauges that will automatically be serialized
    as CSV data to the requested file.

    :param gauges: An iterable of "gauges".  Each gauge is specified
                   by a 2-tuple, with the first element in the tuple a
                   set of fields to be monitored, and the second
                   element a tuple of the 3-space representations of
                   the gauge locations.

    See the Gauges class for an explanation of the other parameters.

    Example::

      p = PointGauges(gauges=((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))),
                             (('p',), ((0.5, 0.5, 0),))),
                      activeTime=(0, 2.5),
                      sampleRate=0.2,
                      fileName='combined_gauge_0_0.5_sample_all.csv')

    This creates a PointGauges object that will monitor the u and v
    fields at the locations [0.5, 0.5, 0] and [1, 0.5, 0], and the p
    field at [0.5, 0.5, 0] at simulation time between = 0 and 2.5 with
    samples taken no more frequently than every 0.2 seconds.  Results
    will be saved to: combined_gauge_0_0.5_sample_all.csv.

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
    """Create a set of line gauges that will automatically be serialized
    as CSV data to the requested file.  The line gauges will gather
    data at every element on the mesh between the two endpoints on
    each line.

    :param gauges: An iterable of "gauges".  Each gauge is specified
                   by a 2-tuple, with the first element in the tuple a
                   set of fields to be monitored, and the second
                   element a list of pairs of endpoints of the gauges
                   in 3-space representation.

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

def LineIntegralGauges(gauges, activeTime=None, sampleRate=0, fileName='line_integral_gauges.csv'):
    """Create a set of line integral gauges that will automatically be
    serialized as CSV data to the requested file.

    :param gauges: An iterable of "gauges".  Each gauge is specified
                   by a 2-tuple, with the first element in the tuple a
                   set of fields to be monitored, and the second
                   element a list of pairs of endpoints of the gauges
                   in 3-space representation.

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

    return Gauges(fields, activeTime, sampleRate, fileName, lines=lines, integrate=True)


class Gauges(AV_base):
    """Monitor fields at specific values.

    This class provides a generic point-wise and line-integral monitor
    that can be instantiated and attached to Proteus simulations by
    including them in the list of Auxiliary Variables in problem
    setup.

    Each Gauges instance may contain one or more fields, which may
    contain one or more locations to monitor.  The monitoring is
    defined over a given time and sample rate, and a filename is also
    supplied.  All results are serialized to a CSV file.

    Parallel Implementation Notes: After the gauge has been attached,
    all processes are partitioned into Gauge Owners and non-Gauge
    Owners.  The calculate method is a "no-op" for non-Owners.  For
    Gauge Owners, all values are computed individually, then
    collectively transmitted to the "root" process, which is the only
    process responsible for serializing gauge results to disk.  This
    code has not been aggressively vetted for parallel correctness or
    scalability.

    """
    def __init__(self, fields, activeTime=None, sampleRate=0, fileName='gauges.csv', points=None, lines=None,
                 integrate=False):
        """Create a set of gauges that will automatically be serialized as
        CSV data to the requested file.


        :param activeTime: If not None, a 2-tuple of start time and
                           end time for which the point gauge is
                           active.
        :param sampleRate: The intervals at which samples should be
                           measured.  Note that this is a rough lower
                           bound, and that the gauge values could be
                           computed less frequently depending on the
                           time integrator.  The default value of zero
                           computes the gauge values at every time
                           step.
        :param fileName: The name of the file to serialize results to.

        Data is currently column-formatted, with 10 characters
        allotted to the time field, and 45 characters allotted to each
        point field.

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
        self.adapted = False

        self.isPointGauge = bool(points)
        self.isLineGauge = bool(lines) and not integrate
        self.isLineIntegralGauge = bool(lines) and integrate

        if not (self.isPointGauge or self.isLineGauge or self.isLineIntegralGauge):
            raise ValueError("Need to provide points or lines")
        if sum((self.isPointGauge, self.isLineGauge, self.isLineIntegralGauge)) > 1:
            raise ValueError("Must be one of point or line gauge but not both")

    def getLocalNearestNode(self, location):
        # determine local nearest node distance
        nearest_node_distance_kdtree, nearest_node_kdtree = self.nodes_kdtree.query(location)
        comm = Comm.get().comm.tompi4py()
        return comm.rank, nearest_node_kdtree, nearest_node_distance_kdtree

    def getLocalElement(self, femSpace, location, node):
        """Given a location and its nearest node, determine if it is on a
        local element.

        Returns None if location is not on any elements owned by this
        process

        """

        # search elements that contain the nearest node
        patchBoundaryNodes=set()
        checkedElements=[]
        for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
            eN = femSpace.mesh.nodeElementsArray[eOffset]
            checkedElements.append(eN)
            patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
            # evaluate the inverse map for element eN
            xi = femSpace.elementMaps.getInverseValue(eN, location)
            # query whether xi lies within the reference element
            if femSpace.elementMaps.referenceElement.onElement(xi):
                return eN
        for node in patchBoundaryNodes:
            for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
                eN = femSpace.mesh.nodeElementsArray[eOffset]
                if eN not in checkedElements:
                    checkedElements.append(eN)
                    # evaluate the inverse map for element eN
                    xi = femSpace.elementMaps.getInverseValue(eN, location)
                    # query whether xi lies within the reference element
                    if femSpace.elementMaps.referenceElement.onElement(xi):
                        return eN
        # no elements found
        return None

    def findNearestNode(self, femSpace, location):
        """Given a gauge location, attempts to locate the most suitable
        process for monitoring information about this location, as
        well as the node on the process closest to the location.

        Returns a 2-tuple containing an identifier for the closest
        'owning' process as well as the local ids of the node and
        nearest element.

        """
        comm = Comm.get().comm.tompi4py()
        comm_rank, nearest_node, nearest_node_distance = self.getLocalNearestNode(location)
        local_element = self.getLocalElement(femSpace, location, nearest_node)

        # determine global nearest node
        haveElement = int(local_element is not None)
        global_have_element, owning_proc = comm.allreduce((haveElement, comm.rank),
                                                          op=MPI.MAXLOC)
        if global_have_element:
            logEvent("Gauges on element at location: [%g %g %g] assigned to %d" % (location[0], location[1], location[2],
                                                                    owning_proc), 3)
        else:
            # gauge isn't on any of the elements, just use nearest node
            global_min_distance, owning_proc = comm.allreduce((nearest_node_distance,comm.rank), op=MPI.MINLOC)
            logEvent("Off-element gauge location: [%g %g %g] assigned to %d" % (location[0], location[1], location[2],
                                                                 owning_proc), 3)
        if comm.rank != owning_proc:
            nearest_node = None

        assert owning_proc is not None
        return owning_proc, nearest_node

    def buildQuantityRow(self, m, femFun, quantity_id, quantity):
        """Builds up contributions to gauge operator from the underlying
        element space

        """

        location, node = quantity

        # search elements that contain the nearest node
        # use nearest node if the location is not found on any elements
        localElement = self.getLocalElement(femFun.femSpace, location, node)
        if localElement is not None:
            for i, psi in enumerate(femFun.femSpace.referenceFiniteElement.localFunctionSpace.basis):
                # assign quantity weights here
                xi = femFun.femSpace.elementMaps.getInverseValue(localElement, location)
                m[quantity_id, femFun.femSpace.dofMap.l2g[localElement, i]] = psi(xi)
        else:
            # just use nearest node for now if we're given a point outside the domain.
            # the ideal thing would be to find the element with the nearest face
            m[quantity_id, node] = 1


    def initOutputWriter(self):
        """Initialize communication strategy for collective output of gauge
        data.

        On the root process in this communicator, create a map of
        quantity owners and the corresponding location in their
        arrays.  This process is responsible for collecting gauge data
        and saving it to disk.

        Gauge data is globally ordered by field, then by location id
        (as ordered by globalMeasuredQuantities)

        """

        numLocalQuantities = sum([len(self.measuredQuantities[field]) for field in self.fields])
        self.localQuantitiesBuf = np.zeros(numLocalQuantities)
        if self.gaugeComm.rank != 0:
            self.globalQuantitiesBuf = None
            self.globalQuantitiesCounts = None
        else:
            if self.adapted:
              if(Profiling.logDir not in self.fileName):
                self.fileName = os.path.join(Profiling.logDir, self.fileName)
              self.file = open(self.fileName, 'a')
            else:
              self.fileName = os.path.join(Profiling.logDir, self.fileName)
              self.file = open(self.fileName, 'w')

            if self.isLineIntegralGauge:
                #Only need to set up mapping for point gauges
                return

            quantityIDs = [0] * self.gaugeComm.size
            numGlobalQuantities = sum([len(self.globalMeasuredQuantities[field]) for field in self.fields])
            # Assign quantity ids to processors
            for field in self.fields:
                for id in range(len(self.globalMeasuredQuantities[field])):
                    location, owningProc = self.globalMeasuredQuantities[field][id]
                    gaugeProc = self.globalGaugeRanks[owningProc]
                    quantityID = quantityIDs[gaugeProc]
                    quantityIDs[gaugeProc] += 1
                    assert gaugeProc >= 0
                    self.globalMeasuredQuantities[field][id] = location, gaugeProc, quantityID
                    logEvent("Gauge for %s[%d] at %e %e %e is at P[%d][%d]" % (field, id, location[0], location[1],
                                                                          location[2], gaugeProc, quantityID), 5)

            logEvent("Quantity IDs:\n%s" % str(quantityIDs), 5)

            # determine mapping from global measured quantities to communication buffers
            self.globalQuantitiesMap = np.zeros(numGlobalQuantities, dtype=int)
            i = 0
            for field in self.fields:
                for location, gaugeProc, quantityID in self.globalMeasuredQuantities[field]:
                    self.globalQuantitiesMap[i] = sum(quantityIDs[:gaugeProc]) + quantityID
                    assert self.globalQuantitiesMap[i] < numGlobalQuantities
                    i += 1

            # a couple consistency checks
            assert sum(quantityIDs) == numGlobalQuantities
            assert all(quantityID > 0 for quantityID in quantityIDs)

            # final ids also equal to the counts on each process
            self.globalQuantitiesCounts = quantityIDs
            self.globalQuantitiesBuf = np.zeros(numGlobalQuantities, dtype=np.double)
            logEvent("Global Quantities Map: \n%s" % str(self.globalQuantitiesMap), 5)

        self.outputWriterReady = True

    def buildGaugeComm(self):
        """Create a communicator composed only of processes that own gauge
        quantities.

        Collective over global communicator.  Builds a local
        communicator for collecting all gauge data.  This communicator
        contains only processes that will contain gauge data.

        """

        comm = Comm.get().comm.tompi4py()

        gaugeOwners = set()

        for field in self.fields:
            for location, owningProc in self.globalMeasuredQuantities[field]:
                gaugeOwners.update((owningProc,))

        self.isGaugeOwner = comm.rank in gaugeOwners
        gaugeComm = comm.Split(color=self.isGaugeOwner)

        logEvent("Gauge owner: %d" % self.isGaugeOwner, 5)
        if self.isGaugeOwner:
            self.gaugeComm = gaugeComm
            gaugeRank = self.gaugeComm.rank
        else:
            self.gaugeComm = None
            gaugeRank = -1
        self.globalGaugeRanks = comm.allgather(gaugeRank)
        logEvent("Gauge ranks: \n%s" % str(self.globalGaugeRanks), 5)


    def addLineGaugePoints(self, line, line_segments):
        """Add all gauge points from each line into self.points
        """
        points = self.points
        new_points = {}
        field, endpoints = line
        comm = Comm.get().comm.tompi4py()

        def addPoint(points, field, point):
            point = tuple(point)
            if point in points:
                if self.isLineIntegralGauge:
                    no_output = points[point]['no_output'] if 'no_output' in points[point] else set()
                    points[point]['no_output'] = no_output.union(set((field,)) - points[point]['fields'])
                points[point]['fields'].update((field,))
            else:
                ignore1, nearestNode, ignore2 = self.getLocalNearestNode(point)
                if self.isLineIntegralGauge:
                    points[point] = {'fields':set((field,)), 'no_output': set((field,)),
                                     'nearest_node': nearestNode,
                                     'owning_proc': comm.rank}
                else:
                    points[point] = {'fields':set((field,)),
                                     'nearest_node': nearestNode,
                                     'owning_proc': comm.rank}
            new_points[point] = points[point]

        for segment in line_segments:
            logEvent("Processing segment [ %e %e %e ] to [ %e %e %e ]" % (
                segment[0][0], segment[0][1], segment[0][2],
                segment[1][0], segment[1][1], segment[1][2]), 5)
            startPoint, endPoint = segment
            # only add both sides of segment to line integral gauges and first segment
            if self.isLineIntegralGauge or all(startPoint == endpoints[0]):
                addPoint(points, field, startPoint)
            addPoint(points, field, endPoint)

        if self.isLineGauge:
            new_points = comm.gather(new_points)
            if comm.rank == 0:
                for new_points_i in new_points:
                    points.update(new_points_i)
                # resort points
                points = OrderedDict(sorted(points.items()))
            self.points = comm.bcast(points)

    def identifyMeasuredQuantities(self):
        """ build measured quantities, a list of fields

        Each field in turn contains a list of gauge locations and their accompanying nearest node
        only local quantities are saved
        """

        self.measuredQuantities = defaultdict(list)
        self.globalMeasuredQuantities = defaultdict(list)
        comm = Comm.get().comm.tompi4py()

        points = self.points

        for point, l_d in points.items():
            if 'nearest_node' not in l_d:
                # TODO: Clarify assumption here about all fields sharing the same element mesh
                field_id = self.fieldNames.index(list(l_d['fields'])[0])
                femSpace = self.u[field_id].femSpace
                owningProc, nearestNode = self.findNearestNode(femSpace, point)
                l_d['nearest_node'] = nearestNode
            else:
                owningProc = l_d['owning_proc']
                # nearestNode only makes sense on owning process
                # so even if we have this information, it's not valid for this point
                if owningProc == comm.rank:
                    nearestNode = l_d['nearest_node']
                else:
                    nearestNode = None
            for field in l_d['fields']:
                self.globalMeasuredQuantities[field].append((point, owningProc))
                if nearestNode is not None:
                    point_id = len(self.measuredQuantities[field])
                    logEvent("Gauge for %s[%d] at %e %e %e is closest to node %d" % (field, point_id, point[0], point[1],
                                                                                point[2], nearestNode), 3)
                    l_d[field] = point_id
                    self.measuredQuantities[field].append((point, nearestNode))

    def buildPointGaugeOperators(self):
        """ Build the linear algebra operators needed to compute the point gauges.

        The operators are all local since the point gauge measurements are calculated locally.
        """

        for field, field_id in zip(self.fields, self.field_ids):
            m = PETSc.Mat().create(PETSc.COMM_SELF)
            m.setSizes([len(self.measuredQuantities[field]),
                        self.u[field_id].femSpace.dim])
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
                logEvent("Gauge for: %s at %e %e %e is on local operator row %d" % (field, location[0], location[1],
                                                                      location[2], quantity_id), 3)
                self.buildQuantityRow(m, femFun, quantity_id, quantity)
            pointGaugesVec = PETSc.Vec().create(comm=PETSc.COMM_SELF)
            pointGaugesVec.setSizes(len(self.measuredQuantities[field]))
            pointGaugesVec.setUp()
            self.pointGaugeVecs.append(pointGaugesVec)

        for m in self.pointGaugeMats:
            m.assemble()

    def pruneDuplicateSegments(self, endpoints, length_segments):
        """ prune duplicate segments across processors

        endpoints - a pair of points in 3-space defining the line
        length_segments - a pair of intersections augmented by length

        this could be optimized
        """

        eps = 1e-4

        comm = Comm.get().comm.tompi4py()

        length_segments = sorted(length_segments)
        length_segments = comm.gather(length_segments)


        if comm.rank != 0:
            selected_segments = None
        else:
            selected_segments = [[] for i in range(len(length_segments))]
            segment_pos = 0
            while segment_pos < (1 - eps):
                # choose the longest line from those that start at segment_pos
                longest_segment = 0, None, None
                for proc_rank, proc_length_segments in enumerate(length_segments):
                    segment_id = 0
                    for segment_id, length_segment in enumerate(proc_length_segments):
                        # ignore segments below current position (they will be discarded)
                        start, end, segment = length_segment
                        if start < (segment_pos - eps):
                            continue
                        # equality test
                        elif start < (segment_pos + eps):
                            segment_length = end - start
                            if segment_length > longest_segment[0]:
                                longest_segment = segment_length, proc_rank, segment
                        else:
                            break
                    # discard any segments that start before our current position
                    proc_length_segments[:] = proc_length_segments[segment_id:]

                segment_length, proc_rank, segment = longest_segment
                if segment_length == 0:
                    print(segment_pos)
                    print('segments')
                    for segment in selected_segments: print(segment)
                    print('length_segments')
                    for length_segment in length_segments: print(length_segment)
                    raise FloatingPointError("Unable to identify next segment while segmenting, are %s in domain?" %
                                             str(endpoints))
                logEvent("Identified best segment of length %g on %d: %s" % (segment_length, proc_rank, str(segment)), 9)
                selected_segments[proc_rank].append(segment)
                segment_pos += segment_length

            err = abs(segment_pos - 1)
            if err > 1e-8:
                msg = "Segmented line %s different from original length by ratio %e\n segments: %s" % (
                    str(endpoints), err, str(selected_segments))
                logEvent(msg, 3)
                if err > 10*eps:
                    raise FloatingPointError(msg)

        logEvent("Selected segments: %s" % str(selected_segments), 9)
        segments = comm.scatter(selected_segments)
        return segments

    def getMeshIntersections(self, line):
        field, endpoints = line
        # get Proteus mesh index for this field
        field_id = self.fieldNames.index(field)
        femFun = self.u[field_id]
        mesh = femFun.femSpace.mesh
        referenceElement = femFun.femSpace.elementMaps.referenceElement
        if referenceElement.dim == 2 and referenceElement.nNodes == 3:
            toPolyhedron = triangleVerticesToNormals
        elif referenceElement.dim == 3 and referenceElement.nNodes == 4:
            toPolyhedron = tetrahedronVerticesToNormals
        else:
            raise NotImplementedError("Unable to compute mesh intersections for this element type")
        intersections = np.asarray(list(getMeshIntersections(mesh, toPolyhedron, endpoints)), dtype=np.double)
        endpoints = np.asarray(endpoints, np.double)
        length = norm(endpoints[1] - endpoints[0])

        length_segments = [(norm(i[0]-endpoints[0])/length, norm(i[1]-endpoints[0])/length, i) for i in intersections]

        segments = self.pruneDuplicateSegments(endpoints, length_segments)
        return segments


    def buildLineIntegralGaugeOperators(self, lines, linesSegments):
        """ Build the linear algebra operators needed to compute the line integral gauges.

        The operators are local to each process, contributions are currently summed in the output functions.
        """

        #create lineIntegralGaugesVec to store contributions to all lines from this process
        self.lineIntegralGaugesVec = PETSc.Vec().create(comm=PETSc.COMM_SELF)
        self.lineIntegralGaugesVec.setSizes(len(lines))
        self.lineIntegralGaugesVec.setUp()

        # create lineIntegralGaugeMats to store coefficients mapping contributions from each field
        # to the line integral gauges
        self.lineIntegralGaugeMats = []

        if not self.isLineIntegralGauge:
            return

        # size of lineIntegralGaugeMats depends on number of local points for each field
        for pointGaugesVec in self.pointGaugeVecs:
            m = PETSc.Mat().create(comm=PETSc.COMM_SELF)
            m.setSizes([len(lines), pointGaugesVec.getSize()])
            m.setType('aij')
            m.setUp()
            self.lineIntegralGaugeMats.append(m)

        # Assemble contributions from each point in each line segment

        for lineIndex, (line, segments) in enumerate(zip(self.lines, linesSegments)):
            field, endpoints = line
            fieldIndex = self.fields.index(field)

            # Trapezoid Rule to calculate coefficients here
            for p1, p2 in segments:
                segmentLength = np.linalg.norm(np.asarray(p2)-np.asarray(p1))
                for point in p1, p2:
                    point_data = self.points[tuple(point)]
                    # only assign coefficients for locally owned points
                    if field in point_data:
                        pointID = point_data[field]
                        self.lineIntegralGaugeMats[fieldIndex].setValue(lineIndex, pointID, segmentLength/2, addv=True)

        for m in self.lineIntegralGaugeMats:
            m.assemble()

    def attachModel(self, model, ar):
        """ Attach this gauge to the given simulation model.
        """
        from scipy import spatial
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
        self.nodes_kdtree = spatial.cKDTree(model.levelModelList[-1].mesh.nodeArray)
        linesSegments = []
        for line in self.lines:
            lineSegments = self.getMeshIntersections(line)
            self.addLineGaugePoints(line, lineSegments)
            linesSegments.append(lineSegments)

        self.identifyMeasuredQuantities()

        self.buildGaugeComm()

        if self.isGaugeOwner:
            self.initOutputWriter()
            self.buildPointGaugeOperators()
            self.buildLineIntegralGaugeOperators(self.lines, linesSegments)
            if self.adapted:
              pass
            else:
              self.outputHeader()
        return self

    def get_time(self):
        """ Returns the current model time"""

        return self.timeIntegration.tLast

    def outputHeader(self):
        """ Outputs a single header for a CSV style file to self.file"""

        assert self.isGaugeOwner
        if self.gaugeComm.rank == 0:
            self.file.write("%10s" % ('time',))
            if self.isPointGauge or self.isLineGauge:
                for field in self.fields:
                    for quantity in self.globalMeasuredQuantities[field]:
                        location, gaugeProc, quantityID = quantity
                        self.file.write(",%12s [%9.5g %9.5g %9.5g]" % (field, location[0], location[1], location[2]))
            elif self.isLineIntegralGauge:
                for line in self.lines:
                    self.file.write(",%12s [%9.5g %9.5g %9.5g] - [%9.5g %9.5g %9.5g]" % (
                        line[0], line[1][0][0], line[1][0][1], line[1][0][2],
                                 line[1][1][0], line[1][1][1], line[1][1][2]))
            self.file.write('\n')

    def outputRow(self, time):
        """ Outputs a single row of currently calculated gauge data to self.file"""

        assert self.isGaugeOwner

        if self.isPointGauge or self.isLineGauge:
            self.localQuantitiesBuf = np.concatenate([gaugesVec.getArray() for gaugesVec in
                                                      self.pointGaugeVecs]).astype(np.double)
            logEvent("Sending local array of type %s and shape %s to root on comm %s" % (
                str(self.localQuantitiesBuf.dtype), str(self.localQuantitiesBuf.shape), str(self.gaugeComm)), 9)
            if self.gaugeComm.rank == 0:
                logEvent("Receiving global array of type %s and shape %s on comm %s" % (
                str(self.localQuantitiesBuf.dtype), str(self.globalQuantitiesBuf.shape), str(self.gaugeComm)), 9)
            self.gaugeComm.Gatherv(sendbuf=[self.localQuantitiesBuf, MPI.DOUBLE],
                                   recvbuf=[self.globalQuantitiesBuf, (self.globalQuantitiesCounts, None),
                                            MPI.DOUBLE], root=0)
            self.gaugeComm.Barrier()
        if self.isLineIntegralGauge:
            lineIntegralGaugeBuf = self.lineIntegralGaugesVec.getArray()
            globalLineIntegralGaugeBuf = lineIntegralGaugeBuf.copy()
            self.gaugeComm.Reduce(lineIntegralGaugeBuf, globalLineIntegralGaugeBuf, op=MPI.SUM)
        else:
            globalLineIntegralGaugeBuf = []

        if self.gaugeComm.rank == 0:
            self.file.write("%25.15e" % time)
            if self.isPointGauge or self.isLineGauge:
                for id in self.globalQuantitiesMap:
                    self.file.write(", %43.18e" % (self.globalQuantitiesBuf[id],))
            if self.isLineIntegralGauge:
                for lineIntegralGauge in globalLineIntegralGaugeBuf:
                    self.file.write(", %80.18e" % (lineIntegralGauge))
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
        logEvent("Gauges calculate called at time %g" % time, 4)
        # check that gauge is in its active time region
        if self.activeTime is not None and (self.activeTime[0] > time or self.activeTime[1] < time):
            return

        # check that gauge is ready to be sampled again
        if self.last_output is not None and time < self.last_output + self.sampleRate:
            return

        for m, dofsVec, gaugesVec in zip(self.pointGaugeMats, self.dofsVecs, self.pointGaugeVecs):
            m.mult(dofsVec, gaugesVec)

        # this could be optimized out... but why?
        self.lineIntegralGaugesVec.zeroEntries()
        for m, dofsVec in zip(self.lineIntegralGaugeMats, self.pointGaugeVecs):
            m.multAdd(dofsVec, self.lineIntegralGaugesVec, self.lineIntegralGaugesVec)

        self.outputRow(time)