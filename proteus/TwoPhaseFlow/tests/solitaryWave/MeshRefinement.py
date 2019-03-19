import numpy as np

class MeshOptions:
    """
    Mesh options for the domain
    """
    def __init__(self, shape):
        self.Shape = shape
        self.constraints = []

    def _addConstraint(self, entity, cons_type, index, variables):
        entities = ['vertex', 'segment', 'facet', 'region', 'global', 'point']
        assert entity in entities, \
            'wrong entity: '+str(entity)
        cons_types = ['fixed', 'around', 'function', 'TFI', 'box', 'boundary']
        assert cons_type in cons_types, \
            'wrong constraint type'
        assert isinstance(index, (list, tuple, int)) or index==None, \
            'must pass integer or list of index'
        assert isinstance(variables, (dict)), \
            'variables must be a dictionary'
        if isinstance(index, int): index = [index];
        self.constraints += [{'entity': entity, 'type': cons_type,
                              'index': index, 'variables':variables}]

    def refineSegment(self, ind, lc):
        """
        Refinement (or coarsening) of mesh inside segment.
        :param ind: list of local index of segment
        :param lc: size of element in segment
        """
        var_dict = {'Lc': lc}
        self._addConstraint(entity='segment', cons_type='fixed',
                            index=ind, variables=var_dict)
    def refineFacet(self, ind, lc):
        """
        Refinement (or coarsening) of mesh inside facet.
        :param ind: list of local index of facets
        :param lc: size of element in facet
        """
        var_dict = {'Lc': lc}
        self._addConstraint(entity='facet', cons_type='fixed',
                            index=ind, variables=var_dict)
    def refineRegion(self, ind, lc):
        """
        Refinement (or coarsening) of mesh inside region/volume.
        :param ind: list of local index of regions
        :param lc: size of element in region
        """
        var_dict = {'Lc': lc}
        self._addConstraint(entity='region', cons_type='fixed',
                            index=ind, variables=var_dict)

    def refineAroundVertex(self, ind, lc_min, lc_max=None, dist_min=None,
                            dist_max=None):
        """
        Refinement (or coarsening) of mesh around vertices.
        tetgen: only lc_min is taken into account
        gmsh: lc_min on segment, lc_max away from vertex, with a transition
              zone between dist_min and dist_max
        :param ind: list of local index of vertices
        :param lc_min: size of element at vertex
        :param lc_max: size of element away from vertex
        :param dist_min: distance away from vertex with lc_min element size
        :param dist_max: distance away from vertex with lc_max element size
        """
        var_dict = {'LcMin': lc_min, 'LcMax': lc_max,
                    'DistMin': dist_min, 'DistMax': dist_max}
        self._addConstraint(entity='vertex', cons_type='around',
                            index=ind, variables=var_dict)

    def refineAroundSegment(self, ind, lc_min, lc_max=None, dist_min=None,
                             dist_max=None):
        """
        Refinement (or coarsening) of mesh around segments.
        tetgen: only lc_min is taken into account
        gmsh: lc_min on segment, lc_max away from segment, with a transition
              zone between dist_min and dist_max
        :param ind: list of local index of segments
        :param lc_min: size of element at segment
        :param lc_max: size of element away from segment
        :param dist_min: distance away from segment with lc_min element size
        :param dist_max: distance away from segment with lc_max element size
        """
        var_dict = {'LcMin': lc_min, 'LcMax': lc_max,
                    'DistMin': dist_min, 'DistMax': dist_max,
                    'coords': coords}
        self._addConstraint(entity='segment', cons_type='around',
                            index=None, variables=var_dict)

    def refineAroundPoint(self, coords, lc_min, lc_max=None, dist_min=None,
                            dist_max=None):
        """
        Refinement (or coarsening) of mesh around vertices.
        tetgen: only lc_min is taken into account
        gmsh: lc_min on segment, lc_max away from vertex, with a transition
              zone between dist_min and dist_max
        :param ind: list of local index of vertices
        :param lc_min: size of element at vertex
        :param lc_max: size of element away from vertex
        :param dist_min: distance away from vertex with lc_min element size
        :param dist_max: distance away from vertex with lc_max element size
        """
        var_dict = {'LcMin': lc_min, 'LcMax': lc_max,
                    'DistMin': dist_min, 'DistMax': dist_max,
                    'coords': coords}
        self._addConstraint(entity='point', cons_type='around',
                            index=None, variables=var_dict)
    
    def refineAroundFacet(self, ind, lc_min, lc_max=None, dist_min=None,
                            dist_max=None):
        """
        Refinement (or coarsening) of mesh around facets.
        tetgen: only lc_min is taken into account
        gmsh: lc_min on segment, lc_max away from facet, with a transition
              zone between dist_min and dist_max
        (!) behaviour can be buggy in gmsh
        :param ind: list of local index of facets
        :param lc_min: size of element at facet
        :param lc_max: size of element away from facet
        :param dist_min: distance away from facet with lc_min element size
        :param dist_max: distance away from facet with lc_max element size
        """
        var_dict = {'LcMin': lc_min, 'LcMax': lc_max,
                    'DistMin': dist_min, 'DistMax': dist_max}
        self._addConstraint(entity='facet', cons_type='around',
                            index=ind, variables=var_dict)

    def refineBox(self, lc_in, lc_out, x_min, x_max, y_min, y_max, z_min=None,
                  z_max=None, restrict=None):
        """
        Refinement (or coarsening) of mesh inside a box.
        (!) for gmsh only
        :param lc_in: size of element inside box
        :param lc_out: size of element outside box
        :param x_min: lower limit of x coordinates of box
        :param x_max: upper limit of x coordinates of box
        :param y_min: lower limit of y coordinates of box
        :param y_max: upper limit of y coordinates of box
        :param z_min: lower limit of z coordinates of box
        :param z_max: upper limit of z coordinates of box
        """
        var_dict = {'VIn': lc_in, 'VOut': lc_out, 'XMin': x_min, 'XMax': x_max,
                    'YMin': y_min, 'YMax': y_max, 'ZMin': z_min, 'ZMax': z_max,
                    'restrict': restrict}
        self._addConstraint(entity='global', cons_type='box',
                            index=None, variables=var_dict)

    def setRefinementFunction(self, function, restrict=None):
        """
        Set a function to make the mesh element size vary.
        (!) for gmsh only. Must use MathEval syntax. Can use x, y and z in
            function
        :param function: function that makes mesh vary (string)
        """
        var_dict = {'function': function, 'restrict': restrict}
        self._addConstraint(entity='global', cons_type='function',
                            index=None, variables=var_dict)

    def setBoundaryLayerEdges(self, hwall_n, hwall_t, ratio=1.1, EdgesList=None,
                         newEdges=None, restrict=None):
        var_dict = {'hwall_n': hwall_n, 'hwall_t': hwall_t, 'ratio': ratio,
                    'newEdges': newEdges, 'restrict': restrict}
        self._addConstraint(entity='segment', cons_type='boundary', index=EdgesList,
                            variables=var_dict)

    def setTransfiniteSegment(self, ind, nb_nodes, prog=1.):
        """
        Sets segment transfinite interpolation. Goes from the segment first
        vertex to its second vertex using the defined progression.
        If TFI should go from segment's second vertex to first vertex,
        use a negative index.
        :param ind: list of local index of segments
        :param nb_nodes: number of nodes on segments
        :param prog: progression parameter for nodes on segments
        """
        var_dict_pos = {'nodes':nb_nodes, 'prog': prog}
        var_dict_neg = {'nodes':nb_nodes, 'prog': 1./prog}
        ind_prog_pos = []
        ind_prog_neg = []
        for ind in indice:
            if ind < 0:
                ind_prog_neg += [abs(ind)]
            else:
                ind_prog_pos += [ind]
        if ind_prog_neg:
            self._addConstraint(entity='segment', cons_type='TFI',
                                index=ind_prog_neg, variables=var_dict_neg)
        if ind_prog_pos:
            self._addConstraint(entity='segment', cons_type='TFI',
                                index=ind_prog_pos, variables=var_dict_pos)


# --------------------------------------------------------------------------- #

def _assembleRefinementOptions(domain):
    domain.MeshOptions.constraints = []
    for shape in domain.shape_list:
        cons = shape.MeshOptions.constraints
        for con in cons:
            domain.MeshOptions.constraints += [con]
            dcon = domain.MeshOptions.constraints[-1]
            if dcon['entity'] == 'vertex':
                dcon['index'] = (np.array(dcon['index'])+shape.start_vertex).tolist()
            if dcon['entity'] == 'segment':
                dcon['index'] = (np.array(dcon['index'])+shape.start_segment).tolist()
            if dcon['entity'] == 'facet':
                dcon['index'] = (np.array(dcon['index'])+shape.start_facet).tolist()
            if dcon['entity'] == 'region':
                dcon['index'] = (np.array(dcon['index'])+shape.start_region).tolist()


def writeGeo(domain, fileprefix, group_names=False, append=False):
    self = domain
    self.geofile = fileprefix
    self.polyfile = fileprefix
    geo = open(self.geofile+'.geo','w')
    pp = {}  # physical points
    pl = {}  # physical lines
    ps = {}  # physical surfaces
    pv = {}  # physical volume
    sN = len(self.segments)

    # Vertices
    geo.write('\n// Points\n')
    z = 0
    for i, v in enumerate(self.vertices):
        if self.nd == 3:
            z = v[2]
        geo.write("Point(%d) = {%g,%g,%g};\n" % (i+1,v[0],v[1], z))
        if self.vertexFlags:
            flag = self.vertexFlags[i]
            if flag in pp:
                pp[flag] += [i+1]
            else:
                pp[flag] = [i+1]
    nb_points = i+1

    lines_dict = {}
    for i in range(nb_points):
        lines_dict[i] = {}

    # Lines
    geo.write('\n// Lines\n')
    # line_list = []
    for i, s in enumerate(self.segments):
        geo.write("Line(%d) = {%d,%d};\n" % (i+1,s[0]+1,s[1]+1))
        # add segments in dictionary
        lines_dict[s[0]][s[1]] = i
        if self.segmentFlags:
            flag = self.segmentFlags[i]
            if flag in pl:
                pl[flag] += [i+1]
            else:
                pl[flag] = [i+1]
    nb_lines = i+1

    # Surfaces
    geo.write('\n// Surfaces\n')
    lines = 0
    lineloop_count = 0
    surface_line_list = []  # need to store for mesh constraints later
    # facet
    for i, f in enumerate(self.facets):
        seg_flag = sN+i+1
        lineloop = []
        lineloops = []
        lineloops_list = []
        line_list = []
        # subfacet
        if self.nd == 3 or (self.nd == 2 and i not in self.holes_ind):
            for j, subf in enumerate(f):
                lineloop = []
                # vertices in facet
                for k, ver in enumerate(subf):
                    if ver in lines_dict[subf[k-1]].keys():
                        lineloop += [lines_dict[subf[k-1]][ver]+1]
                    elif subf[k-1] in lines_dict[ver].keys():
                        # reversed
                        lineloop += [-(lines_dict[ver][subf[k-1]]+1)]
                    else:
                        ind = seg_flag+lines
                        lines += 1
                        geo.write('Line(%d) = {%d,%d};\n' % (ind, subf[k-1]+1, ver+1))
                        lineloop += [ind]
                line_list += lineloop
                geo.write('Line Loop(%d) = {%s};\n' % (lineloop_count+1, str(lineloop)[1:-1]))
                lineloops += [lineloop_count+1]
                lineloop_count += 1
            surface_line_list += [line_list]
            geo.write('Plane Surface(%d) = {%s};\n' % (i+1, str(lineloops)[1:-1]))
            if self.facetFlags:
                flag = self.facetFlags[i]
                if flag in ps:
                    ps[flag] += [i+1]
                else:
                    ps[flag] = [i+1]
        nb_lines += lines

    # Volumes
    geo.write('\n// Volumes\n')
    for i, V in enumerate(self.volumes):
        surface_loops = []
        if i not in self.holes_ind:
            for j, sV in enumerate(V):
                lineloop_count += 1
                geo.write('Surface Loop(%d) = {%s};\n' % (lineloop_count, str((np.array(sV)+1).tolist())[1:-1]))
                surface_loops += [lineloop_count]
            geo.write('Volume(%d) = {%s};\n' % (i+1, str(surface_loops)[1:-1]))
            flag = self.regionFlags[i]
            if self.regionFlags:
                if flag in pv:
                    pv[flag] += [i+1]
                else:
                    pv[flag] = [i+1]

    # Physical Groups
    geo.write('\n// Physical Groups\n')
    if self.boundaryTags:
        inv_bt = {v: k for k, v in self.boundaryTags.iteritems()}
    for flag in pp:
        ind = pp[flag]
        if self.boundaryTags and group_names is True:
            flag = '"'+inv_bt[flag]+'", '+str(flag)
        geo.write("Physical Point({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))
    for flag in pl:
        ind = pl[flag]
        if self.boundaryTags and group_names is True:
            flag = '"'+inv_bt[flag]+'", '+str(flag)
        geo.write("Physical Line({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))
    for flag in ps:
        ind = ps[flag]
        if self.boundaryTags and group_names is True:
            flag = '"'+inv_bt[flag]+'", '+str(flag)
        geo.write("Physical Surface({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))
    for flag in pv:
        ind = pv[flag]
        if self.boundaryTags and group_names is True:
            flag = '"'+inv_bt[flag]+'", '+str(flag)
        geo.write("Physical Volume({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))

    # Other
    mesh = self.MeshOptions
    geo.write('\n// ----------------\n')
    geo.write('\n// Other Operations\n')

    def idx2str(ind_list):
        return str((np.array(ind_list)+1).tolist())[1:-1] 

    def write_restrict_entity(entity, ind_list, nf):
        geo.write('Field[{0}] = Restrict; Field[{0}].IField = {1};\n'
                    .format(nf, nf-1))
        if entity == 'segment':
            write_restrict_segment(ind_list, nf)
        elif entity == 'facet':
            write_restrict_facet(ind_list, nf)
            # also refine segments of facets
        elif entity == 'region' and self.nd == 3:
            write_restrict_volume(ind_list, nf)     

    def write_restrict_segment(ind_list, nf):
        geo.write('Field[{0}].EdgesList = {{{1}}};\n'.format(nf, idx2str(ind_list)))

    def write_restrict_facet(ind_list, nf):
        geo.write('Field[{0}].FacesList = {{{1}}};\n'.format(nf, idx2str(ind_list)))
        sur_list = []
        for i in ind_list:
            sur_list += surface_line_list[i]
        geo.write('Field[{0}].EdgesList = {{{1}}};\n'
                    .format(nf, str(sur_list)[1:-1]))

    def write_restrict_volume(ind_list, nf):
        geo.write('Field[{0}].RegionsList = {{{1}}};\n'.format(nf, idx2str(ind_list)))
        for vol in ind_list:
            vol = np.array(self.volumes[i])
            i_list = [i for i in ind_list]
            faces_list = []
            edges_list = []
            for i in i_list:
                faces_list += [face+1 for subvol in self.volumes[i] for face in subvol]
            for i in faces_list:
                edges_list += surface_line_list[i]
            geo.write('Field[{0}].FacesList = {{{1}}};\n'
                        .format(nf, str(faces_list)[1:-1]))
            geo.write('Field[{0}].EdgesList = {{{1}}};\n'
                    .format(nf, str(edges_list)[1:-1]))

    def write_restrict(restrict_list, nf):
        geo.write('Field[{0}] = Restrict; Field[{0}].IField = {1};\n'
                    .format(nf, nf-1))
        for restrict in restrict_list:
            restrict_ent = restrict[0]
            restrict_ind = restrict[1] # index list
            if restrict_ent == 'segment':
                write_restrict_segment(restrict_ind, nf)
            elif restrict_ent == 'facet':
                write_restrict_facet(restrict_ind, nf)

    geo.write('\n// Fields\n')
    field_list = []
    nf = 1  # ID of next field to be defined
    for c in self.MeshOptions.constraints:
        if c['index']:
            ind = str((np.array(c['index'])+1).tolist())[1:-1] # index list
        v = c['variables'] # variables dictionary
        if c['type'] == 'fixed':
            # MathEval
            geo.write('Field[{0}] = MathEval; Field[{0}].F = "{1}";\n'
                        .format(nf, v['Lc']))
            nf += 1
            # Restrict field to entity
            write_restrict_entity(c['entity'], c['index'], nf)
            field_list += [nf]
            nf += 1
        if c['type'] == 'around':
            # Attractor
            geo.write("Field[{0}] = Attractor;\n".format(nf))
            if c['entity'] == 'vertex':
                geo.write("Field[{0}].NodesList = {{{1}}};\n".format(nf, ind))
            elif c['entity'] == 'segment':
                geo.write("Field[{0}].NNodesByEdge = 100;\n"
                            "Field[{0}].EdgesList = {{{1}}};\n"
                            .format(nf, ind))
            elif c['entity'] == 'facet':
                geo.write("Field[{0}].FacesList = {{{1}}};\n".format(nf, ind))
            elif c['entity'] == 'point':
                p = v['coords']
                if self.nd == 3:
                    z = p[2]
                elif self.nd == 2:
                    z = 0
                nb_points += 1
                geo.write("Point(%d) = {%g,%g,%g};\n" % (nb_points,p[0],p[1], z))
                geo.write("Field[{0}].NodesList = {{{1}}};\n".format(nf, nb_points))
            nf += 1
            # Threshold
            geo.write("Field[{0}] = Threshold; Field[{0}].IField = {1};\n"
                        "Field[{0}].LcMin = {2};\n"
                        .format(nf, nf-1, v['LcMin']))
            if v['LcMax']:
                geo.write("Field[{0}].LcMax = {1};\n".format(nf, v['LcMax']))
            if v['DistMin']:
                geo.write("Field[{0}].DistMin = {1};\n".format(nf, v['DistMin']))    
            if v['DistMax']:
                geo.write("Field[{0}].DistMax = {1};\n".format(nf, v['DistMax']))
            field_list += [nf]
            nf += 1
        elif c['type'] == 'TFI':
            if c['entity'] == 'segment':
                geo.write('Transfinite Line {{{0}}} = {1} Using Progression {2};\n'
                            .format(ind, v['nodes'], v['prog']))
        elif c['type'] == 'function':
            geo.write('Field[{0}] = MathEval;\n'
                        'Field[{0}].F = "{1}";\n'.format(nf, v['function']))
            if v['restrict'] is not None:
                nf += 1
                write_restrict(v['restrict'], nf)
            field_list += [nf]
            nf += 1
        elif c['type'] == 'box':
            geo.write('Field[{0}] = Box;\n'
                        'Field[{0}].VIn = {1}; Field[{0}].VOut = {2};\n'
                        'Field[{0}].XMin = {3}; Field[{0}].XMax = {4};\n'
                        'Field[{0}].YMin = {5}; Field[{0}].YMax = {6};\n'
                        .format(nf, v['VIn'], v['VOut'], v['XMin'], v['XMax'],
                                v['YMin'], v['YMax']))
            if self.nd == 3:
                geo.write('Field[{0}].ZMin = {1}; Field[{0}].ZMax = {2};\n'
                            .format(nf, v['ZMin'], v['ZMax']))
            if v['restrict'] is not None:
                nf += 1
                write_restrict(v['restrict'], nf)
            field_list += [nf]
            nf += 1
        elif c['type'] == 'boundary':
            edges =[]
            if v['newEdges'] is not None:
                for p in c['newEdges']:
                    if self.nd == 3:
                        z = p[2]
                    elif self.nd == 2:
                        z = 0
                    nb_points += 1
                    geo.write("Point(%d) = {%g,%g,%g};\n" % (nb_points,p[0],p[1], z))
                nb_lines += 1
                geo.write("Line(%d) = {%d, %d};\n" % (nb_lines, nb_points-1, nb_points))
                edges += [nb_lines]
            geo.write('Field[{0}] = BoundaryLayer;\n'
                      'Field[{0}].hwall_n = {1};\n'
                      'Field[{0}].ratio = {2};\n'
                      .format(nf, v['hwall_n'], v['ratio']))
            if c['index']:
                edges += [e+1 for e in c['index']]
            if edges:
                geo.write('Field[{0}].EdgesList = {{{1}}};\n'
                          .format(nf, str(edges)[1:-1]))
            field_list += [nf]
            nf += 1

    if append:
        pass
    else:
        geo.write('\n// Background Mesh\n')
        if nf == 1:
            # no other fields defined => constant background field
            geo.write(("Field[1] = MathEval; Field[1].F = '{0}';\n"
                        "Background Field = 1;\n".format(mesh.he)))
        else:
            geo.write("Field[{0}] = Min;\n"
                        "Field[{0}].FieldsList = {{{1}}};\n"
                        "Background Field = {0};\n"
                        .format(nf, str(field_list)[1:-1]))

    if self.MeshOptions.LcMax is not None:
        geo.write('Mesh.CharacteristicLengthMax = {0};\n'.format(self.MeshOptions.LcMax))

    # the following line does not work with refinement when 2 same entities
    # are defined..
    geo.write("Coherence;\n") # remove duplicates

    geo.close()


