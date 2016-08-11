import numpy as np
from proteus.MeshTools import Edge
from proteus.Profiling import logEvent

def msh2triangle(fileprefix):
    mshfile = open(fileprefix+'.msh', 'r')
    nodes = []
    edges_msh = []
    triangles = []
    tetrahedra = []
    triangle_nb = 0
    edge_nb = 0
    nd = 2
    switch = None
    switch_count = -1
    logEvent('msh2triangle: getting nodes and elements')
    for line in mshfile:
        if 'Nodes' in line:
            switch = 'nodes'
            switch_count = -1
        if 'Elements' in line:
            switch = 'elements'
            switch_count = -1
        if switch == 'nodes' and switch_count >= 0:
            words = line.split()
            if switch_count == 0:
                node_nb = int(words[0])
            else:
                nid = int(words[0])
                if nd == 2:
                    x, y, z = float(words[1]), float(words[2]), 0
                elif nd == 3:
                    x, y, z = float(words[1]), float(words[2]), float(words[3])
                nodes += [[nid, x, y, z, 0]]
        if switch == 'elements' and switch_count >= 0:
            words = line.split()
            if switch_count == 0:
                el_nb = int(words[0])
            else:
                el_id = int(words[0])
                el_type = int(words[1])
                nb_tags = int(words[2])
                if nb_tags == 2:
                    flag = int(words[3])
                else:
                    flag = 0
                s = 3+nb_tags # starting index on words for element info
                if el_type == 1: # segment
                    edge_nb += 1
                    edges_msh += [[edge_nb, int(words[s]), int(words[s+1]), flag]]
                elif el_type == 2: # triangle
                    triangle_nb += 1
                    triangles += [[triangle_nb, int(words[s]), int(words[s+1]), int(words[s+2]), flag]]
                elif el_type == 15: # node
                    nodes[el_id-1][4] = flag
        switch_count += 1
    mshfile.close()

    # construct ALL edges with flags and add flags to nodes
    edges_dict = {}
    triangles = np.array(triangles)
    edge_nb = 0
    edges = []

    logEvent('msh2triangle: constructing edges')
    for triangle in triangles[:,1:4]:  # take only vertices index
        for i in range(len(triangle)):
            edge = Edge(edgeNumber=edge_nb, nodes=[triangle[i-1], triangle[i]])
            edge_exist = bool(edges_dict.get(edge.nodes))
            if not edge_exist:
                edge_nb += 1
                edges_dict[edge.nodes] = edge
                edges += [[edge_nb, edge.nodes[0], edge.nodes[1], 0]]
    logEvent('msh2triangle: updating edges and nodes flags')
    edges = np.array(edges)
    for edge in edges_msh:
        edge_nodes = [edge[1], edge[2]]
        edge_nodes.sort()
        edge_nodes = tuple(edge_nodes)
        edge_class = edges_dict.get(edge_nodes)
        edges[edge_class.N, 3] = edge[3]
        # ! edge nodes are indexed from 1 with gmsh
        if nodes[edge[1]-1][-1] == 0:  # update node flags
            nodes[edge[1]-1][-1] = edge[3]
        if nodes[edge[2]-1][-1] == 0:  # update node flags
            nodes[edge[1]-1][-1] = edge[3]


    logEvent('msh2triangle: writing .node .ele .edge files')
    header = '{0:d} {1:d} 0 1'.format(node_nb, nd)

    if nd == 2:
        nodes = np.array(nodes)
        nodes = np.delete(nodes, 3, 1)
        fmt = ['%d', '%f', '%f', '%d']
    elif nd == 3:
       fmt = ['%d', '%f', '%f', '%f', '%d']
    np.savetxt(fileprefix+'.node', nodes, fmt=fmt, header=header, comments='')

    header = '{0:d} 1'.format(edge_nb)
    np.savetxt(fileprefix+'.edge', edges, fmt='%d', header=header, comments='')

    if nd == 2:
        header = '{0:d} 3 1'.format(triangle_nb)
        np.savetxt(fileprefix+'.ele', triangles, fmt='%d', header=header, comments='')

    logEvent('msh2triangle: finished converting .msh to triangle files')
