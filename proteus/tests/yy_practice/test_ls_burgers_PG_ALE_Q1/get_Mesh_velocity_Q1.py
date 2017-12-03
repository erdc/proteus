from proteus import *
import math
import numpy as np
import cmath


smoothing_times = 0
min_area_in_previous_step = 1.0
min_angle_in_previous_step = math.pi * 0.25
max_angle_in_previous_step = math.pi * 0.5


def get_angle(x1, x2, x3):
    z1 = complex(x1[0] - x2[0], x1[1] - x2[1])
    z2 = complex(x3[0] - x2[0], x3[1] - x2[1])
    try:
        z2 /= z1
    except ZeroDivisionError:
        import pdb
        pdb.set_trace()
        print z1, z2, x1, x2, x3

    # abs since we do not know the orientation
    return abs(cmath.phase(z2))

def get_triangle_area(e_nodes):
    detJ = (e_nodes[1][0] - e_nodes[0][0]) * (e_nodes[2][1] -
                                              e_nodes[0][1]) - (e_nodes[1][1] - e_nodes[0][1]) * (e_nodes[2][0] -
                                                                                                  e_nodes[0][0])
    area = 0.5*abs(detJ)
    
    return area
    
def get_angle_area(e_nodes):
    
    # one 
    area = get_triangle_area(e_nodes[[0,1,2],:])+get_triangle_area(e_nodes[[0,2,3],:])
        
    # use for-loop starting from -1
    a1 = get_angle(e_nodes[3], e_nodes[0], e_nodes[1])
    a2 = get_angle(e_nodes[0], e_nodes[1], e_nodes[2])
    a3 = get_angle(e_nodes[1], e_nodes[2], e_nodes[3])
    a4 = get_angle(e_nodes[2], e_nodes[3], e_nodes[0])
    
    min_angle = min([a1, a2, a3, a4])
    max_angle = max([a1, a2, a3, a4])

    return min_angle, max_angle, area


def get_mesh_info(node_coord, element_nodes):

    min_angle, max_angle, min_area = math.pi, 0, 1000

    for e in xrange(element_nodes.shape[0]):
        e_min_angle, e_max_angle, e_area = get_angle_area(
            node_coord[element_nodes[e]])
        # print e, e_min_angle, e_max_angle, e_area

        if e_min_angle < 0.0 or e_area < 0.0 or e_max_angle < 0.0:
            print "element is ", e
            print node_coord[element_nodes[e]]
            return e_min_angle, e_max_angle, e_area

        min_angle = min([e_min_angle, min_angle])
        max_angle = max([e_max_angle, max_angle])
        min_area = min([e_area, min_area])

    return min_angle, max_angle, min_area


def get_smoothed_coord(_old_coord, node_star_offset,
                       node_star_array, fixed_node_index,
                       N, _new_coord):

    _new_coord[:] = 0.0

    for _ in range(N):
        for i in xrange(_old_coord.shape[0]):
            if i not in fixed_node_index:
                n_neightbor = 0
                star_center = np.zeros((3,), 'd')
                for j in range(node_star_offset[i], node_star_offset[i + 1]):
                    if node_star_array[j] != i:
                        star_center[:] += _old_coord[node_star_array[j]]
                        n_neightbor += 1

                _new_coord[i] = star_center / n_neightbor
            else:
                _new_coord[i] = _old_coord[i]
        # swap two object
        _new_coord, _old_coord = _old_coord, _new_coord


def get_mesh_velocity(node_coord,
                      node_star_offset,
                      node_star_array,
                      element_nodes,
                      fixed_node_index,
                      dt,
                      initial_area,
                      _mesh_velocity_at_node):
    global smoothing_times, min_area_in_previous_step, max_angle_in_previous_step, min_angle_in_previous_step
    # Lagrangian mesh
    _node_coord_lagrange = np.copy(node_coord)
    _node_coord_lagrange += dt * _mesh_velocity_at_node


    _node_coord_smoothed_before = np.copy(_node_coord_lagrange)
    _node_coord_smoothed_after = np.copy(_node_coord_lagrange)

    # mesh info
    _min_angle, _max_angle, _min_area = get_mesh_info(
        _node_coord_smoothed_before, element_nodes)

    print _min_angle, _max_angle, _min_area

    if smoothing_times == 0:
        min_area_in_previous_step = initial_area

    # use smoothed coord to modify lang coord
    while True:
        if _min_area < 0.0 or _min_angle < 0.0:
            # If N=0, since Lagrangian mesh has more weight, maybe it has dead
            # loop.
            N = 8
            # assume previous mesh is regular
            _node_coord_smoothed_after[:] = node_coord
        elif _min_area < 0.8 * min_area_in_previous_step or _max_angle > 1.2 * max_angle_in_previous_step or _min_angle < 0.8 * min_area_in_previous_step or math.pi - 0.3 < _max_angle or _min_angle < 0.3:
            N = 4
        else:
            N = 2

        # continue to smooth previous smoothed mesh
        _node_coord_smoothed_before[:] = _node_coord_smoothed_after
        get_smoothed_coord(_node_coord_smoothed_before, node_star_offset,
                           node_star_array, fixed_node_index,
                           N,
                           _node_coord_smoothed_after)

        smoothing_times += N
        print '>>' * 10, 'smoothing mesh...', smoothing_times

        # Use big weight on Lagrangian mesh since the problem is well-posed.
        # Otherwise, the mesh velocity should be regular enough.
        _node_coord_smoothed_before[:] = _node_coord_lagrange
        _node_coord_smoothed_before *= 0.5
        _node_coord_smoothed_before += 0.5 * _node_coord_smoothed_after

        _min_angle, _max_angle, _min_area = get_mesh_info(
            _node_coord_smoothed_before, element_nodes)

        # Give good check rule to avoid dead loop
        if _min_area > 0.4 * min_area_in_previous_step and math.pi - 0.2 > _max_angle or _min_angle > 0.2:
            break

    # save for the next step
    min_angle_in_previous_step, max_angle_in_previous_step, min_area_in_previous_step = _min_angle, _max_angle, _min_area

    # if no smoothing happened, it is equal to Lagrangian coordinate mesh.
    #_node_coord_lagrange *= 0.2
    #_node_coord_lagrange += 0.8 * _node_coord_smoothed_after
    # use final mesh to get velocity
    #_node_coord_lagrange -= node_coord
    #_node_coord_lagrange /= dt
    #_mesh_velocity_at_node[:] = _node_coord_lagrange

    _mesh_velocity_at_node[:] = (_node_coord_smoothed_before - node_coord) / dt

    #_mesh_velocity_at_node[:] = 0.0
