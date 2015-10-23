#!/usr/bin/env python


def read_from_hdf5(hdfFile,label,dof_map=None,has_h5py=True):
    """
    Just grab the array stored in the node with label label and return it
    If dof_map is not none, use this to map values in the array
    If dof_map is not none, this determines shape of the output array
    """
    assert hdfFile != None, "requires hdf5 for heavy data"
    if not has_h5py:
        vals = hdfFile.getNode(label).read()
    else:
        vals = hdfFile[label][:]
    if dof_map is not None:
        dof = vals[dof_map]
    else:
        dof = vals

    return dof
