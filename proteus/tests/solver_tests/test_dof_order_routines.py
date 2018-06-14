#!/usr/bin/env python
"""
Test module for routines that manage dof order
"""
from proteus import LinearSolvers as LS
import numpy as np
import pytest

@pytest.mark.parametrize("ownership_range,num_components",[
    ((5,25), 3),
    ((6,17), 4),
    ])
def test_interlaced_dof_order(ownership_range,
                              num_components):
    interlaced = LS.InterlacedDofOrderType()

    num_equations = ownership_range[1] - ownership_range[0] +1
    velocity, pressure = interlaced.create_DOF_lists(ownership_range,
                                                     num_equations,
                                                     num_components)
    idx_range = np.arange(ownership_range[0],ownership_range[1]+1)
    pressure_idx = np.arange(ownership_range[0],
                             ownership_range[1],
                             num_components)
    velocity_mask = np.ones(len(idx_range),dtype='bool')

    velocity_mask[pressure_idx -ownership_range[0]] = 0
    assert np.array_equal(pressure_idx, pressure)
    assert np.array_equal(idx_range[velocity_mask],
                          velocity)

@pytest.mark.parametrize("ownership_range,num_components,n_DOF_pressure",[
    ((5,25), 3, 7),
    ((6,17), 4, 3),
    ])
def test_blocked_dof_order(ownership_range,
                           num_components,
                           n_DOF_pressure):
    array_start = ownership_range[0]
    array_end = ownership_range[1]

    blocked = LS.BlockedDofOrderType(n_DOF_pressure)
    num_equations = array_end - array_start + 1

    velocity, pressure = blocked.create_DOF_lists(ownership_range,
                                                  num_equations,
                                                  num_components)
    pressure_idx = np.arange(array_start,
                             array_start + n_DOF_pressure)
    velocity_idx = np.arange(pressure_idx[-1] + 1,
                             array_end + 1)
    assert np.array_equal(pressure_idx, pressure)
    assert np.array_equal(velocity_idx, velocity)
