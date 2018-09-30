import pytest
import os
import sys
import inspect
import numpy
import tables
import pickle

import proteus.test_utils.TestTools as TestTools
from proteus.iproteus import *

TestTools.addSubFolders( inspect.currentframe() )
import cavity2d
import twp_navier_stokes_cavity_2d_so
import twp_navier_stokes_cavity_2d_p
import twp_navier_stokes_cavity_2d_n

@pytest.fixture(scope='module')
def load_cavity_problem():
    reload(cavity2d)
    reload(twp_navier_stokes_cavity_2d_so)
    reload(twp_navier_stokes_cavity_2d_p)
    reload(twp_navier_stokes_cavity_2d_n)
    pList = [twp_navier_stokes_cavity_2d_p]
    nList = [twp_navier_stokes_cavity_2d_n]
    so = twp_navier_stokes_cavity_2d_so
    so.sList = [default_s]
    yield pList, nList, so
