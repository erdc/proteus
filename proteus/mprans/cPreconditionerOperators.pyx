cdef extern from "PreconditionerOperators.h" namespace "proteus":
    cdef cppclass RANS2P_Op_Builder:
        RANS2P_Op_Builder() except +
        void attachTwoPhaseInvScaledMassOperator()

cdef class cRANS2P_Op_Builder:
    cdef RANS2P_Op_Builder op_instance
    def __cinit__(self):
        self.op_instance = RANS2P_Op_Builder()
    def attachTwoPhaseInvScaledMassOperator(self):
        self.op_instance.attachTwoPhaseInvScaledMassOperator()