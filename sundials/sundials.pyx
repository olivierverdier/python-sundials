# -*- coding: utf-8 -*-
from libc.stdlib cimport malloc, free
from SundialsLib cimport *

class SundialsError(Exception):
    pass

cdef class SolverStatistics:
    """Class to encapsulate statistics for solver"""
    def __repr__(self):
        s = []
        s.append("nsteps = %d"      % self.nsteps)
        s.append("nfevals = %d"     % self.nfevals)
        s.append("njevals = %d"     % self.njevals)
        s.append("ngevals = %d"     % self.ngevals)
        s.append("netfails = %d"    % self.netfails)
        s.append("nniters = %d"     % self.nniters)
        s.append("nncfails = %d"    % self.nncfails)
    
        return "(%s)" % (", ".join(s))
    
    def __str__(self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())

include "NVECTOR.pxi"
include "CVODE.pxi"
include "IDA.pxi"
include "KINSOL.pxi"
