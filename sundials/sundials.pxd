# -*- coding: utf-8 -*-
ctypedef double realtype

cdef class N_Vector_Serial:
    cdef void *thisptr
    cdef Py_ssize_t __shape[1]
    
    cdef int getsize(N_Vector_Serial self)
    cdef double *getdata(N_Vector_Serial self)
    cdef double getitem(N_Vector_Serial self, int index)
    cdef setitem(N_Vector_Serial self, int index, double value)
    cpdef copy(N_Vector_Serial self)

cdef double *N_Vector_data(void *vector)
cdef int N_Vector_length(void *vector)

cdef struct cv_userdata:
    void *RHSF
    void *ROOTF
    void *JACF
    void *SW
    void *Y
    
cdef class SolverStatistics:
    cdef readonly long int nsteps
    cdef readonly long int nfevals
    cdef readonly long int njevals
    cdef readonly long int ngevals
    cdef readonly long int netfails
    cdef readonly long int nniters
    cdef readonly long int nncfails

cdef class CVodeSettings:
    cdef public object lmm
    cdef public object iter
    cdef public int maxord
    cdef public long int mxsteps
    cdef public realtype tstop
    cdef public realtype hmax
    cdef public realtype reltol
    cdef public realtype abstol
    cdef public object abstolv
    
    cdef public object RHS
    cdef public object ROOT
    cdef public object SW
    cdef public object JAC

cdef class CVodeSolver

cdef class CVodeIterator:
    cdef CVodeSolver solver
    
    cdef public realtype dt
    cdef readonly realtype t
    
    cdef realtype tret
    cdef realtype troot
    cdef N_Vector_Serial yroot
    cdef realtype hlast
    
    cdef bint root
    cdef bint stop
    cdef bint last
    
    cdef N_Vector_Serial Next(self)

cdef class CVodeSolver:
    cdef void *thisptr
    cdef void *RhsFn
    cdef void *RootFn
    cdef void *JacFn
    
    cdef realtype t0
    cdef N_Vector_Serial y
    cdef N_Vector_Serial yv
    cdef cv_userdata userdata
    cdef public CVodeSettings settings
    
    cdef object RHSF
    cdef object ROOTF
    cdef object JACF
    cdef public object SW
    
    cpdef __handleRoot(self, realtype event_time, N_Vector_Serial y)
    cpdef CVodeIterator iter(self, realtype t0, realtype dt)
    cpdef N_Vector_Serial step(self, realtype tf)

cdef struct ida_userdata:
    void *RESF
    void *ROOTF
    void *JACF
    void *SW
    void *Y
    void *YDOT

cdef class IDASettings:
    cdef public int maxord
    cdef public long int mxsteps
    cdef public realtype tstop
    cdef public realtype hmax
    cdef public bint suppressalg
    cdef public bint lsoff
    cdef public icopt
    cdef public object idv
    cdef public realtype reltol
    cdef public realtype abstol
    cdef public object abstolv
    
    cdef public object RES
    cdef public object ROOT
    cdef public object SW
    cdef public object JAC

cdef class IDASolver

cdef class IDAIterator:
    cdef IDASolver solver
    
    cdef public realtype dt
    cdef readonly realtype t
    
    cdef realtype tret
    cdef realtype troot
    cdef N_Vector_Serial yroot
    cdef N_Vector_Serial ydotroot
    cdef realtype hlast
    
    cdef bint root
    cdef bint stop
    cdef bint last
    
    cdef Next(self)
    
cdef class IDASolver:
    cdef void *thisptr
    cdef void *ResFn
    cdef void *RootFn
    cdef void *JacFn
    
    cdef realtype t0
    
    cdef N_Vector_Serial y
    cdef N_Vector_Serial yv
    cdef N_Vector_Serial ydot
    cdef N_Vector_Serial ydotv
    
    cdef object RESF
    cdef object ROOTF
    cdef object JACF
    cdef public object SW
    
    cdef ida_userdata userdata
    cdef public IDASettings settings
    
    cpdef GetConsistentIC(self, realtype direction)
    cpdef __handleRoot(self, realtype event_time, N_Vector_Serial yroot, N_Vector_Serial ydotroot)
    cpdef IDAIterator iter(self, realtype t0, realtype dt)
    cpdef step(self, realtype tf)

cdef struct kinsol_userdata:
    void *F
    void *U

cdef class KINSOLStatistics:
    cdef readonly long int nfevals
    cdef readonly long int nniters
    
cdef class KINSOLSettings:
    cdef readonly int NVAR
    cdef readonly int NEQ
    cdef public int printfl
    cdef public object strategy
    cdef public long int mxiter
    cdef public long int msbset
    cdef public double fnormtol
    cdef public double scsteptol
    
    cdef public object F
    cdef public N_Vector_Serial x0
    cdef public N_Vector_Serial cstr
    
cdef class KINSOLSolver:
    cdef void *thisptr
    cdef kinsol_userdata userdata
    cdef N_Vector_Serial u
    cdef N_Vector_Serial u_scale
    cdef N_Vector_Serial f_scale
    cdef public KINSOLSettings settings
    cdef public N_Vector_Serial x