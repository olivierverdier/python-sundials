#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
Cython Wrapper for interfacing Python with CVode and IDA (Sundials Version 2.4.0)
Claus Fuhrer,        Lund University        October 2009
Christian Andersson, Lund University        Februari 2010

see also Jon Olav Vik: 
http://codespeak.net/pipermail/cython-dev/2009-June/005947.html

"""

# ==============================================
#  external definitions from Sundial headers
# ==============================================

cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype
    ctypedef bint booleantype # should be bool instead of bint, but there is a bug in Cython
    # This bug is fixed in http://trac.cython.org/cython_trac/ticket/227

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void* content
    ctypedef _generic_N_Vector* N_Vector
    N_Vector N_VNew_Serial(long int vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)
    
    
cdef extern from "nvector/nvector_serial.h":
    cdef struct _N_VectorContent_Serial:
        long int length
        realtype* data
    ctypedef _N_VectorContent_Serial* N_VectorContent_Serial
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data)

cdef extern from "sundials/sundials_direct.h":
    cdef struct _DlsMat:
        int type
        int M
        int N
        int ldim
        int mu
        int ml
        int s_mu
        realtype *data
        int ldata
        realtype **cols
    ctypedef _DlsMat* DlsMat
    cdef realtype* DENSE_COL(DlsMat A, int j)

cdef extern from "cvode/cvode.h":
    void* CVodeCreate(int lmm, int iter)
    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *f_data)
    int CVodeInit(void* cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
    int CVodeReInit(void* cvode_mem, realtype t0, N_Vector y0)
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    int CVodeSetStopTime(void* cvode_mem, realtype tstop)
    int CVodeGetIntegratorStats(void* cvode_mem, long int *nsteps, long int *nfevals,
                                long int *nlinsetups, long int *netfails, int *qlast, int *qcur,
                                realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)
    int CVodeSetMaxOrd(void * cvode_mem, int maxord)
    int CVodeSetMaxNumSteps(void * cvode_mem, long int mxsteps)
    int CVodeSetMaxStep(void* cvode_mem, realtype hmax)
    int CVodeSetInitStep(void * cvode_mem, realtype hin)
    
    # Error in def?
    void CVodeFree(void *cvode_mem)
    
    int CVodeStep "CVode"(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, 
        int itask) nogil
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky) nogil
    
    int CVodeSetUserData(void *cvode_mem,void *user_data)
    
    # functions for discontinuity handling
    ctypedef int (*CVRootFn)(realtype tt, N_Vector yy, realtype *gout, void *user_data)
    int CVodeRootDirection(void *cvode_mem, int *rootdir)
    int CVodeSetNoInactiveRootWarn(void *cvode_mem)
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
    
    #Functions for retrieving statistics
    int CVodeGetLastOrder(void * cvode_mem,int *qlast)
    int CVodeGetCurrentOrder(void * cvode_mem,int *qcurrent)
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps) #Number of steps
    int CVodeGetLastStep(void *cvode_mem, realtype *hlast) # Last internal step
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nrevals) #Number of function evals
    int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals) #Number of jac evals
    int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nrevalsLS) #Number of res evals due to jac evals
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals) #Number of root evals
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails) #Number of local error test failures
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters) #Number of nonlinear iteration
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails) #Number of nonlinear conv failures
    
cdef extern from "cvode/cvode_dense.h":
    int CVDense(void *cvode_mem, long int N)
    ctypedef int (*CVDlsDenseJacFn)(int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn djac)
    
cdef extern from "ida/ida.h":
    ctypedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    void* IDACreate()
    int IDAInit(void* ida_mem, IDAResFn res, realtype t0, N_Vector y0, N_Vector yp0)
    int IDAReInit(void* ida_mem, realtype t0, N_Vector y0, N_Vector yp0)
    int IDASetStopTime(void* ida_mem, realtype tstop)
    int IDASetMaxNumSteps(void * cvode_mem, long int mxsteps)
    int IDASetMaxOrd(void * cvode_mem, int maxord)
    int IDASetMaxStep(void* ida_mem, realtype hmax)
    void IDAFree(void *cvode_mem)
    int IDAGetIntegratorStats(void* ida_mem,long int  *nsteps, long int *nrevals, 
                            long int *nlinsetups, long int *netfails, int *klast, 
                            int *kcur, realtype *hinused, realtype *hlast, 
                            realtype *hcur, realtype *tcur)
    int IDASolve(void* ida_mem, realtype tout,realtype  *tret, N_Vector yret, 
                            N_Vector ypret, int itask)
    int IDAGetSolution(void *ida_mem, realtype t, N_Vector yret, N_Vector ypret)
    int IDASetUserData(void *ida_mem,void *user_data)
    # functions to control the error test
    int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol)
    int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol)
    int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
    int IDASetId(void *ida_mem, N_Vector id)
    
    # functions for discontinuity handling
    ctypedef int (*IDARootFn)(realtype tt, N_Vector yy, N_Vector yp, realtype *gout, void *user_data)
    int IDASetRootDirection(void *ida_mem, int *rootdir)
    int IDASetNoInactiveRootWarn(void *ida_mem)
    int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g)
    int IDAGetRootInfo(void *ida_mem, int *rootsfound)
    int IDACalcIC(void *ida_men, int icopt, realtype tout1)
    int IDAGetConsistentIC(void *ida_mem, N_Vector y0, N_Vector yp0)
    int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
    
    #Functions for retrieving statistics
    int IDAGetLastOrder(void *ida_mem,int *qlast) #Last order used
    int IDAGetCurrentOrder(void *ida_mem,int *qcurrent) #Order that is about to be tried
    int IDAGetNumSteps(void *ida_mem, long int *nsteps) #Number of steps
    int IDAGetNumResEvals(void *ida_mem, long int *nrevals) #Number of res evals
    int IDAGetLastStep(void *ida_mem, realtype *hlast) # last step size
    int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals) #Number of jac evals
    int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS) #Number of res evals due to jac evals
    int IDAGetNumGEvals(void *ida_mem, long int *ngevals) #Number of root evals
    int IDAGetNumErrTestFails(void *ida_mem, long int *netfails) #Number of local error test failures
    int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters) #Number of nonlinear iteration
    int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails) #Number of nonlinear conv failures

cdef extern from "ida/ida_dense.h":
    int IDADense(void *ida_mem, long int N)
    ctypedef int (*IDADlsDenseJacFn)(int Neq, realtype tt, realtype cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn djac)

cdef extern from "kinsol/kinsol_direct.h":
    int KINDense(void *kinmem, int N)
    
cdef extern from "kinsol/kinsol.h":
    void *KINCreate()
    void KINFree(void *kinmem)
    
    ctypedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void *user_data)
    int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl)
    
    int KINSetUserData(void *kinmem, void *user_data)
    int KINSetConstraints(void *kinmem, N_Vector constraints)
    int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
    int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
    int KINSetNumMaxIters(void *kinmem, long int mxiter)
    int KINSetMaxSetupCalls(void *kinmem, long int msbset)
    
    int KINSol(void *kinmem, N_Vector uu, int strategy, N_Vector u_scale, N_Vector f_scale)
    
    int KINGetNumFuncEvals(void *kinmem, long int *nfevals)
    int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters)
    
#===============================================================
# Constants
#===============================================================
cdef enum:
    #a) CVODE in
    CV_RHS_IND        = 0   # Index to user data rhs handling
    CV_RHSF_IND       = 0   # Index to user data rhs
    CV_JAC_IND        = 1   # Index to user data jacobian
    CV_ROOT_IND       = 1   # Index to user data root handling
    CV_ROOTF_IND      = 0   # Index to user data root function
    CV_SW_IND         = 1   # Index to user data root switches
    CV_ADAMS = 1
    CV_BDF   = 2
    CV_FUNCTIONAL = 1
    CV_NEWTON     = 2
    CV_SS = 1
    CV_SV = 2
    CV_WF = 3
    CV_NORMAL         = 1
    CV_ONE_STEP       = 2
    CV_NORMAL_TSTOP   = 3
    CV_ONE_STEP_TSTOP = 4

    #b) CVODE out
    CV_SUCCESS = 0
    CV_TSTOP_RETURN = 1
    CV_ROOT_RETURN    = 2   # CVSolve succeeded and found one or more roots.
    
    CV_TOO_MUCH_WORK        = -1  # The solver took mxstep internal steps but still could not reach tout
    CV_TOO_MUCH_ACC         = -2  # The solver could not satisfy the accuracy
    CV_ERR_FAILURE          = -3  # Error test failures occurred too many times
    CV_CONV_FAILURE         = -4  # Convergence test failures occurred too many times
    CV_LINIT_FAIL           = -5  # The linear solver's initialization function failed.
    CV_LSETUP_FAIL          = -6  # The linear solver's setup function failed in an unrecoverable manner.
    CV_LSOLVE_FAIL          = -7  # The linear solver's solve function failed in an unrecoverable manner.
    CV_RHSFUNC_FAIL         = -8  # The right-hand side function failed in an unrecoverable manner.
    CV_FIRST_RHSFUNC_ERR    = -9  # The right-hand side function had a recoverable error at the first call.
    CV_REPTD_RHSFUNC_ERR    = -10 # Convergence test failures occurred too many times
    CV_UNREC_RHSFUNC_ERR    = -11 # The right-hand function had a recoverable error, but no recovery was possible
    CV_RTFUNC_FAIL          = -12 # The rootfinding function failed in an unrecoverable manner.
    
    CV_MEM_FAIL       = -20 # Memory failure
    CV_MEM_NULL       = -21 # The cvode_mem argument was NULL.
    CV_ILL_INPUT      = -22 # Inputs to CVode was illegal
    CV_NO_MALLOC      = -23 # The CVODE memory was not allocated by a call to CVodeInit.
    CV_BAD_K          = -24 # k is not in the range 0, 1,..., qu.
    CV_BAD_T          = -25 # t is not in the interval [tn - hu, tn].
    CV_BAD_DKY        = -26 # The dky argument was NULL.
    CV_TOO_CLOSE      = -27 # The initial time t0 and the final time tout are too close
    
    #c) IDA in
    IDA_NORMAL         = 1   # Solver returns at specified output time.
    IDA_ONE_STEP       = 2   # Solver returns after each successful step.
    IDA_RES_IND        = 0   # Index to user data residual handling
    IDA_RESF_IND       = 0   # Index to user data residual
    IDA_JAC_IND        = 1   # Index to user data jacobian
    IDA_ROOT_IND       = 1   # Index to user data root handling
    IDA_ROOTF_IND      = 0   # Index to user data root function
    IDA_SW_IND         = 1   # Index to user data root switches
    IDA_YA_YDP_INIT    = 1   # See IDA Documentation 4.5.4
    IDA_Y_INIT         = 2   # See IDA Documentation 4.5.4

    #d) IDA out
    IDA_SUCCESS        = 0   # Successful function return.   
    IDA_TSTOP_RETURN   = 1   # IDASolve succeeded by reaching the specified stopping point.
    IDA_ROOT_RETURN    = 2   # IDASolve succeeded and found one or more roots.
    
    IDA_MEM_NULL        = -1    # The ida_mem argument was NULL
    IDA_ILL_INPUT       = -2    # Inputs to IDA was illegal
    IDA_NO_MALLOC       = -3    # The allocation function IDAInit has not been called.
    IDA_TOO_MUCH_WORK   = -4    # The solver took mxstep internal steps but could not reach tout.
    IDA_TOO_MUCH_ACC    = -5    # The solver could not satisfy the accuracy
    IDA_ERR_FAIL        = -6    # Error test failures occurred too many times
    IDA_CONV_FAIL       = -7    # Convergence test failures occurred too many times
    IDA_LINIT_FAIL      = -8    # The linear solver's initialization function failed
    IDA_LSETUP_FAIL     = -9    # The linear solver's setup function failed in an unrecoverable manner
    IDA_LSOLVE_FAIL     = -10   # The linear solver's solve function failed in an unrecoverable manner
    IDA_RES_FAIL        = -11   # The user's residual function returned a nonrecoverable error flag.
    IDA_CONSTR_FAIL     = -12   # The inequality constraints were violated and the solver was unable to recover.
    IDA_REP_RES_ERR     = -13   # Unable recover from multiple residual recoverable error flags
    IDA_MEM_FAIL        = -14   # A memory allocation request has failed.
    IDA_BAD_T           = -15   # t is not in the interval [tn - hu, tn].
    IDA_BAD_EWT         = -16   # Some component of the error weight vector is zero
    IDA_FIRST_RES_FAIL  = -17   # Unable recover from first residual recoverable error flag
    IDA_LINESEARCH_FAIL = -18   # The linesearch algorithm failed to find a solution
    IDA_NO_RECOVERY     = -19   # Unable recover from residual recoverable error flag
    IDA_RTFUNC_FAIL     = -20   # The rootfinding function failed.
    
    KIN_NONE                    = 0
    KIN_LINESEARCH              = 1

    KIN_SUCCESS                 = 0
    KIN_INITIAL_GUESS_OK        = 1
    KIN_STEP_LT_STPTOL          = 2

    KIN_WARNING                 = 99

    KIN_MEM_NULL                = -1
    KIN_ILL_INPUT               = -2
    KIN_NO_MALLOC               = -3
    KIN_MEM_FAIL                = -4
    KIN_LINESEARCH_NONCONV      = -5
    KIN_MAXITER_REACHED         = -6
    KIN_MXNEWT_5X_EXCEEDED      = -7
    KIN_LINESEARCH_BCFAIL       = -8
    KIN_LINSOLV_NO_RECOVERY     = -9
    KIN_LINIT_FAIL              = -10
    KIN_LSETUP_FAIL             = -11
    KIN_LSOLVE_FAIL             = -12

    KIN_SYSFUNC_FAIL            = -13
    KIN_FIRST_SYSFUNC_ERR       = -14
    KIN_REPTD_SYSFUNC_ERR       = -15

    KINDLS_SUCCESS              =  0
    KINDLS_MEM_NULL             = -1
    KINDLS_LMEM_NULL            = -2
    KINDLS_ILL_INPUT            = -3
    KINDLS_MEM_FAIL             = -4
    KINDLS_JACFUNC_UNRECVR      = -5
    KINDLS_JACFUNC_RECVR        = -6
