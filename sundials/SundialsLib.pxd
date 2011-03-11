# -*- coding: utf-8 -*-
cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype
    ctypedef bint booleantype

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void *content
        
    ctypedef _generic_N_Vector *N_Vector
    N_Vector N_VNew_Serial(long int vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)

cdef extern from "nvector/nvector_serial.h":
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data)
    
    cdef struct _N_VectorContent_Serial:
        long int length
        realtype *data
        
    ctypedef _N_VectorContent_Serial *N_VectorContent_Serial
    
cdef extern from "cvode/cvode.h":
    int CV_ADAMS
    int CV_BDF
    int CV_FUNCTIONAL
    int CV_NEWTON
    int CV_NORMAL
    int CV_ONE_STEP

    int CV_SUCCESS
    int CV_TSTOP_RETURN
    int CV_ROOT_RETURN

    int CV_WARNING

    int CV_TOO_MUCH_WORK
    int CV_TOO_MUCH_ACC
    int CV_ERR_FAILURE
    int CV_CONV_FAILURE

    int CV_LINIT_FAIL
    int CV_LSETUP_FAIL
    int CV_LSOLVE_FAIL
    int CV_RHSFUNC_FAIL
    int CV_FIRST_RHSFUNC_ERR
    int CV_REPTD_RHSFUNC_ERR
    int CV_UNREC_RHSFUNC_ERR
    int CV_RTFUNC_FAIL

    int CV_MEM_FAIL
    int CV_MEM_NULL
    int CV_ILL_INPUT
    int CV_NO_MALLOC
    int CV_BAD_K
    int CV_BAD_T
    int CV_BAD_DKY
    int CV_TOO_CLOSE
    
    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data)
    ctypedef int (*CVRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data)
    
    void *CVodeCreate(int lmm, int iter)
    int CVodeStep "CVode"(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask) nogil
    int CVodeSetUserData(void *cvode_mem, void *user_data)
    int CVodeSetMaxOrd(void *cvode_mem, int maxord)
    int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
    int CVodeSetMaxHnilWarns(void *cvode_mem, int mxhnil)
    int CVodeSetStabLimDet(void *cvode_mem, booleantype stldet)
    int CVodeSetInitStep(void *cvode_mem, realtype hin)
    int CVodeSetMinStep(void *cvode_mem, realtype hmin)
    int CVodeSetMaxStep(void *cvode_mem, realtype hmax)
    int CVodeSetStopTime(void *cvode_mem, realtype tstop)
    int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef)
    int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor)
    int CVodeSetMaxConvFails(void *cvode_mem, int maxncf)
    int CVodeSetNonlinConvCoef(void *cvode_mem, realtype nlscoef)
    int CVodeSetIterType(void *cvode_mem, int iter)
    int CVodeSetRootDirection(void *cvode_mem, int *rootdir)
    int CVodeSetNoInactiveRootWarn(void *cvode_mem)
    int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
    int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0)
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
    int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky)
    int CVodeGetWorkSpace(void *cvode_mem, long int *lenrw, long int *leniw)
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps)
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevals)
    int CVodeGetNumLinSolvSetups(void *cvode_mem, long int *nlinsetups)
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails)
    int CVodeGetLastOrder(void *cvode_mem, int *qlast)
    int CVodeGetCurrentOrder(void *cvode_mem, int *qcur)
    int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred)
    int CVodeGetActualInitStep(void *cvode_mem, realtype *hinused)
    int CVodeGetLastStep(void *cvode_mem, realtype *hlast)
    int CVodeGetCurrentStep(void *cvode_mem, realtype *hcur)
    int CVodeGetCurrentTime(void *cvode_mem, realtype *tcur)
    int CVodeGetTolScaleFactor(void *cvode_mem, realtype *tolsfac)
    int CVodeGetErrWeights(void *cvode_mem, N_Vector eweight)
    int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals)
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
    int CVodeGetIntegratorStats(void *cvode_mem, long int *nsteps,
                                long int *nfevals, long int *nlinsetups,
                                long int *netfails, int *qlast,
                                int *qcur, realtype *hinused, realtype *hlast,
                                realtype *hcur, realtype *tcur)
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters)
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails)
    int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, long int *nncfails)
    char *CVodeGetReturnFlagName(int flag)
    void CVodeFree(void **cvode_mem)
    
    int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals)
    int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nrevalsLS)

cdef extern from "sundials/sundials_direct.h":
    cdef struct _DlsMat:
        pass
        
    ctypedef _DlsMat *DlsMat
    cdef realtype* DENSE_COL(DlsMat A, int j)

cdef extern from "cvode/cvode_dense.h":
    int CVDense(void *cvode_mem, int N)
    ctypedef int (*CVDlsDenseJacFn)(int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn djac)
    
cdef extern from "ida/ida.h":
    int  IDA_NORMAL
    int  IDA_ONE_STEP
    int  IDA_YA_YDP_INIT
    int  IDA_Y_INIT
    int  IDA_SUCCESS
    int  IDA_TSTOP_RETURN
    int  IDA_ROOT_RETURN
    int  IDA_WARNING
    int  IDA_MEM_NULL
    int  IDA_ILL_INPUT
    int  IDA_NO_MALLOC
    int  IDA_TOO_MUCH_WORK
    int  IDA_TOO_MUCH_ACC
    int  IDA_ERR_FAIL
    int  IDA_CONV_FAIL
    int  IDA_LINIT_FAIL
    int  IDA_LSETUP_FAIL
    int  IDA_LSOLVE_FAIL
    int  IDA_RES_FAIL
    int  IDA_CONSTR_FAIL
    int  IDA_REP_RES_ERR
    int  IDA_MEM_FAIL
    int  IDA_BAD_T
    int  IDA_BAD_EWT
    int  IDA_FIRST_RES_FAIL
    int  IDA_LINESEARCH_FAIL
    int  IDA_NO_RECOVERY
    int  IDA_RTFUNC_FAIL
    
    ctypedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    ctypedef int (*IDARootFn)(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *user_data)
    
    void *IDACreate()
    int IDASetUserData(void *ida_mem, void *user_data)
    int IDASetMaxOrd(void *ida_mem, int maxord)
    int IDASetMaxNumSteps(void *ida_mem, long int mxsteps)
    int IDASetInitStep(void *ida_mem, realtype hin)
    int IDASetMaxStep(void *ida_mem, realtype hmax)
    int IDASetStopTime(void *ida_mem, realtype tstop)
    int IDASetNonlinConvCoef(void *ida_mem, realtype epcon)
    int IDASetMaxErrTestFails(void *ida_mem, int maxnef)
    int IDASetMaxNonlinIters(void *ida_mem, int maxcor)
    int IDASetMaxConvFails(void *ida_mem, int maxncf)
    int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
    int IDASetId(void *ida_mem, N_Vector id)
    int IDASetConstraints(void *ida_mem, N_Vector constraints)
    int IDASetRootDirection(void *ida_mem, int *rootdir)
    int IDASetNoInactiveRootWarn(void *ida_mem)
    int IDAInit(void *ida_mem, IDAResFn res, realtype t0, N_Vector yy0, N_Vector yp0)
    int IDAReInit(void *ida_mem, realtype t0, N_Vector yy0, N_Vector yp0)
    int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol)
    int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol)
    int IDASetNonlinConvCoefIC(void *ida_mem, realtype epiccon)
    int IDASetMaxNumStepsIC(void *ida_mem, int maxnh)
    int IDASetMaxNumJacsIC(void *ida_mem, int maxnj)
    int IDASetMaxNumItersIC(void *ida_mem, int maxnit)
    int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
    int IDASetStepToleranceIC(void *ida_mem, realtype steptol)
    int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g)
    int IDACalcIC(void *ida_mem, int icopt, realtype tout1)
    int IDASolve(void *ida_mem, realtype tout, realtype *tret, N_Vector yret, N_Vector ypret, int itask)
    int IDAGetSolution(void *ida_mem, realtype t,  N_Vector yret, N_Vector ypret)
    int IDAGetWorkSpace(void *ida_mem, long int *lenrw, long int *leniw)
    int IDAGetNumSteps(void *ida_mem, long int *nsteps)
    int IDAGetNumResEvals(void *ida_mem, long int *nrevals)
    int IDAGetNumLinSolvSetups(void *ida_mem, long int *nlinsetups)
    int IDAGetNumErrTestFails(void *ida_mem, long int *netfails)
    int IDAGetNumBacktrackOps(void *ida_mem, long int *nbacktr)
    int IDAGetConsistentIC(void *ida_mem, N_Vector yy0_mod, N_Vector yp0_mod)
    int IDAGetLastOrder(void *ida_mem, int *klast)
    int IDAGetCurrentOrder(void *ida_mem, int *kcur)
    int IDAGetActualInitStep(void *ida_mem, realtype *hinused)
    int IDAGetLastStep(void *ida_mem, realtype *hlast)
    int IDAGetCurrentStep(void *ida_mem, realtype *hcur)
    int IDAGetCurrentTime(void *ida_mem, realtype *tcur)
    int IDAGetTolScaleFactor(void *ida_mem, realtype *tolsfact)
    int IDAGetErrWeights(void *ida_mem, N_Vector eweight)
    int IDAGetEstLocalErrors(void *ida_mem, N_Vector ele)
    int IDAGetNumGEvals(void *ida_mem, long int *ngevals)
    int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals)
    int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS)
    int IDAGetRootInfo(void *ida_mem, int *rootsfound)
    int IDAGetIntegratorStats(void *ida_mem, long int *nsteps, 
                                          long int *nrevals, long int *nlinsetups, 
                                          long int *netfails, int *qlast, int *qcur, 
                                          realtype *hinused, realtype *hlast, realtype *hcur, 
                                          realtype *tcur)
    int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters)
    int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails)
    int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters,  long int *nncfails)
    char *IDAGetReturnFlagName(int flag)
    void IDAFree(void **ida_mem)

cdef extern from "ida/ida_dense.h":
    int IDADense(void *ida_mem, int Neq)
    
    ctypedef int (*IDADlsDenseJacFn)(int Neq, realtype tt, realtype cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn djac)

cdef extern from "kinsol/kinsol.h":
    int KIN_SUCCESS
    int KIN_INITIAL_GUESS_OK
    int KIN_STEP_LT_STPTOL
    int KIN_WARNING
    int KIN_MEM_NULL
    int KIN_ILL_INPUT
    int KIN_NO_MALLOC
    int KIN_MEM_FAIL
    int KIN_LINESEARCH_NONCONV
    int KIN_MAXITER_REACHED
    int KIN_MXNEWT_5X_EXCEEDED
    int KIN_LINESEARCH_BCFAIL
    int KIN_LINSOLV_NO_RECOVERY
    int KIN_LINIT_FAIL
    int KIN_LSETUP_FAIL
    int KIN_LSOLVE_FAIL
    int KIN_SYSFUNC_FAIL
    int KIN_FIRST_SYSFUNC_ERR
    int KIN_REPTD_SYSFUNC_ERR
    int KIN_ETACHOICE1
    int KIN_ETACHOICE2
    int KIN_ETACONSTANT
    int KIN_NONE
    int KIN_LINESEARCH
    
    ctypedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void *user_data)
    
    void *KINCreate()
    int KINSetUserData(void *kinmem, void *user_data)
    int KINSetPrintLevel(void *kinmemm, int printfl)
    int KINSetNumMaxIters(void *kinmem, long int mxiter)
    int KINSetNoInitSetup(void *kinmem, booleantype noInitSetup)
    int KINSetNoResMon(void *kinmem, booleantype noNNIResMon)
    int KINSetMaxSetupCalls(void *kinmem, long int msbset)
    int KINSetMaxSubSetupCalls(void *kinmem, long int msbsetsub)
    int KINSetEtaForm(void *kinmem, int etachoice)
    int KINSetEtaConstValue(void *kinmem, realtype eta)
    int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha)
    int KINSetResMonParams(void *kinmem, realtype omegamin, realtype omegamax)
    int KINSetResMonConstValue(void *kinmem, realtype omegaconst)
    int KINSetNoMinEps(void *kinmem, booleantype noMinEps)
    int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep)
    int KINSetMaxBetaFails(void *kinmem, long int mxnbcf)
    int KINSetRelErrFunc(void *kinmem, realtype relfunc)
    int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
    int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
    int KINSetConstraints(void *kinmem, N_Vector constraints)
    int KINSetSysFunc(void *kinmem, KINSysFn func)
    int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl)
    int KINSol(void *kinmem, N_Vector uu, int strategy, N_Vector u_scale, N_Vector f_scale)
    int KINGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw)
    int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters)
    int KINGetNumFuncEvals(void *kinmem, long int *nfevals)
    int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails) 
    int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr)
    int KINGetFuncNorm(void *kinmem, realtype *fnorm)
    int KINGetStepLength(void *kinmem, realtype *steplength)
    char *KINGetReturnFlagName(int flag)
    void KINFree(void **kinmem)

cdef extern from "kinsol/kinsol_direct.h":
    int KINDense(void *kinmem, int N)