class KINSOLError(SundialsError):
    MAP = { \
            KIN_MEM_NULL                : "The KINSOL memory block was not initialized",
            KIN_ILL_INPUT               : "An input argument to KINInit has an illegal value.",
            KIN_NO_MALLOC               : "The KINSOL memory was not allocated",
            KIN_MEM_FAIL                : "A memory allocation request has failed.",
            KIN_LINESEARCH_NONCONV      : "The line search algorithm was unable converge",
            KIN_MAXITER_REACHED         : "The maximum number of nonlinear iterations has been reached.",
            KIN_MXNEWT_5X_EXCEEDED      : "Five consecutive steps have been taken above mxnewtstep",
            KIN_LINESEARCH_BCFAIL       : "The line search algorithm was unable to satisfy the 'beta-condition'",
            KIN_LINSOLV_NO_RECOVERY     : "The user-supplied routine psolve encountered a recoverable error",
            KIN_LINIT_FAIL              : "The linear solver initialization routine (linit) encountered an error.",
            KIN_LSETUP_FAIL             : "The user-supplied routine pset encountered an unrecoverable error.",
            KIN_LSOLVE_FAIL             : "Error in psolve or  linear solver routine",
            KIN_STEP_LT_STPTOL          : "Failed stopping tolerance on scaled step length",
            KIN_SYSFUNC_FAIL            : "The system function failed in an unrecoverable manner.",
            KIN_FIRST_SYSFUNC_ERR       : "The system function failed recoverably at the first call.",
            KIN_REPTD_SYSFUNC_ERR       : "The system function had repeated recoverable errors.",
          }
         
    def __init__(self, value = 0, msg = None):
        self.value = value
        self.msg = msg
        
    def __str__(self):
        if self.value == 0:
            return repr(self.msg)
        else:
            return repr(self.MAP[self.value])

cdef class KINSOLStatistics:
    """Class to encapsulate statistics for solver"""
    def __repr__(self):
        s = []
        s.append("nfevals = %d"      % self.nfevals)
        s.append("nniters = %d"      % self.nniters)
        
        return "(%s)" % (", ".join(s))
    
    def __str__(self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())
        
cdef class KINSOLSettings:
    """Class to encapsulate settings for KINSOL"""
    def __init__(self, **kwargs):
        # set default
        self.NVAR = -1
        self.NEQ = -1
        self.printfl = 0
        self.strategy = 'linesearch'
        self.mxiter = 200
        self.msbset = 1
        self.fnormtol = float('inf')
        self.scsteptol = float('inf')
        
        self.F = None
        self.x0 = None
        self.cstr = None
        
        # process arguments
        self.update(**kwargs)
        
        if self.strategy not in ('none','linesearch'):
            raise KINSOLError(0, "'strategy' parameter either 'none' or 'linesearch'")
            
        if self.F is None:
            raise KINSOLError(0, "Function must be set")
        
        if self.x0 is None:
            raise KINSOLError(0, "x0 values must be set")
        else:
            self.NVAR = len(self.x0)
            
        if self.cstr is None:
            if self.NEQ == -1:
                raise KINSOLError(0, "cstr or NEQ value must be set")
        else:
            if self.NEQ == -1:
                self.NEQ = self.cstr.size
            else:
                if self.NEQ != self.cstr.size:
                    raise KINSOLError(0, "cstr and NEQ must match")
            
    def __repr__(self):
        s = []
        s.append("printfl = %d" % self.printfl)
        s.append("strategy = '%s'" % self.strategy)
        s.append("mxiter = %d" % self.mxiter)
        s.append("msbset = %d" % self.msbset)
        s.append("fnormtol = %e" % self.fnormtol)
        s.append("scsteptol = %e" % self.scsteptol)
        
        return "(%s)" % (", ".join(s))
    
    def __str__(self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())
        
    def update(self, **kwargs):
        cdef int size, i
        cdef double dval
        
        for name, value in kwargs.iteritems():
            if name == 'NEQ':
                self.NEQ = value
                
            elif name == 'x0':
                size = len(value)
                self.x0 = N_Vector_Serial(size)
                for i in range(size):
                    self.x0[i] = value[i]
            
            elif name == 'cstr':
                size = len(value)
                self.cstr = N_Vector_Serial(size)
                for i in range(size):
                    dval = value[i]
                    if dval != -2. and dval != -1. and dval != 0. and \
                       dval != 1. and dval != 1.:
                        raise KINSOLSettings(0, "KINSOLSettings: Value for cstr either -2,-1,0,1,2")
                        
                    self.cstr[i] = value[i]
                    
            else:
                try:
                    setattr(self, name, value)
                except AttributeError:
                    raise KINSOLSettings(0, "KINSOLSettings: Unknown parameter: '%s'" % name)

cdef int func(N_Vector u, N_Vector f, void *user_data):
    cdef object[double, ndim=1, mode="c"] buf
    
    cdef kinsol_userdata *ud = <kinsol_userdata *>user_data
    
    cdef double *u_data = (<N_VectorContent_Serial>u.content).data
    cdef double *f_data = (<N_VectorContent_Serial>f.content).data
    
    cdef int u_size = (<N_VectorContent_Serial>u.content).length
    cdef int f_size = (<N_VectorContent_Serial>f.content).length
    
    cdef double dval
    cdef int i
    
    # copy vars to U vector
    cdef N_Vector_Serial U = <N_Vector_Serial>ud.U
    for i in range(u_size):
        U.setitem(i, u_data[i])
    
    # evaluate python callback
    y = (<object>ud.F)(U)
    
    # write back function value
    try:
        dval = y
    except TypeError:
        try:
            buf = y
        except TypeError:
            for i in range(f_size):
                f_data[i] = y[i]
        else:
            for i in range(f_size):
                f_data[i] = buf[i]
    else:
        f_data[0] = dval
    
    return 0
    
cdef class KINSOLSolver:
    """Class to wrap KINSOL solver"""
    def __init__(self, settings = None, **kwargs):
        self.thisptr = KINCreate()
        if self.thisptr == NULL:
            raise KINSOLError(KIN_MEM_FAIL)
        
        if settings is None:
            self.settings = KINSOLSettings(**kwargs)
        else:
            self.settings = settings
            
            if kwargs:
                settings.update(**kwargs)
        
        self.userdata = kinsol_userdata()
    
    def stat(self):
        """
        Return solver statistics
        """
        cdef KINSOLStatistics stat = KINSOLStatistics()
    
        KINGetNumFuncEvals(self.thisptr, &stat.nfevals)
        KINGetNumNonlinSolvIters(self.thisptr, &stat.nniters)
        
        return stat
        
    def solve(self):
        cdef KINSOLSettings settings = self.settings
        cdef N_Vector_Serial u_scale, f_scale
        cdef int i, size, flag, strategy
        
        # set userdata
        self.userdata.F = <void *>settings.F
        self.u = N_Vector_Serial(settings.NVAR)
        self.userdata.U = <void *>self.u
        
        flag = KINSetUserData(self.thisptr, <void *>&self.userdata)
        if flag != KIN_SUCCESS:
            raise KINSOLError(flag)
        
        # setup constraint
        if not settings.cstr is None:
            flag = KINSetConstraints(self.thisptr, <N_Vector>settings.cstr.thisptr)
        else:
            flag = KINSetConstraints(self.thisptr, NULL)
            
        if flag != KIN_SUCCESS:
            raise KINSOLError(flag)
        
        # set tolerances if given
        if settings.fnormtol != float('inf'):
            flag = KINSetFuncNormTol(self.thisptr, settings.fnormtol)
            if flag != KIN_SUCCESS:
                raise KINSOLError(flag)
        
        if settings.scsteptol != float('inf'):
            flag = KINSetScaledStepTol(self.thisptr, settings.scsteptol)
            if flag != KIN_SUCCESS:
                raise KINSOLError(flag)
        
        # set 'mxiter'
        flag = KINSetNumMaxIters(self.thisptr, settings.mxiter)
        if flag != KIN_SUCCESS:
            raise KINSOLError(flag)
        
        # set 'msbset'
        flag = KINSetMaxSetupCalls(self.thisptr, settings.msbset)
        if flag != KIN_SUCCESS:
            raise KINSOLError(flag)
            
        # setup solver
        self.x = N_Vector_Serial(settings.NVAR)
        for i in range(settings.NVAR):
            self.x[i] = settings.x0[i]
            
        flag = KINInit(self.thisptr, <KINSysFn>func, <N_Vector>settings.x0.thisptr)
        if flag != KIN_SUCCESS:
            raise KINSOLError(flag)
            
        # setup matrix type
        flag = KINDense(self.thisptr, settings.NEQ)
        if flag < 0:
            raise KINSOLError(0, "KINDense function failed!")
        
        # no scaling
        self.u_scale = N_Vector_Serial(settings.NEQ)
        for i in range(settings.NEQ):
            self.u_scale[i] = 1.
        
        self.f_scale = N_Vector_Serial(settings.NEQ)
        for i in range(settings.NEQ):
            self.f_scale[i] = 1.
            
        # run solver
        if settings.strategy == "none":
            strategy = KIN_NONE
        else:
            strategy = KIN_LINESEARCH
        
        flag = KINSol(self.thisptr, <N_Vector>self.x.thisptr, strategy,
                      <N_Vector>self.u_scale.thisptr, <N_Vector>self.f_scale.thisptr)
        
        if flag not in (KIN_SUCCESS, KIN_INITIAL_GUESS_OK):
            raise KINSOLError(flag)
        
        return self.x
                                        
    def __dealloc__(self):
        if self.thisptr != NULL:
            KINFree(&self.thisptr)
