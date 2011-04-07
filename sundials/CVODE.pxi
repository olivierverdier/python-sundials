class CVodeError(SundialsError):
    def __init__(self, value = 0, msg = None):
        self.value = value
        self.msg = msg
        
    def __str__(self):
        if self.value >= 0:
            return repr(self.msg)
        else:
            return CVodeGetReturnFlagName(self.value)

class CVodeRootException(CVodeError):
    def __init__(self, realtype t, y, SW):
        self.t = t
        self.y = y
        self.SW = SW
        
    def __str__(self):
        args = self.t, repr(self.SW)
        return "CVodeRoot(t = %f, SW = %s)" % args
        
cdef class CVodeSettings:
    """Class to encapsulate settings for CVode"""
    def __init__(self, **kwargs):
        # set default
        self.lmm = "adams"
        self.iter = "functional"
        self.maxord = -1
        self.mxsteps = 500
        self.tstop = float('+inf')
        self.hmax = float('+inf')
        
        self.reltol = -1.0
        self.abstol = -1.0
        self.abstolv = None
        
        self.RHS = None
        self.ROOT = None
        self.SW = None
        self.JAC = None
        
        # process arguments
        self.update(**kwargs)
        
        if self.lmm not in ('adams','bdf'):
            raise CVodeError(0, "CVodeSettings: 'lmm' parameter either 'adams' or 'bdf'")
        
        if self.iter not in ('functional', 'newton'):
            raise CVodeError(0, "CVodeSettings: 'iter' parameter either 'functional' or 'newton'")
            
        if self.maxord == -1:
            if self.lmm == "adams":
                self.maxord = 12
            else:
                self.maxord = 5
                
        if self.reltol == -1.0:
            raise CVodeError(0, "CVodeSettings: 'reltol' parameter must be given")
        
        if self.abstol == -1.0 and self.abstolv is None:
            raise CVodeError(0, "CVodeSettings: 'abstol' parameter must be given")
        
        if not self.SW is None:
            self.SW = [bool(value) for value in self.SW]
    
    def __repr__(self):
        s = []
        s.append("lmm = '%s'" % self.lmm)
        s.append("iter = '%s'" % self.iter)
        s.append("maxord = %d" % self.maxord)
        s.append("mxsteps = %d" % self.mxsteps)
        s.append("tstop = %g" % self.tstop)
        s.append("hmax = %g" % self.hmax)
        s.append("reltol = %g" % self.reltol)
        s.append("abstol = %s" % (self.abstolv or self.abstol))
        
        return "(%s)" % (", ".join(s))
    
    def __str__(self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())
        
    def update(self, **kwargs):
        cdef int size, i
        
        for name, value in kwargs.iteritems():
            if name == 'abstol':
                try:
                    size = len(value)
                except TypeError:
                    self.abstol = value
                else:
                    self.abstolv = N_Vector_Serial(size)
                    for i in range(size):
                        self.abstolv[i] = value[i]
                        
            else:
                try:
                    setattr(self, name, value)
                except AttributeError:
                    raise CVodeError(0, "CVodeSettings: Unknown parameter: '%s'" % name)
    
cdef int cv_rhs(realtype t, N_Vector yv, N_Vector yvdot, void* user_data) with gil:
    """
    Wraps  Python rhs-callback function to obtain CVode required interface
    see also ctypedef statement above
    """
    cdef object[double, ndim=1, mode="c"] buf
    cdef cv_userdata *ud = <cv_userdata *>user_data
    
    cdef double *data = (<N_VectorContent_Serial>yvdot.content).data
    cdef int i, size = (<N_VectorContent_Serial>yvdot.content).length
    
    # copy data
    cdef N_Vector_Serial Y = <N_Vector_Serial>ud.Y
    for i in range(Y.getsize()):
        Y.setitem(i, (<N_VectorContent_Serial>yv.content).data[i])
    
    if not ud.SW == NULL:
        ydot = (<object>ud.RHSF)(t, Y, <object>ud.SW)
    else:
        ydot = (<object>ud.RHSF)(t, Y)
    
    try:
        buf = ydot
    except TypeError:
        for i in range(size):
            data[i] = ydot[i]
    else:
        for i in range(size):
            data[i] = buf[i]
    
    return 0
        
cdef int cv_root(realtype t, N_Vector yv, realtype *gout,  void* user_data) with gil:
    """
    Wraps  Python root-callback function to obtain CV required interface
    see also ctypedef statement above
    """
    cdef object[double, ndim=1, mode="c"] buf
    cdef cv_userdata *ud = <cv_userdata *>user_data
    cdef int i
    
    # copy data
    cdef N_Vector_Serial Y = <N_Vector_Serial>ud.Y
    for i in range(Y.getsize()):
        Y.setitem(i, (<N_VectorContent_Serial>yv.content).data[i])
        
    out = (<object>ud.ROOTF)(t, Y, <object>ud.SW)
    
    try:
        buf = out
    except TypeError:
        for i in range(len(out)):
            gout[i] = out[i]
    else:
        for i in range(len(out)):
            gout[i] = buf[i]
            
    return 0

cdef int cv_jac(int Neq, realtype t, N_Vector yv, N_Vector fy, DlsMat Jac,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)  with gil:
    """
    Wraps Python jacobian-callback function to obtain CV required interface.
    """
    cdef object[double, ndim=2, mode="c"] buf
    cdef cv_userdata *ud = <cv_userdata *>user_data
    cdef realtype* col_i
    cdef int i, j
    
    # copy data
    cdef N_Vector_Serial Y = <N_Vector_Serial>ud.Y
    for i in range(Y.getsize()):
        Y.setitem(i, (<N_VectorContent_Serial>yv.content).data[i])
    
    if not ud.SW == NULL:
        jac = (<object>ud.JACF)(t, Y, <object>ud.SW)
    else:
        jac = (<object>ud.JACF)(t, Y)
    
    try:
        buf = jac
    except TypeError:
        # list/tuple etc.
        for i in range(Neq):
            col_i = DENSE_COL(Jac, i)
            for j in range(Neq):
                col_i[j] = jac[j][i]
    else:
        # fast c array
        for i in range(Neq):
            col_i = DENSE_COL(Jac, i)
            for j in range(Neq):
                col_i[j] = buf[j,i]
                
    return 0

cdef class CVodeIterator:
    """
    Single step iterator
    """
    def __init__(self, solver, realtype t0, realtype dt):
        self.solver = solver
        self.dt = dt
        self.t = t0 + dt
        
        # ensure first step is run
        self.tret = t0
        self.hlast = 0.
        
        self.root = False
        self.stop = False
        self.last = False
        
    def __iter__(self):
        return self
    
    cdef N_Vector_Serial Next(self):
        """
        Interpolate result at time t and return
        solution vector. The time must not be
        past the current internal time of the solver.
        """
        cdef N_Vector_Serial y
        cdef realtype tcur
        cdef int flag
        
        y = self.solver.y
        
        while True:
            if self.root and self.t >= self.troot:
                self.root = False
                self.solver.__handleRoot(self.troot, self.yroot)
                
            elif (self.tret - self.hlast) <= self.t <= self.tret:
                tcur = self.t
                self.t += self.dt
                
                if self.stop:
                    # last step
                    self.last = True
                
                flag = CVodeGetDky(self.solver.thisptr, tcur, 0, <N_Vector>y.thisptr)
                if flag != CV_SUCCESS:
                    raise CVodeError(flag)
                
                break
                
            else:
                flag = CVodeStep(self.solver.thisptr, self.t, <N_Vector>y.thisptr,
                         &self.tret, CV_ONE_STEP)
            
                if flag < 0:
                    raise CVodeError(flag)
                
                CVodeGetLastStep(self.solver.thisptr, &self.hlast)
                
                if flag == CV_ROOT_RETURN:
                    self.root = True
                    self.troot = self.tret
                    self.yroot = y.copy()
                
                elif flag == CV_TSTOP_RETURN:
                    self.stop = True
                    
                continue
        
        return y
        
    def __next__(self):
        cdef realtype tcur
        
        if self.last:
            raise StopIteration()
        
        # save current step
        tcur = self.t
        
        return tcur, self.Next()
                
cdef class CVodeSolver:
    """Class to wrap CVode"""
    def __init__(self, settings = None, **kwargs):
        if settings is None:
            self.settings = CVodeSettings(**kwargs)
        else:
            self.settings = settings
            
            if kwargs:
                settings.update(**kwargs)
        
        self.RhsFn = <void *>cv_rhs
        self.RootFn = <void *>cv_root
        self.JacFn = <void *>cv_jac
            
    def __dealloc__(self):
        if self.thisptr != NULL:
            CVodeFree(&self.thisptr)
        
    def stat(self):
        """
        Return solver statistics
        """
        cdef SolverStatistics stat = SolverStatistics()
    
        CVodeGetNumSteps(self.thisptr, &stat.nsteps)
        CVodeGetNumRhsEvals(self.thisptr, &stat.nfevals)
        CVDlsGetNumJacEvals(self.thisptr, &stat.njevals)
        CVodeGetNumGEvals(self.thisptr, &stat.ngevals)
        CVodeGetNumErrTestFails(self.thisptr, &stat.netfails)
        CVodeGetNumNonlinSolvIters(self.thisptr, &stat.nniters)
        CVodeGetNumNonlinSolvConvFails(self.thisptr, &stat.nncfails)
        
        return stat
        
    def init(self, realtype t0, y0, **kwargs):
        """
        Initialize the solver or restart the solver
        """
        cdef CVodeSettings settings
        cdef N_Vector_Serial abstolv
        cdef int i, flag, lmm, iter
        
        # update settings
        settings = self.settings
        if kwargs:
            settings.update(**kwargs)
        
        # store initial values
        self.t0 = t0
        try:
            # create from buffer type: numpy/array etc.
            self.y = N_Vector_Serial_Copy_Array(y0)
        except TypeError:
            # create from list/tuple etc.
            y0 = [float(value) for value in y0]
            
            self.y = N_Vector_Serial(len(y0))
            for i in range(self.y.size):
                self.y[i] = y0[i]
        
        # Initialization
        if self.thisptr == NULL:
            lmm = 1 if settings.lmm == 'adams' else 2
            iter = 1 if settings.iter == 'functional' else 2
            
            self.thisptr = CVodeCreate(lmm, iter)
            if self.thisptr == NULL:
                raise CVodeError(CV_MEM_FAIL)
            
            flag = CVodeInit(self.thisptr, <CVRhsFn>self.RhsFn, self.t0, <N_Vector>self.y.thisptr)
            if flag != CV_SUCCESS:
                raise CVodeError(flag)
            
            flag = CVodeSetMaxOrd(self.thisptr, settings.maxord)
            if flag != CV_SUCCESS:
                raise CVodeError(flag)
                
            flag = CVDense(self.thisptr, self.y.size)
            if flag != CV_SUCCESS:
                raise CVodeError(flag)
            
            self.RHSF = settings.RHS
            self.userdata = cv_userdata(<void*>settings.RHS, <void*>NULL,
                                        <void*>NULL, <void*>NULL, <void*>NULL)
            
            if not settings.ROOT is None:
                self.ROOTF = settings.ROOT
                self.userdata.ROOTF = <void*>settings.ROOT
            
                if not settings.SW is None:
                    self.SW = settings.SW
                    self.userdata.SW = <void*>settings.SW
                else:
                    raise CVodeError(0, "Expected switch vector")
                
                flag = CVodeRootInit(self.thisptr, len(settings.SW), <CVRootFn>self.RootFn)
                if flag != CV_SUCCESS:
                    raise CVodeError(flag)
            
            if not settings.JAC is None:
                self.JACF = settings.JAC
                self.userdata.JACF = <void*>settings.JAC
                
                flag = CVDlsSetDenseJacFn(self.thisptr, <CVDlsDenseJacFn>self.JacFn)
                if flag != CV_SUCCESS:
                    raise CVodeError(flag)
                    
            self.yv = self.y.copy()
            self.userdata.Y = <void*>self.yv
        
        # Reinitialization
        else:
            flag = CVodeReInit(self.thisptr, self.t0, <N_Vector>self.y.thisptr)
            if flag != CV_SUCCESS:
                raise CVodeError(flag)
            
            if self.RHSF != settings.RHS:
                self.RHSF = settings.RHS
                self.userdata.RHSF = <void*>settings.RHS
            
            if self.ROOTF != settings.ROOT:
                self.ROOTF = settings.ROOT
                self.userdata.ROOTF = <void*>settings.ROOT
            
            if self.SW != settings.SW:
                self.SW = settings.SW
                self.userdata.SW = <void*>settings.SW
            
            if self.JACF != settings.JAC:
                self.JACF = settings.JAC
                self.userdata.JACF = <void*>settings.JAC
        
        if settings.abstolv is None:
            flag = CVodeSStolerances(self.thisptr, settings.reltol, settings.abstol)
        else:
            abstolv = settings.abstolv
            flag = CVodeSVtolerances(self.thisptr, settings.reltol, <N_Vector>abstolv.thisptr)
            
        if flag != CV_SUCCESS:
            raise CVodeError(flag)
        
        flag = CVodeSetMaxNumSteps(self.thisptr, settings.mxsteps)
        if flag != CV_SUCCESS:
            raise CVodeError(flag)
                
        flag = CVodeSetMaxStep(self.thisptr, settings.hmax)
        if flag != CV_SUCCESS:
            raise CVodeError(flag)
        
        if settings.tstop != float('+inf'):
            flag = CVodeSetStopTime(self.thisptr, settings.tstop)
            if flag != CV_SUCCESS:
                raise CVodeError(flag)
        
        flag = CVodeSetUserData(self.thisptr, <void*>&self.userdata)
        if flag != CV_SUCCESS:
            raise CVodeError(flag)
        
    cpdef __handleRoot(self, realtype event_time, N_Vector_Serial y):
        """
        Create root expection with info on time
        and solution vector.
        """
        cdef int *event_info
        cdef int i, size = len(self.settings.SW)
        
        # Allocate memory for the event_info
        event_info = <int *>malloc(size * sizeof(int))
        CVodeGetRootInfo(self.thisptr, event_info)
        
        SW = [bool(event_info[i]) for i in range(size)]
        
        raise CVodeRootException(event_time, y, SW)
    
    cpdef CVodeIterator iter(self, realtype t0, realtype dt):
        """
        Run solver to next root or to 
        tstop if set
        """
        return CVodeIterator(self, t0, dt)
                
    cpdef N_Vector_Serial step(self, realtype tf):
        """
        Run solver to time tf, next root change
        or to tstop if set
        """
        cdef int flag
        cdef realtype tret
        
        flag = CVodeStep(self.thisptr, tf, <N_Vector>self.y.thisptr, &tret, CV_NORMAL)
        
        if flag < 0:
            raise CVodeError(flag)
        
        if flag == CV_ROOT_RETURN:
            self.__handleRoot(tret, self.y.copy())
        
        return self.y.copy()
