cdef class N_Vector_Serial:
    """
    Wrapper arround sundials N_Vector_Serial
    """
    # buffer interface
    cdef __cythonbufferdefaults__ = {"ndim": 1, "mode": "c"}
    
    def __init__(N_Vector_Serial self, int size):
        assert size > 0
        self.thisptr = N_VNew_Serial(size)
        
    def __dealloc__(N_Vector_Serial self):
        if not self.thisptr == NULL:
            N_VDestroy_Serial(<N_Vector>self.thisptr)
    
    cdef int getsize(N_Vector_Serial self):
        cdef N_Vector v = <N_Vector>self.thisptr
        return (<N_VectorContent_Serial>v.content).length
    
    cdef double *getdata(N_Vector_Serial self):
        cdef N_Vector v = <N_Vector>self.thisptr
        return (<N_VectorContent_Serial>v.content).data
        
    cdef double getitem(N_Vector_Serial self, int index):
        cdef N_Vector v = <N_Vector>self.thisptr
        return (<N_VectorContent_Serial>v.content).data[index]
    
    cdef setitem(N_Vector_Serial self, int index, double value):
        cdef N_Vector v = <N_Vector>self.thisptr
        (<N_VectorContent_Serial>v.content).data[index] = value
        
    def __getbuffer__(N_Vector_Serial self, Py_buffer* buffer, int flags):
        self.__shape[0] = self.getsize()
        
        buffer.buf = <void *>self.getdata()
        buffer.obj = self
        buffer.len = self.getsize() * sizeof(double)
        buffer.readonly = 0
        buffer.format = <char*>"d"
        buffer.ndim = 1
        buffer.shape = <Py_ssize_t *>&self.__shape
        buffer.strides = NULL
        buffer.suboffsets = NULL
        buffer.itemsize = sizeof(double)
        buffer.internal = NULL
        
    def __releasebuffer__(N_Vector_Serial self, Py_buffer* buffer):
        pass
    
    property dtype:
        def __get__(N_Vector_Serial self):
            return "d"
    
    property size:
        def __get__(N_Vector_Serial self):
            return self.getsize()
    
    property itemsize:
        def __get__(N_Vector_Serial self):
            return sizeof(double)
    
    property nbytes:
        def __get__(N_Vector_Serial self):
            return sizeof(double) * self.getsize()
    
    def __repr__(N_Vector_Serial self):
        cdef double *data = self.getdata()
        cdef int size = self.getsize()
        cdef int index
        cdef char *typestr = "%g"
            
        ret = [None] * size
        for index in range(size):
            ret[index] = typestr  % data[index]
        
        return "(%s)" % ",".join(ret)
    
    def __str__(N_Vector_Serial self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())
        
    def __len__(N_Vector_Serial self):
        return self.getsize()
    
    def __getitem__(N_Vector_Serial self, int index):
        assert index < self.getsize()
        return self.getitem(index)
    
    def __setitem__(N_Vector_Serial self, int index, double value):
        assert index < self.getsize()
        self.setitem(index, value)
    
    cpdef copy(N_Vector_Serial self):
        cdef N_Vector_Serial ret = N_Vector_Serial.__new__(N_Vector_Serial, 0)
        cdef double *data = self.getdata()
        cdef int i, size = self.getsize()
        
        ret.thisptr = N_VNew_Serial(size)
        
        for i in range(size):
            ret.setitem(i, data[i])
        
        return ret
    
cdef double *N_Vector_data(void *vector):
    cdef N_Vector v = <N_Vector>vector
    return (<N_VectorContent_Serial>v.content).data

cdef int N_Vector_length(void *vector):
    cdef N_Vector v = <N_Vector>vector
    return (<N_VectorContent_Serial>v.content).length
    
def N_Vector_Serial_Copy_Array(object[double, ndim=1] arg):
    # fast copy construct bypassing __init__
    cdef N_Vector_Serial ret = N_Vector_Serial.__new__(N_Vector_Serial, 0)
    cdef int i, size
    
    size = len(arg)
    ret.thisptr = N_VNew_Serial(size)
    
    for i in range(size):
        ret.setitem(i, arg[i])
    
    return ret
