# -*- coding: utf-8 -*-
#cython: boundscheck=False
import os
import time

if os.name == "nt":
    timer = time.clock
else:
    timer = time.time
    
from sundials cimport *
from sundials import CVodeRootException

cdef int f(realtype t, void *yv, void *yvdot, void* user_data):
    cdef double *y = N_Vector_data(yv)
    cdef double *ydot = N_Vector_data(yvdot)
    
    ydot[0] = -0.04*y[0] + 1.0e4*y[1]*y[2]
    ydot[2] = 3.0e7*y[1]*y[1]
    ydot[1] = -ydot[0] - ydot[2]
    
    return 0
    
def test_ODE():
    cdef CVodeSolver solver
    cdef CVodeIterator iter
    cdef N_Vector_Serial y
    cdef realtype t, t0, dt
    
    t0 = 0
    y0 = [1.,0.,0.]
    solver = CVodeSolver(RHS = object(),
                        lmm = "bdf", iter = "newton", 
                        abstol = [1.0e-8, 1.0e-14, 1.0e-6], reltol = 1.0e-4)
    
    solver.RhsFn = <void *>f
    
    solver.init(t0,y0)
    
    start = timer() 
    
    t = 0.4
    while t < 4e11:
        try:
            y = solver.step(t)
        except CVodeRootException, info:
            #print 'root'
            #print "%.4e  %.6e  %.6e  %.6e" % (info.t, info.y[0], info.y[1], info.y[2])
            continue
        
        print "%.4e  %.6e  %.6e  %.6e" % (t, y[0], y[1], y[2])
        t *= 10
    
    print 'timer : ', timer() - start
    print solver.stat()

    
if __name__ == "__main__":
    test_ODE()
