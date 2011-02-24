# -*- coding: utf-8 -*-
"""
/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/29 22:21:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 */
 """
import os
import time

if os.name == "nt":
    timer = time.clock
else:
    timer = time.time
    
from array import array

from sundials import *

def f(t,y,sw = None, ret = array('d', [0.,0.,0.])):
    ret[0] = -0.04*y[0] + 1.0e4*y[1]*y[2]
    ret[2] = 3.0e7*y[1]*y[1]
    ret[1] = -ret[0] - ret[2]
    return ret
    
def rootf(t,y,sw, ret = array('d', [0.,0.])):
    ret[0] = y[0] - 0.0001
    ret[1] = y[2] - 0.01
    return ret

def jacf(t,y,sw, ret = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]):
    ret[0][0] = -0.04
    ret[0][1] = 1.0e4*y[2]
    ret[0][2] = 1.0e4*y[1]
    ret[1][0] = 0.04
    ret[1][1] = -1.0e4*y[2] - 6.0e7*y[1]
    ret[1][2] = -1.0e4*y[1]
    ret[2][1] = 6.0e7*y[1]
    return ret
    
solver = CVodeSolver(RHS = f, ROOT = rootf, SW = [False,False], JAC = jacf,
               lmm = "bdf", iter = "newton", 
               abstol = [1.0e-8, 1.0e-14, 1.0e-6], reltol = 1.0e-4)

t0 = 0
y0 = [1., 0., 0.]
solver.init(t0,y0)

start = timer()

t = 0.4
while t < 4e11:
    try:
        y = solver.step(t)
    except CVodeRootException, info:
        print 'root'
        print "%.4e  %.6e  %.6e  %.6e" % (info.t, info.y[0], info.y[1], info.y[2])
        continue
    
    print "%.4e  %.6e  %.6e  %.6e" % (t, y[0], y[1], y[2])
    t *= 10

print 'timer : ', timer() - start
print
print solver.stat()
