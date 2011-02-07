# -*- coding: utf-8 -*-
"""
/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2009/01/21 21:46:40 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * This simple example problem for IDA, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:
 *
 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.
 *
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01.
 *
 * The problem is solved with IDA using IDADENSE for the linear
 * solver, with a user-supplied Jacobian. Output is printed at
 * t = .4, 4, 40, ..., 4e10.
 * -----------------------------------------------------------------
 */
"""
from array import array

from sundials import *

def f(t, y, ydot, sw, ret = array('d', [0.,0.,0.])):
    ret[0] = -.04*y[0] + 1.e4*y[1]*y[2]
    ret[1] = -ret[0] - 3.e7*y[1]*y[1] - ydot[1]
    ret[0] -=  ydot[0]
    ret[2] = y[0] + y[1] + y[2] - 1.
    
    return ret

def r(t, y, ydot, sw, ret = array('d', [0.,0.])):
    ret[0] = y[0] - 0.0001
    ret[1] = y[2] - 0.01
    
    return ret

def jacf(c,t,y,ydot,sw, ret = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]):
    ret[0][0] = -0.04 - c
    ret[1][0] = 0.04
    ret[2][0] = 1.
    ret[0][1] = 1.0e4*y[2]
    ret[1][1] = -1.0e4*y[2] - 6.0e7*y[1] - c
    ret[2][1] = 1.
    ret[0][2] = 1.0e4*y[1]
    ret[1][2] = -1.0e4 * y[1]
    ret[2][2] = 1.
    return ret
    
t0 = 0
y0 = [1.,0.,0.]
yd0 = [-0.04, 0.04, 0.]

solver = IDA(RES = f, ROOT = r, SW = [False,False], JAC = jacf, mxsteps = 5000,
             abstol = [1.0e-8, 1.0e-14, 1.0e-6], reltol = 1.0e-4)

solver.init(t0,y0,yd0)

t = t0 + .4
while t < 4e11:
    try:
        y, ydot = solver.step(t)
    except IDARootException, info:
        print 'root:'
        print "%.4e  %.4e  %.4e  %.4e" % (info.t, info.y[0], info.y[1], info.y[2])
        continue
        
    print "%.4e  %.4e  %.4e  %.4e" % (t, y[0], y[1], y[2])
    t *= 10

print
print solver.stat()
