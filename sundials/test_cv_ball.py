# -*- coding: utf-8 -*-
"""
Bouncing ball example from OpenModelical manual
"""
import array
import matplotlib
import pylab

from sundials import *

e = 0.7
g = 9.81
TINY = 1e-16

def f(t,y,sw):
    if sw[0]:
        return [0., 0.]
    else:
        return [-g, y[0]]
    
def rootf(t,y,sw):
    return [y[1] + TINY]

t0 = 0
y0 = array.array('d', [0.,1.])
solver = CVodeSolver(RHS = f, ROOT = rootf, SW = [False],
               abstol = 1.0e-6, reltol = 1.0e-6)

solver.init(t0,y0)

dt = .01
iter = solver.iter(t0, dt)
tres = [t0]
hres = [y0[1]]
while True:
    try:
        t, y = next(iter)
        tres.append(t)
        hres.append(y[1])
            
    except CVodeRootException as info:
        if abs(info.y[0]) < 0.01:
            solver.SW[0] = True
            
        tres.append(info.t)
        hres.append(0.)
        
        solver.init(info.t, [-e*info.y[0], 0.])
    
    if t > 3.0:
        break

pylab.plot(tres,hres, linewidth=1.0)
pylab.xlabel('time (s)')
pylab.ylabel('height (m)')
pylab.title('Bouncing ball demo')
pylab.grid(True)
pylab.show()
