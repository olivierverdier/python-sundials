# -*- coding: utf-8 -*-
import os
import time

if os.name == "nt":
    timer = time.clock
else:
    timer = time.time
    
#import numpy
from array import array

from sundials import *

def f(t,y,sw, ret = array('d', [0.])):
    ret[0] = 1.
    return ret
    
def rootf(t,y,sw, ret = array('d', [0.])):
    ret[0] = y[0] - .5
    return ret

t0 = 0
y0 = [0.]
solver = CVode(RHS = f, ROOT = rootf, SW = [False],
              abstol = 1.0e-6, reltol = 1.0e-6)

solver.init(t0,y0)
print solver.settings


while True:
    try:
        y = solver.step(5.)
    except CVodeRootException, info:
        print 'root t = ', info.t
        continue
    else:
        print 't = 5., y = ', y
        break

print

solver.init(t0, y0)

dt = .01
iter = solver.iter(t0, dt)
start = timer()
while True:
    try:
        t, y = next(iter)
    except CVodeRootException, info:
        print 'root found at t =', info.t
    
    #print "t = %f, y = %f" % (t, y[0])
    
    if iter.t > 5.0:
        break

print 'timer : ', timer() - start
print
print solver.stat()
print
