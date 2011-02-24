# -*- coding: utf-8 -*-
# -----------------------------------------------------------------
# $Revision: 1.2 $
# $Date: 2008/12/17 19:38:48 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# -----------------------------------------------------------------
# Example (serial):
#
# This example solves a nonlinear system from.
#
# Source: "Handbook of Test Problems in Local and Global Optimization",
#             C.A. Floudas, P.M. Pardalos et al.
#             Kluwer Academic Publishers, 1999.
# Test problem 4 from Section 14.1, Chapter 14: Ferraris and Tronconi
# 
# This problem involves a blend of trigonometric and exponential terms.
#    0.5 sin(x1 x2) - 0.25 x2/pi - 0.5 x1 = 0
#    (1-0.25/pi) ( exp(2 x1)-e ) + e x2 / pi - 2 e x1 = 0
# such that
#    0.25 <= x1 <=1.0
#    1.5 <= x2 <= 2 pi
# 
# The treatment of the bound constraints on x1 and x2 is done using
# the additional variables
#    l1 = x1 - x1_min >= 0
#    L1 = x1 - x1_max <= 0
#    l2 = x2 - x2_min >= 0
#    L2 = x2 - x2_max >= 0
# 
# and using the constraint feature in KINSOL to impose
#    l1 >= 0    l2 >= 0
#    L1 <= 0    L2 <= 0
# 
# The Ferraris-Tronconi test problem has two known solutions.
# The nonlinear system is solved by KINSOL using different 
# combinations of globalization and Jacobian update strategies 
# and with different initial guesses (leading to one or the other
# of the known solutions).
#
#
# Constraints are imposed to make all components of the solution
# positive.
# -----------------------------------------------------------------

from math import e, pi, sin, exp
from array import array
from sundials import KINSOLSolver

NVAR    = 6
NEQ     = 6

x1_min = .25
x1_max = 1.0
x2_min = 1.5
x2_max = 2. * pi
    
def eval_f(u, ret = array('d', range(NEQ))):
    x1, x2 = u[0], u[1]
    l1 = u[2]
    L1 = u[3]
    l2 = u[4]
    L2 = u[5]
  
    ret[0] = 0.5*sin(x1*x2) - 0.25 * x2/pi - 0.5*x1
    ret[1] = (1. - 0.25/pi) * ( exp(2. * x1) - e) + e*x2 / pi - 2.*e*x1
    ret[2] = l1 - x1 + x1_min
    ret[3] = L1 - x1 + x1_max
    ret[4] = l2 - x2 + x2_min
    ret[5] = L2 - x2 + x2_max

    return ret
  
XINIT = [
        x1_min,
        x2_min,
        x1_min - x1_min,
        x1_min - x1_max,
        x2_min - x2_min,
        x2_min - x2_max,
]

CONSTR = [
        0.,     # no constraint on x1
        0.,     # no constraint on x2
        1.,     # l1 = x1 - x1_min >= 0
        -1.,    # L1 = x1 - x1_max <= 0
        1.,     # l2 = x2 - x2_min >= 0
        -1.,    # L2 = x2 - x22_min <= 0
]
                
# this init. guess should take us to (0.29945; 2.83693) */
solver = KINSOLSolver(F = eval_f, x0 = XINIT, cstr = CONSTR,
                strategy = 'linesearch',
                msbset = 1,
                fnormtol = 1.e-5,
                scsteptol = 1.e-5)

try:
    solver.solve()
except:
    print 'error ', sys.exc_value 
    
print 'solution x = ', solver.x
print 
print 'expected x = (0.29945, 2.83693)'
print
print solver.stat()
