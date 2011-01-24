**python-sundials** is a Cython wrapper for the Sundials solver suite.
The wrapper is based on code posten on the cython-dev mailing list by
Mr. Jon Olav Vik.

**Highlights**

 * CVODE - Solver for stiff and nonstiff ordinary differential equation
 * IDA - Solver for the solution of differential-algebraic equation (DAE) systems
 * KINSOL - solver for nonlinear algebraic systems based on Newton-Krylov solver technology

The CVODE and IDA solvers support root finding and the solver throws an exception
on finding a root.

There is also an example of implementing the explicit equation in Cython for speed.