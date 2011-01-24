#!/usr/bin/python
# -*- utf8 -*-
"""
Copyright (C) 2010 Runar Tenfjord <runar.tenfjord@tenko.no>
"""
import sys
import os
import glob
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#sys.argv.append('build_ext')
#sys.argv.extend(['sdist','--formats=gztar,zip'])
sys.argv.append('bdist_wininst')

classifiers = '''\
Environment :: Console
Development Status :: 4 - Beta
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering :: Mathematics
'''

try:
    setup(
      name = 'python-sundials',
      version = '0.1',
      description = 'Sundials solver library wrapper',
      long_description = '''\
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
''',
        
    classifiers = [value for value in classifiers.split("\n") if value],
        
    author='Runar Tenfjord',
    author_email = 'runar.tenfjord@google.com',
    license = 'BSD',
    download_url='http://pypi.python.org/pypi/python-sundials/',
    url = 'http://code.google.com/p/python-sundials/',
    platforms = ['any'],
    
    ext_modules=[   
        Extension("sundials",
                  sources = ["sundials/sundials.pyx",],
                  depends = glob.glob('sundials/*.pxi'),
                  include_dirs = ['sundials', 'sundials/include'],
                  library_dirs=['sundials'],
                  libraries = ['sundials_cvode','sundials_ida',
                               'sundials_kinsol', 'sundials_nvecserial'])
        ],
        
      cmdclass = {'build_ext': build_ext}
    )
except:
    print('Traceback\n:%s\n' % str(sys.exc_info()[-2]))
else:
    print('\n')
    
input('Press enter to continue')