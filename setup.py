#!/usr/bin/env python

import sys
import textwrap

from distutils.core import setup
from distutils.core import Extension


min_numpy_version = '1.4.1'
min_scipy_version = '0.10.0'
min_cython_version = '0.15'
min_mpi4py_version = '1.0'

# Some functions for showing errors and warnings.
def _print_admonition(kind, head, body):
    tw = textwrap.TextWrapper(
        initial_indent='   ', subsequent_indent='   ')

    print(".. %s:: %s" % (kind.upper(), head))
    for line in tw.wrap(body):
        print(line)

def exit_with_error(head, body=''):
    _print_admonition('error', head, body)
    sys.exit(1)

def print_warning(head, body=''):
    _print_admonition('warning', head, body)


# Check for Python
if not (sys.version_info[0] >= 2 and sys.version_info[1] >= 6):
    exit_with_error("You need Python 2.6.x or Python 2.7.x to install lsdmap package!")

if (sys.version_info[0] >= 3 and sys.version_info[1] >= 0):
    exit_with_error("You need Python 2.6.x or Python 2.7.x to install lsdmap package!")

# Check for required Python packages
def check_import(pkgname, pkgver):
    try:
        mod = __import__(pkgname)
    except ImportError:
            exit_with_error(
                "You need %(pkgname)s %(pkgver)s or greater to run lsdmap!"
                % {'pkgname': pkgname, 'pkgver': pkgver} )
    else:
        if mod.__version__ < pkgver:
            exit_with_error(
                "You need %(pkgname)s %(pkgver)s or greater to run lsdmap!"
                % {'pkgname': pkgname, 'pkgver': pkgver} )

    print(( "* Found %(pkgname)s %(pkgver)s package installed."
            % {'pkgname': pkgname, 'pkgver': mod.__version__} ))
    globals()[pkgname] = mod

check_import('numpy', min_numpy_version)
check_import('scipy', min_scipy_version)
check_import('cython', min_cython_version)
check_import('mpi4py', min_mpi4py_version)

import numpy as np
from Cython.Distutils import build_ext

ext_modules = [Extension(
    name='lsdmap/util/cqcprot',
    sources=["lsdmap/util/cqcprot.pyx", "lsdmap/util/qcprot.c"],
    include_dirs = [np.get_include()],  # .../site-packages/numpy/core/include
    extra_compile_args=["-O3","-ffast-math"],
    )]

setup(name='lsdmap',
      version='2.0',
      packages=['lsdmap', 'lsdmap.util'],
      scripts = ['bin/lsdmap'],
      ext_modules = ext_modules,
      cmdclass = {'build_ext': build_ext},
      license='LICENSE.txt',
      description='LSDMap package',
      long_description=open('README.txt').read(),
     )
