import os
import sys
import re
import textwrap

from distutils.core import setup
from distutils.core import Extension
from distutils.sysconfig import get_python_lib

min_numpy_version = '1.4.1'
min_scipy_version = '0.9.0'
min_mpi4py_version = '1.0'
min_cython_version = '0.20'

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
        if len(mod.__version__)>6:
                mod_ver = mod.__version__[:6]
        else:
                mod_ver = mod.__version__

        def mycmp(version1, version2):
            def normalize(v):
                return [int(x) for x in re.sub(r'(\.0+)*$','', v).split(".")]
            return cmp(normalize(version1), normalize(version2))

        # code for mycmp() taken from http://stackoverflow.com/questions/1714027/version-number-comparison
        if mycmp(mod_ver,pkgver) < 0:
            exit_with_error(
                "You need %(pkgname)s %(pkgver)s or greater to run lsdmap!"
                % {'pkgname': pkgname, 'pkgver': pkgver} )

    print(( "* Found %(pkgname)s %(pkgver)s package installed."
            % {'pkgname': pkgname, 'pkgver': mod.__version__} ))
    globals()[pkgname] = mod

check_import('numpy', min_numpy_version)
check_import('scipy', min_scipy_version)
check_import('mpi4py', min_mpi4py_version)
check_import('cython', min_cython_version)

import numpy as np
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# Handle cython modules
from Cython.Distutils import build_ext
from Cython.Build import cythonize
cmdclass = {'build_ext': build_ext}

ext_modules = [Extension(
    name='lsdmap/util/pyqcprot',
    sources=["lsdmap/util/pyqcprot.pyx"],
    include_dirs=[numpy_include],
    ), Extension(
    name='lsdmap/util/util',
    sources=["lsdmap/util/util.pyx"],
    include_dirs=[numpy_include],
    )]

setup(name='lsdmap',
      packages=['lsdmap', 'lsdmap.mpi', 'lsdmap.rw', 'lsdmap.util', 'lsdmap.rbf', 'dmdmd', 'dmdmd.tools'],
      scripts = ['bin/lsdmap','bin/dmdmd', 'bin/rbffit','bin/reweighting','bin/selection','bin/p_mdrun'],
      ext_modules = cythonize(ext_modules),
      cmdclass = cmdclass,
      license='LICENSE.txt',
      description='LSDMap package',
      long_description=open('README.md').read(),
     )
