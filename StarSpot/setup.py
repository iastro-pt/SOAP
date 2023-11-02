#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy

#Change those variables to meet your system installation
scisoft_gsl_include_dir = ['/usr/local/scisoft/packages/gsl/include/']
macport_gsl_include_dir = ['/opt/local/include/']
linux_gsl_include_dir = ['/usr/local/include/', '/usr/include/']

scisoft_gsl_lib_dir = ['/usr/local/scisoft/packages/gsl/lib/']
macport_gsl_lib_dir = ['/opt/local/lib/', '/opt/local/lib64/']
linux_gsl_lib_dir = ['/usr/lib64', '/usr/lib/']

module1 = Extension('starspot',
                    include_dirs=['.'] + [numpy.get_include()] + scisoft_gsl_include_dir + macport_gsl_include_dir + linux_gsl_include_dir,
                    library_dirs=scisoft_gsl_lib_dir + macport_gsl_lib_dir + linux_gsl_lib_dir,
                    libraries=['gsl', 'gslcblas'],
                    sources=['starspotmodule.c', 'starspot.c']
                    )

setup(name='starspot',
      version='2.0',
      description='Star spot simulator',
      author='X. Dumusque',
      author_email='x.dumusque@gmail.com',
      url='https://www.cfa.harvard.edu/~xdumusqu',
      ext_modules=[module1])
