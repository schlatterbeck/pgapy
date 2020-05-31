#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2005-20 Dr. Ralf Schlatterbeck Open Source Consulting.
# Reichergasse 131, A-3411 Weidling.
# Web: http://www.runtux.com Email: office@runtux.com
# ****************************************************************************
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ****************************************************************************

from distutils.core import setup, Extension
try :
    from Version    import VERSION
except :
    VERSION = None
from textwrap       import dedent
from os             import path, environ
from glob           import glob

description = []
with open ('README.rst') as f :
    description = f.read ()

license = 'BSD License'
pgapack_sources = []
stub = 'pgapack/source/mpi_stub.c'
pgapack_sources.extend (fn for fn in glob ('pgapack/source/*.c') if fn != stub)

module_from_pgapack_submodule = Extension \
    ( 'pga'
    , sources = ['pgamodule.c'] + pgapack_sources + [stub]
    , extra_compile_args = ['-fPIC']
    , define_macros = [('FAKE_MPI', '1')]
    , include_dirs  = ['.', 'pgapack/fakempi', 'pgapack/include']
    , depends       = ['pgapack/include/pgapack.h', 'pgapack/fakempi/mpi.h']
    )

# example config for default pga home when installing pga from source
# contributed by Márk Váradi. You need to comment the module1 below if
# you want to use this setting.
BASE     = '/usr/local/pga/'	# default pga home in Linux machines
DEB_BASE = '/usr/include/pgapack-serial' # Default on Debian since lenny
module_from_source = Extension \
    ( 'pga'
    , sources       = ['pgamodule.c']
    , define_macros = [('WL', '32')]
    , include_dirs  = ['.', path.join (BASE, 'include'), DEB_BASE]
    , libraries     = [':libpgaO.a'] # you might need to adapt name of pga lib
    , library_dirs  = [path.join (BASE, 'lib/linux')]
    )

# default config since debian lenny, serial version (installation in /usr):
# (comment following lines if you want to use a config above)
module_from_install = Extension \
    ( 'pga'
    , sources       = ['pgamodule.c']
    , define_macros = [('WL', '32')]
    , include_dirs  = ['.', '/usr/include/pgapack-serial']
    , libraries     = ['pgapack-serial']
    )

module1 = environ.get ('PGA_MODULE', 'module_from_pgapack_submodule')
module1 = locals ()[module1]

setup \
    ( name             = 'PGAPy'
    , version          = VERSION
    , description      = 'Python wrapper for pgapack, the parallel genetic '\
                         'algorithm library'
    , long_description = ''.join (description)
    , ext_modules      = [module1]
    , data_files       = [ ( 'share/pgapy/examples'
                           , ['examples/xor.py', 'examples/magic_square.py']
                           )
                         ]
    , author           = "Ralf Schlatterbeck"
    , author_email     = "rsc@runtux.com"
    , url              = "http://pgapy.sourceforge.net/"
    , classifiers      = \
        [ 'Development Status :: 5 - Production/Stable'
        , 'Intended Audience :: Developers'
        , 'Intended Audience :: Education'
        , 'Intended Audience :: Science/Research'
        , 'License :: OSI Approved :: ' + license
        , 'Operating System :: OS Independent'
        , 'Programming Language :: C'
        , 'Programming Language :: Python'
        , 'Topic :: Scientific/Engineering :: Artificial Intelligence'
        , 'Topic :: Software Development :: Libraries :: Python Modules'
# Would be nice if distutils supported the following category -- it
# doesn't according to "python setup.py register --list-classifiers"
#       , 'Topic :: Software Development :: Algorithms :: Genetic Algorithms'
        ]
    )
