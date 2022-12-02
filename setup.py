#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2005-22 Dr. Ralf Schlatterbeck Open Source Consulting.
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

import sys
from setuptools import setup, Extension
if sys.version_info.major > 2 :
    from subprocess     import run, PIPE
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

# This is used to override the module to install from the module
# definitions below.
modulename = environ.get ('PGA_MODULE', 'module_from_pgapack_submodule')

# The default pgapack module: Use the version from the submodule pgapack
# that comes with pgapy.
module_from_pgapack_submodule = Extension \
    ( 'pga'
    , sources = ['pgamodule.c'] + pgapack_sources + [stub]
    , define_macros = [('FAKE_MPI', '1')]
    , include_dirs  = ['.', 'pgapack/fakempi', 'pgapack/include']
    , depends       = ['pgapack/include/pgapack.h', 'pgapack/fakempi/mpi.h']
    )

# example config for default pga home when installing pga from source
# contributed by Márk Váradi. You need to set the environment-variable
# 'PGA_MODULE' to 'module_from_source' if you want to use this setting.
BASE     = '/usr/local/pga/'	# default pga home in Linux machines
DEB_BASE = '/usr/include/pgapack-serial' # Default on Debian since lenny
module_from_source = Extension \
    ( 'pga'
    , sources       = ['pgamodule.c']
    , define_macros = []
    , include_dirs  = ['.', path.join (BASE, 'include'), DEB_BASE]
    , libraries     = [':libpgaO.a'] # you might need to adapt name of pga lib
    , library_dirs  = [path.join (BASE, 'lib/linux')]
    )

# default config since debian lenny, serial version (installation in /usr):
# Set environment variable PGA_MODULE to 'module_from_install' if you
# want to use that version.
module_from_install = Extension \
    ( 'pga'
    , sources       = ['pgamodule.c']
    , define_macros = []
    , include_dirs  = ['.', '/usr/include/pgapack-serial']
    , libraries     = ['pgapack-serial']
    )

# Parallel version. You need to specify PGA_PARALLEL_VARIANT in the
# environment, the default is 'mpich'.
# Note: The parallel version below is specific to python3 on debian-like
# distributions. Let me know if you want to build a parallel version on
# other architectures and we can figure out how to.
if modulename == 'module_from_parallel_install' and sys.version_info.major > 2 :
    parallel_variant = environ.get ('PGA_PARALLEL_VARIANT', 'mpich')
    result = run (['dpkg-architecture', '-qDEB_BUILD_GNU_TYPE'], stdout = PIPE)
    arch_dir = result.stdout.decode ('utf-8').rstrip ()
    arch_dir = path.join ('/usr/include', arch_dir)
    include_by_variant = dict \
        ( mpich   = path.join (arch_dir, 'mpich')
        , lam     = '/usr/include/lam'
        , openmpi = path.join (arch_dir, 'openmpi')
        )
    module_from_parallel_install = Extension \
        ( 'pga'
        , sources       = ['pgamodule.c']
        , define_macros = []
        , include_dirs  = \
            [ '.'
            , '/usr/include/pgapack-' + parallel_variant
            , include_by_variant [parallel_variant]
            ]
        , libraries     = ['pgapack-' + parallel_variant]
        )

module1 = locals ()[modulename]

setup \
    ( name             = 'PGAPy'
    , version          = VERSION
    , description      = 'Python wrapper for pgapack, the parallel genetic '\
                         'algorithm library'
    , long_description_content_type = 'text/x-rst'
    , long_description = ''.join (description)
    , ext_modules      = [module1]
    , data_files       = [ ( 'share/pgapy/examples'
                           , [ 'examples/cards.py'
                             , 'examples/cards_mutate.py'
                             , 'examples/constraint.py'
                             , 'examples/dtlz2.py'
                             , 'examples/fourbar.py'
                             , 'examples/gears.py'
                             , 'examples/hello_world_char.py'
                             , 'examples/hello_world_int.py'
                             , 'examples/himmelblau.py'
                             , 'examples/magic_prio.py'
                             , 'examples/magic_square.py'
                             , 'examples/minfloat.py'
                             , 'examples/multi.py'
                             , 'examples/namefull.py'
                             , 'examples/one_max.py'
                             , 'examples/sort_numbers.py'
                             , 'examples/twobar.py'
                             , 'examples/vibr.py'
                             , 'examples/xor.py'
                             ]
                           )
                         , ( 'share/pgapy/examples/gp'
                           , [ 'examples/gp/gp.py'
                             , 'examples/gp/opt_integral.py'
                             , 'examples/gp/opt_parity3.py'
                             , 'examples/gp/opt_xor.py'
                             , 'examples/gp/README.rst'
                             ]
                           )
                         , ( 'share/pgapy/examples/sequence'
                           , [ 'examples/sequence/croes.tsp'
                             , 'examples/sequence/oliver30.tsp'
                             , 'examples/sequence/plot_tour.py'
                             , 'examples/sequence/tsp.py'
                             ]
                           )
                         , ( 'share/pgapy'
                           , ['pgapack/docs/user_guide.pdf']
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
