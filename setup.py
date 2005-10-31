#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# Copyright (C) 2005 Dr. Ralf Schlatterbeck Open Source Consulting.
# Reichergasse 131, A-3411 Weidling.
# Web: http://www.runtux.com Email: office@runtux.com
# All rights reserved
# ****************************************************************************
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# ****************************************************************************


from distutils.core import setup, Extension
from Version        import VERSION

license = 'GNU Library or Lesser General Public License (LGPL)'

description = \
    """\
        pgapack, the parallel genetic algorithm library (see
        ftp://info.mcs.anl.gov/pub/pgapack/README) is a powerfull genetic
        algorithm library by D. Levine, Mathematics and Computer Science
        Division Argonne National Laboratory. The library is written in
        C. PGAPy wraps this library for use with Python.

        You really, absolutely need to read the documentation that comes
        with pgapack -- and of course you need the pgapack library.
    """

module1 = Extension \
    ( 'pga'
    , sources       = ['pgamodule.c']
    , define_macros = [('WL', '32')]
    , include_dirs  = ['.']
    , libraries     = ['pgaO'] # you might need to adapt name of pga lib
    )

setup \
    ( name             = 'PGAPy'
    , version          = VERSION
    , description      = 'Python wrapper for pgapack, the parallel genetic '\
                         'algorithm library'
    , long_description = description
    , ext_modules      = [module1]
    , author           = "Ralf Schlatterbeck"
    , author_email     = "rsc@runtux.com"
    , classifiers      = \
        [ 'Development Status :: 3 - Alpha'
        , 'Intended Audience :: Developers'
        , 'Intended Audience :: Education'
        , 'Intended Audience :: Science/Research'
        , 'License :: OSI Approved :: ' + license
        , 'Operating System :: OS Independent'
        , 'Programming Language :: C'
        , 'Programming Language :: Python'
        , 'Topic :: Software Development :: Libraries :: Python Modules'
        ]
    )
