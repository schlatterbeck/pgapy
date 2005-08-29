#!/usr/bin/env python

from distutils.core import setup, Extension

module1 = Extension \
    ( 'pga'
    , sources       = ['pgamodule.c']
    , define_macros = [('WL', '32')]
    , libraries     = ['pgaO']
    )

setup \
    ( name        = 'PGA'
    , version     = '0.1'
    , description = 'Wrapper for pgapack genetic algorithm library'
    , ext_modules = [module1]
    )
