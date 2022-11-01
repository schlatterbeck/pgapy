# Copyright (C) 2022 Dr. Ralf Schlatterbeck Open Source Consulting.
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

import unittest
import pytest
import pga
import sys
import os
from filecmp import cmp
from io      import StringIO

# Import from examples
sys.path.insert (1, "examples")
from cards        import main as cards_main
from cards_mutate import main as cards_mutate_main
from constraint   import main as constraint_main
from dtlz2        import main as dtlz2_main
from fourbar      import main as fourbar_main
from gears        import main as gears_main
from hello_world  import main as hello_world_main
from himmelblau   import main as himmelblau_main

class Test_PGA (unittest.TestCase):
    out_name    = 'test/output.out'

    @property
    def basename (self):
        assert self.test_name.startswith ('test_')
        return self.test_name [5:]
    # end def basename

    @property
    def out_name (self):
        return os.path.join ('test', self.basename + '.out')
    # end def out_name

    @property
    def out_options (self):
        return ['-O', self.out_name]
    # end def out_options

    @pytest.fixture (autouse = True)
    def run_clean_test (self, request):
        self.test_name = request.node.name
        self.cleanup ()
        yield
        # Only clean up if successful, see pytest_runtest_makereport
        # in conftest.py
        if request.node.rep_call.passed:
            self.cleanup ()
    # end def run_clean_test

    def cleanup (self):
        try:
            os.unlink (self.out_name)
        except OSError:
            pass
    # end def cleanup

    def compare (self):
        fn = os.path.join ('test', self.basename + '.data')
        assert cmp (fn, self.out_name, shallow = False)
    # end def compare

    def test_cards (self):
        cards_main (self.out_options + ['-R', '2'])
        self.compare ()
    # end def test_cards

    def test_cards_mutate (self):
        cards_mutate_main (self.out_options)
        self.compare ()
    # end def test_cards_mutate

    def test_constraint (self):
        constraint_main (self.out_options)
        self.compare ()
    # end def test_constraint

    def test_dtlz2 (self):
        dtlz2_main (self.out_options + ['-R', '42'])
        self.compare ()
    # end def test_dtlz2

    def test_fourbar (self):
        fourbar_main (self.out_options + ['-R', '42'])
        self.compare ()
    # end def test_fourbar

    def test_hello_world (self):
        hello_world_main (self.out_options + ['-R', '42'])
        self.compare ()
    # end def test_hello_world

    def test_gears (self):
        gears_main (self.out_options + ['-R', '42'])
        self.compare ()
    # end def test_gears

    def test_gears_check (self):
        sio = StringIO ()
        so  = sys.stdout
        sys.stdout = sio
        arglist = '-l 17 -u 90 -n 950 -d 150 -c 23,49,83,86'.split ()
        gears_main (self.out_options + arglist)
        sys.stdout = so
        r = 'Factor: 6.333333\nGear Error:  0.004670060%\n'
        self.assertEqual (r, sio.getvalue ())
    # end def test_gears_check

    def test_himmelblau (self):
        himmelblau_main (self.out_options + ['-R', '1'])
        self.compare ()
    # end def test_himmelblau

# end class Test_PGA
