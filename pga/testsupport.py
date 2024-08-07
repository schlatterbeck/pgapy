#!/usr/bin/python3
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
#    documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ****************************************************************************

import os
import pytest
from . import pga
from filecmp import cmp

class PGA_Test_Instrumentation:
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
    def data_name (self):
        return os.path.join ('test', self.basename + '.data')
    # end def data_name

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
        # When test is skipped by calling pytest.skip, the hook is not
        # run and rep_call doesn't exist, so to be more robust here we
        # check if the rep_call attribute exists. No cleanup is
        # necessary if the test was skipped...
        rep_call = getattr (request.node, 'rep_call', None)
        if rep_call and rep_call.passed:
            self.cleanup ()
    # end def run_clean_test

    def cleanup (self):
        if pytest.mpi_rank == 0:
            try:
                os.unlink (self.out_name)
            except OSError:
                pass
    # end def cleanup

    def compare (self):
        if pytest.mpi_rank == 0:
            assert cmp (self.data_name, self.out_name, shallow = False)
    # end def compare
# end class PGA_Test_Instrumentation

@pytest.fixture (scope = "session", autouse = True)
def pga_setup_test (request):
    """ Since MPI_Init may be only called *once* per program and since
        the tests are all run by a single invocation of python we need
        to initialize MPI once in the beginning and finalize it after
        all tests have been run.
    """
    pga.MPI_Init ([])
    request.addfinalizer (pga.MPI_Finalize)
    class T (pga.PGA):
        def __init__ (self):
            super ().__init__ (float, 10)
    # end class T
    t = T ()
    pytest.mpi_rank   = t.mpi_rank
    pytest.mpi_n_proc = t.mpi_n_proc
    del t
    del T
# end def pga_setup_test

@pytest.hookimpl (tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport (item, call):
    # execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result ()
    # set a report attribute for each phase of a call, which can
    # be "setup", "call", "teardown"
    setattr (item, "rep_" + rep.when, rep)
# end def pytest_runtest_makereport

