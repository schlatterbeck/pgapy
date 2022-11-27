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

import pytest
import pga
import sys

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
    pytest.mpi_rank = t.mpi_rank
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
