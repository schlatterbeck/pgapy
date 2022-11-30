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
import os
from filecmp import cmp
from io      import StringIO

# Import from examples
sys.path.insert (1, "examples")
from cards            import main as cards_main
from cards_mutate     import main as cards_mutate_main
from constraint       import main as constraint_main
from dtlz2            import main as dtlz2_main
from fourbar          import main as fourbar_main
from gears            import main as gears_main
from hello_world_char import main as hello_world_char_main
from hello_world_int  import main as hello_world_int_main
from himmelblau       import main as himmelblau_main
from magic_prio       import main as magic_prio_main
from magic_square     import main as magic_square_main
from minfloat         import main as minfloat_main
from multi            import main as multi_main
from one_max          import main as one_max_main
from sort_numbers     import main as sort_numbers_main
from twobar           import main as twobar_main
from vibr             import main as vibr_main
from namefull         import main as namefull_main

skip_tsplib = skip_fann = lambda fun, *args, **kw: fun

try:
    from xor      import main as xor_main
except ImportError as err:
    skip_fann = pytest.mark.skip (reason = str (err))

sys.path.insert (1, "examples/sequence")
try:
    from tsp      import main as tsp_main
except ImportError as err:
    skip_tsplib = pytest.mark.skip (reason = str (err))

sys.path.insert (1, "examples/gp")
from opt_xor      import main as gp_xor_main
from opt_parity3  import main as gp_parity3_main
from opt_integral import main as gp_integral_main

class _Test_PGA:
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
        if request.node.rep_call.passed:
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
# end class _Test_PGA

class Test_PGA_Fast (_Test_PGA):
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

    def test_dtlz2_scaled (self):
        opt= '-R 42 -r --scale=10'.split ()
        dtlz2_main (self.out_options + opt)
        self.compare ()
    # end def test_dtlz2_scaled

    def test_fourbar (self):
        fourbar_main (self.out_options + ['-R', '42'])
        self.compare ()
    # end def test_fourbar

    def test_hello_world_char (self):
        hello_world_char_main (self.out_options + ['-R', '42', '-i'])
        self.compare ()
    # end def test_hello_world_char

    def test_hello_world_int (self):
        hello_world_int_main (self.out_options + ['-R', '42'])
        self.compare ()
    # end def test_hello_world_int

    def test_gears (self):
        gears_main (self.out_options + ['-R', '42'])
        self.compare ()
    # end def test_gears

    def test_gears_check (self, capfd):
        if pytest.mpi_rank != 0:
            return
        arglist = '-l 17 -u 90 -n 950 -d 150 -c 23,49,83,86'.split ()
        gears_main (self.out_options + arglist)
        r = 'Factor: 6.333333\nGear Error:  0.004670060%\n'
        captured = capfd.readouterr ()
        assert r == captured.out
    # end def test_gears_check

    def test_himmelblau (self):
        himmelblau_main (self.out_options + ['-R', '1'])
        self.compare ()
    # end def test_himmelblau

    def test_magic_prio (self):
        # Use length 3, doesn't find a solution otherwise :-)
        magic_prio_main (self.out_options + '-m 0.111111'.split ())
        self.compare ()
    # end def test_magic_prio

    def test_magic_prio_ed (self):
        # Use length 3, doesn't find a solution otherwise :-)
        opt = '-m 0.111111 --use-euclidian-gene-distance'.split ()
        magic_prio_main (self.out_options + opt)
        self.compare ()
    # end def test_magic_prio_ed

    def test_magic_square (self):
        magic_square_main (self.out_options + '-l 4'.split ())
        self.compare ()
    # end def test_magic_square

    def test_magic_square_ed (self):
        opt = '-l 4 --use-euclidian-gene-distance'.split ()
        magic_square_main (self.out_options + opt)
        self.compare ()
    # end def test_magic_square_ed

    def test_magic_square_cm (self):
        ca = '-l 3 -m 0.111 --use-custom-mutation'.split ()
        magic_square_main (self.out_options + ca)
        self.compare ()
    # end def test_magic_square_cm

    def test_minfloat (self):
        minfloat_main (self.out_options)
        self.compare ()
    # end def test_minfloat

    def test_multi (self):
        multi_main (self.out_options)
        self.compare ()
    # end def test_multi

    def test_one_max (self):
        one_max_main (self.out_options + '-R 23 -v'.split ())
        self.compare ()
    # end def test_one_max

    def test_sort_numbers (self):
        sort_numbers_main (self.out_options + '-R 23'.split ())
        self.compare ()
    # end def test_sort_numbers

    def test_twobar (self):
        twobar_main (self.out_options + '-R 42'.split ())
        self.compare ()
    # end def test_twobar

    def test_vibr (self):
        vibr_main (self.out_options + '-R 42'.split ())
        self.compare ()
    # end def test_vibr

    @skip_fann
    def test_xor (self):
        xor_main (self.out_options)
        self.compare ()
    # end def test_xor

    @skip_fann
    def test_xor_binary (self):
        xor_main (self.out_options + '-b -m 100'.split ())
        self.compare ()
    # end def test_xor_binary

    @skip_fann
    def test_xor_gray (self):
        xor_main (self.out_options + '-g -m 100'.split ())
        self.compare ()
    # end def test_xor_gray

    @skip_tsplib
    def test_tsp_croes (self):
        tsp_main (self.out_options + ['examples/sequence/croes.tsp'])
        self.compare ()
    # end def test_tsp_croes

    @skip_tsplib
    def test_tsp_croes_lk (self, capfd):
        if pytest.mpi_rank != 0:
            return
        a = '--lin-kernighan examples/sequence/croes.tsp'.split ()
        tsp_main (self.out_options + a)
        with open (self.data_name, 'r') as f:
            r = f.read ()
        captured = capfd.readouterr ()
        assert r == captured.out
    # end def test_tsp_croes_lk

    def test_gp_xor (self):
        gp_xor_main (self.out_options + '-R 2 -D'.split ())
        self.compare ()
    # end def test_gp_xor

    def test_gp_xor_verbose (self):
        gp_xor_main (self.out_options + '-R 2 -v'.split ())
        if pytest.mpi_rank != 0:
            return
        with open (self.data_name, 'r') as f:
            dat = f.read ()
        with open (self.out_name, 'r') as f:
            out = f.read ()
        d_seen = False
        for d, o in zip (dat.split ('\n'), out.split ('\n')):
            if d_seen:
                if d == '}':
                    assert o == '}'
                    break
                d_a, d_b = d.split (None, 1)
                o_a, o_b = o.split (None, 1)
                if 'label' in d_b:
                    assert d_b == o_b
                else:
                    assert o_b.startswith ('->')
            else:
                if d.startswith ('digraph'):
                    assert o.startswith ('digraph')
                    d_seen = True
                continue
    # end def test_gp_xor_verbose

    def test_gp_integral (self):
        p = '-T 7 -p 500 -f np.cos(x)+2*x+1'.split ()
        gp_integral_main (self.out_options + p)
        self.compare ()
    # end def test_gp_integral

    def test_getter_setter (self):
        if pytest.mpi_rank != 0:
            return
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (float, 10)
        # end class T
        t = T ()
        assert t.crossover_prob == 0.85
        assert t.crossover_bounce_back == 0
        assert t.crossover_bounded == 0
        assert t.crossover_SBX_eta == 2.0
        assert t.crossover_SBX_once_per_string == 0
        assert t.DE_variant == pga.PGA_DE_VARIANT_RAND
        assert t.DE_num_diffs == 1
        assert t.DE_scale_factor == 0.9
        assert t.DE_aux_factor == 0.95
        assert t.DE_crossover_prob == 0.9
        assert t.DE_jitter == 0
        assert t.DE_dither == 0
        assert t.DE_probability_EO == 0.5
        assert t.epsilon_exponent == 0
        assert t.epsilon_generation == 0
        assert t.epsilon_theta == 0
        assert t.eval_count == 0
        assert t.fitness_cmax == 1.01
        assert t.GA_iter == 0
        assert t.max_fitness_rank == 1.2
        assert t.max_GA_iter == 1000
        assert t.multi_obj_precision == 14
        assert t.mutation_and_crossover == 0
        assert t.mutation_bounce_back == 0
        assert t.mutation_bounded == 0
        assert t.mutation_poly_eta == 100.0
        assert t.mutation_or_crossover == 1
        assert t.mutation_only == 0
        assert t.mutation_poly_value == -1.0
        assert t.mutation_prob == 0.1
        assert t.num_constraint == 0
        assert t.num_replace == 10
        assert t.pop_size == 100
        assert t.print_frequency == 10
        assert t.p_tournament_prob == 0.6
        assert t.randomize_select == 0
        assert t.restart == 0
        assert t.restart_frequency == 50
        assert t.rtr_window_size == 5
        assert t.string_length == 10
        assert t.sum_constraints == 0
        assert t.tournament_size == 2
        assert t.tournament_with_replacement == 1
        assert t.truncation_proportion == 0.5
        assert t.uniform_crossover_prob == 0.6
        assert t.mutation_value == 0.1
        assert t.num_eval == 1
        assert t.mpi_rank == 0
        t.crossover_prob = 0.85
        t.epsilon_exponent = 7.5
        t.multi_obj_precision = 13
        t.p_tournament_prob = 0.9
        t.uniform_crossover_prob = 0.9
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10)
        # end class T
        t = T ()
        assert t.mutation_value == 1
        class T (pga.PGA):
            def __init__ (self):
                d = dict \
                    ( crossover_bounce_back = True
                    , crossover_bounded     = True
                    )
                super ().__init__ (bool, 10, num_eval = 2, **d)
        # end class T
        t = T ()
        assert t.mutation_value is None
        t.set_evaluation (0, pga.PGA_OLDPOP, 0, 1)
        t.set_evaluation_up_to_date (0, pga.PGA_OLDPOP, 1)
        assert 0.5 <= t.random_uniform (0.5, 0.7) <= 0.7
        assert -1000 <= t.random_gaussian (0, 0.01) <= 1000
    # end def test_getter_setter

    def test_init_param (self):
        if pytest.mpi_rank != 0:
            return
        initp = [(i, .5) for i in range (10)]
        d = dict \
            ( crossover_SBX_once_per_string = 1
            , crossover_SBX_eta = 0.5
            , max_similarity = 99
            , DE_num_diffs = 2
            , DE_dither_per_individual = 1
            , DE_aux_factor = 0.99
            , DE_probability_EO = 0.8
            , uniform_crossover_prob = 0.8
            , p_tournament_prob = 0.9
            , max_fitness_rank = 1.1
            , fitness_cmax = 1.1
            , init_percent = initp
            , restart = True
            , restart_frequency = 100
            , mutation_bounded = True
            , mutation_or_crossover = True
            , mutation_and_crossover = True
            , mutation_value = 0.2
            , mutation_poly_eta = 80.0
            , mutation_poly_value = 5
            , fitness_type = pga.PGA_FITNESS_RAW
            , fitness_min_type = pga.PGA_FITNESSMIN_CMAX
            , rtr_window_size = 23
            , sum_constraints = False
            , epsilon_generation = 23
            , epsilon_exponent = 7.5
            , epsilon_theta = 33
            , multi_obj_precision = 10
            , randomize_select = True
            )
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (float, 10, **d)
        t = T ()
        assert t.crossover_SBX_once_per_string == 1
        assert t.crossover_SBX_eta == 0.5
        assert t.max_similarity == 99
        assert t.DE_num_diffs == 2
        assert t.DE_dither_per_individual == 1
        assert t.DE_aux_factor == 0.99
        assert t.DE_probability_EO == 0.8
        assert t.uniform_crossover_prob == 0.8
        assert t.p_tournament_prob == 0.9
        assert t.max_fitness_rank == 1.1
        assert t.fitness_cmax == 1.1
        assert t.restart == True
        assert t.restart_frequency == 100
        assert t.mutation_bounded == True
        assert t.mutation_or_crossover == True
        assert t.mutation_and_crossover == False
        assert t.mutation_value == 0.2
        assert t.mutation_poly_eta == 80.0
        assert t.mutation_poly_value == 5
        assert t.fitness_type == pga.PGA_FITNESS_RAW
        assert t.fitness_min_type == pga.PGA_FITNESSMIN_CMAX
        assert t.rtr_window_size == 23
        assert t.sum_constraints == False
        assert t.epsilon_generation == 23
        assert t.epsilon_exponent == 7.5
        assert t.epsilon_theta == 33
        assert t.multi_obj_precision == 10
        assert t.get_iteration () == 0
        assert t.randomize_select == 1
    # end def test_init_param

    def test_init_param_intmutation (self):
        if pytest.mpi_rank != 0:
            return
        d = dict (mutation_type = pga.PGA_MUTATION_PERMUTE)
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10, **d)
        t = T ()
        assert t.mutation_type == pga.PGA_MUTATION_PERMUTE
    # end def test_init_param_intmutation

    def test_missing_eval (self):
        """ This needs special care with mpi: the rank-0 process will
            not raise the exception because it never calls evaluate.
            But for a single process and two processes the rank-0
            process *will* call evaluate.
        """
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (float, 10)
        t = T ()
        if pytest.mpi_n_proc > 2 and pytest.mpi_rank == 0:
            t.run ()
            return
        with pytest.raises (NotImplementedError):
            t.run ()
    # end def test_missing_eval

    def test_set_allele (self):
        if pytest.mpi_rank != 0:
            return
        # BIN
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (bool, 10)
        t = T ()
        t.set_allele (0, pga.PGA_OLDPOP, 0, True)
        assert t.get_allele (0, pga.PGA_OLDPOP, 0)

        # FLOAT
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (float, 10)
        t = T ()
        t.set_allele (0, pga.PGA_OLDPOP, 0, 4711.0815)
        assert t.get_allele (0, pga.PGA_OLDPOP, 0) == 4711.0815

        # INT
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10)
        t = T ()
        t.set_allele (0, pga.PGA_OLDPOP, 0, 4711)
        assert t.get_allele (0, pga.PGA_OLDPOP, 0) == 4711

        # CHAR
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (bytes, 10)
        t = T ()
        t.set_allele (0, pga.PGA_OLDPOP, 0, b'A')
        assert t.get_allele (0, pga.PGA_OLDPOP, 0) == b'A'

        # USER
        class T (pga.PGA):
            def initstring (self, p, pop):
                pass
            def __init__ (self):
                super ().__init__ (tuple, 10)
        t = T ()
        with pytest.raises (ValueError):
            t.set_allele (0, pga.PGA_OLDPOP, 0, 'A')
        with pytest.raises (ValueError):
            t.get_allele (0, pga.PGA_OLDPOP, 0)
    # def test_set_allele

    def test_bin_gray (self):
        if pytest.mpi_rank != 0:
            return
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (bool, 10)
        t = T ()
        pop = pga.PGA_OLDPOP
        t.encode_int_as_binary (0, pop, 0, 4, 17)
        t.encode_int_as_binary (0, pop, 5, 9, 23)
        i = t.get_int_from_binary (0, pop, 0, 9)
        assert i == 23 + (17 << 5)
        l = 2.0
        u = 100.0
        f = t.get_real_from_binary (0, pop, 0, 9, l, u)
        assert abs (f - (u - l) * i / 1023 + l) < 1e12
        t.encode_real_as_binary (0, pop, 0, 9, l + 1, u + 1, f + 1)
        assert t.get_int_from_binary (0, pop, 0, 9) == 23 + (17 << 5)
        t.encode_int_as_gray_code (0, pop, 0, 4, 17)
        t.encode_int_as_gray_code (0, pop, 5, 9, 23)
        t.get_int_from_gray_code (0, pop, 0, 4) == 17
        t.get_int_from_gray_code (0, pop, 5, 9) == 23
        t.encode_real_as_gray_code (0, pop, 0, 9, l, u, 50.0)
        f = t.get_real_from_gray_code (0, pop, 0, 9, l, u)
        assert abs (f - 50) <= (1 / 1023) * (u - l)
    # end def test_bin_gray

    def test_int_params (self, capfd):
        if pytest.mpi_rank != 0:
            return
        d = dict (mutation_value = 2, random_seed = 23)
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10, **d)
        t = T ()
        assert t.mutation_value == 2
        t.print_context ()
        with open (self.data_name, 'r') as f:
            r = f.read ()
        captured = capfd.readouterr ().err
        lines = []
        for line in captured.split ('\n'):
            if line.strip ().startswith ('PrintString'):
                continue
            if line.strip ().startswith ('Stopping'):
                continue
            if line.strip ().startswith ('Gene distance'):
                continue
            lines.append (line)
        captured = '\n'.join (lines)
        assert r == captured
    # end def test_int_params

    def test_eval_misuse (self):
        if pytest.mpi_n_proc > 1:
            return
        d = dict (num_eval = 3)
        class T (pga.PGA):
            def evaluate (self, p, pop):
                return 2
            def __init__ (self):
                super ().__init__ (int, 10, **d)
        t = T ()
        with pytest.raises (ValueError):
            t.run ()
        class T (pga.PGA):
            def evaluate (self, p, pop):
                return 1, 2, 3, 4
            def __init__ (self):
                super ().__init__ (int, 10, **d)
        t = T ()
        with pytest.raises (ValueError):
            t.run ()
    # end def test_eval_misuse

    def test_check_allele (self):
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10)
        t = T ()
        # wrong population:
        with pytest.raises (ValueError):
            t.get_allele (0, 4711, 0)
        # wrong population index p:
        with pytest.raises (ValueError):
            t.get_allele (500, pga.PGA_OLDPOP, 0)
        # wrong allele index
        with pytest.raises (ValueError):
            t.get_allele (0, pga.PGA_OLDPOP, 10)
    # end def test_check_allele

    def test_print_string (self):
        if pytest.mpi_rank != 0:
            return
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10)
        t = T ()
        for i in range (10):
            t.set_allele (0, pga.PGA_OLDPOP, i, i)
        with open (self.out_name, 'w') as f:
            t.print_string (f, 0, pga.PGA_OLDPOP)
        with open (self.data_name, 'r') as f:
            r = f.read ()
        self.compare ()
    # end def test_print_string

    def test_das_dennis (self):
        if pytest.mpi_rank != 0:
            return
        dd = pga.das_dennis (3, 12, 0.05, [1, 1, 1])
        assert len (dd [0]) == 3
        assert len (dd) == 91
        with pytest.raises (ValueError):
            dd = pga.das_dennis (3, 12, 0.05, [1, 1, 1, 5])
        with pytest.raises (ValueError):
            dd = pga.das_dennis (3, 12, 0.05, 1)
    # end def test_das_dennis

    def test_refpoints (self):
        if pytest.mpi_rank != 0:
            return
        d = dict (num_eval = 3, num_constraint = 0, reference_points = 2)
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10, **d)
        with pytest.raises (ValueError):
            t = T ()
        d.update (reference_points = [])
        with pytest.raises (ValueError):
            t = T ()
        d.update (reference_points = [[1, 2, 3, 4]])
        with pytest.raises (ValueError):
            t = T ()
        d.update (reference_points = [['a', 2, 3]])
        with pytest.raises (ValueError):
            t = T ()
    # end def test_refpoints

    def test_print_option_hamming (self):
        if pytest.mpi_rank != 0:
            return
        d = dict (print_options = [pga.PGA_REPORT_HAMMING])
        class T (pga.PGA):
            def __init__ (self):
                super ().__init__ (int, 10, **d)
        with pytest.raises (ValueError):
            t = T ()
        d = dict (print_options = [''])
        with pytest.raises (TypeError):
            t = T ()
    # end def test_print_option_hamming

# end class Test_PGA_Fast

class Test_PGA_Slow (_Test_PGA):

    def test_gp_parity3 (self):
        gp_parity3_main (self.out_options + '-R 25'.split ())
        self.compare ()
    # end def test_gp_parity3

    def test_gp_integral_multi (self):
        p = '-r 1 -f np.sin(k*x) --multi'.split ()
        gp_integral_main (self.out_options + p)
        self.compare ()
    # end def test_gp_integral_multi

    def test_namefull (self):
        namefull_main (self.out_options + '-R 42'.split ())
        self.compare ()
    # end def test_namefull

    @skip_tsplib
    def test_tsp_linhp318_lk (self):
        opt = '-m 20 -k 0.01 examples/sequence/tsplib/linhp318.tsp'.split ()
        tsp_main (self.out_options + opt)
        self.compare ()
    # end def test_tsp_linhp318_lk

# end class Test_PGA_Slow
