#!/usr/bin/python3

import pga
import sys
import weakref
import warnings
import numpy as np
from gp import Terminal, Genetic_Programming, Function, Const
from gp import F_add, F_sub, F_mul, F_div, F_sin, F_cos, F_sqrt
from random import Random
from argparse import ArgumentParser
from scipy.optimize import curve_fit
from scipy.integrate import fixed_quad

# In some cases the automatically-generated expression will produce
# warnings when curve-fitting, we ignore those.
msg = 'Covariance of the parameters could not be estimated'
warnings.filterwarnings ("ignore", message = msg)

class PGA_Random (Random) :

    def __init__ (self, pga_instance) :
        self.pga_instance = weakref.ref (pga_instance)
        super ().__init__ (1)
    # end def __init__

    def getstate (self, *args, **kw) :
        raise NotImplementedError ("Getting the state is not possible")
    # end def getstate

    def setstate (self, *args, **kw) :
        raise NotImplementedError ("Setting the state is not possible")
    # end def setstate

    def random (self) :
        return self.pga_instance ().random01 ()
    # end def random

# end class PGA_Random

def good_solution ():
    """ Build a good solution to integral (sin (k * x))
    """
    tree = F_sub ()
    tree.add_child (F_div (Terminal ('a'), Terminal ('k')))
    fd = F_div ()
    tree.add_child (fd)
    fd.add_child (F_cos (F_mul (Terminal ('x'), Terminal ('k'))))
    fd.add_child (Terminal ('k'))
    return tree
# end def good_solution

def medium_solution ():
    """ Build a medium solution to integral (sin (k * x))
    """
    tree = F_div ()
    sub  = F_sub ()
    tree.add_child (sub)
    tree.add_child (Terminal ('k'))
    sub.add_child \
        (F_cos (F_sub (F_mul (Terminal ('x'), Terminal ('k')), Terminal ('c'))))
    sub.add_child \
        (F_mul (F_sub (Terminal ('a'), Terminal ('c')), F_sin (Terminal ('b'))))
    return tree
# end def good_solution

class Find_Integral (pga.PGA, Genetic_Programming):

    def __init__ (self, args):
        self.args    = args
        self.randpop = []
        self.popsize = args.popsize
        self.random  = PGA_Random (self)
        self.param   = {}
        terms = [Terminal ('x')]
        if args.use_constant_1:
            terms.append (Const ('1'))
        for k in range (self.args.n_terminals):
            name = chr (ord ('a') + k)
            terms.append (Terminal (name))
            self.param [name] = None
        lower = 1 + bool (args.use_constant_1)
        if 'k' in self.args.function:
            terms.append (Terminal ('k'))
        funcs = \
            [F_add, F_sub, F_mul, F_div, F_sin, F_cos, F_sqrt]
        Genetic_Programming.__init__ (self, funcs, terms)
        if 'k' in self.args.function:
            self.prange = self.terminals [lower:-1]
        else:
            self.prange = self.terminals [lower:]
        d = dict \
            ( maximize        = False
            , pop_size        = self.popsize
            , num_replace     = self.popsize - self.popsize // 10
            , tournament_size = args.tournament_size
            , mutation_prob   = 0.01
            , crossover_prob  = 1.0
            , max_GA_iter     = args.generations
            , random_seed     = args.random_seed
            , print_options   = [pga.PGA_REPORT_STRING]
            , print_frequency = 1
            , tournament_with_replacement = False
            #, no_duplicates   = True
            )
        if self.args.multiobjective:
            d ['num_replace']    = self.popsize
            d ['num_eval']       = 2
            d ['num_constraint'] = 0
        pga.PGA.__init__ (self, Function, 10, **d)
        npoints = 50
        self.X = np.arange (0, np.pi + np.pi / npoints, np.pi / npoints)
        def func (x, k = 1):
            return eval (self.args.function)
        self.Y = []
        r = [fixed_quad (func, 0, x, n = 9) [0] for x in self.X]
        self.Y.append (np.array (r))
        if 'k' in self.args.function:
            for j in range (4):
                r = [fixed_quad (func, 0, x, args = (j+2,), n = 9) [0]
                      for x in self.X
                    ]
                self.Y.append (np.array (r))
    # end def __init__

    def eval_solution (self):
        tree = good_solution ()
        print (tree)
        s = self.evaluate (0, 0, tree)
        self.print_string (sys.stdout, 0, 0, tree)
        print (s)
        tree = medium_solution ()
        print (tree)
        s = self.evaluate (0, 0, tree)
        self.print_string (sys.stdout, 0, 0, tree)
        print (s)
    # end def eval_solution

    def eval (self, x, *args):
        for n, t in enumerate (self.prange):
            try:
                t.value = args [n]
            except IndexError:
                t.value = self.param [t.name]
        if 'k' in self.args.function:
            self.terminals [-1].value = self.k
        if np.isscalar (x):
            self.terminals [0].value = x
            v = self.individuum.eval ()
        else:
            v = []
            for xx in x:
                self.terminals [0].value = xx
                v.append (self.individuum.eval ())
            v = np.array (v)
        return v
    # end def eval

    def evaluate (self, p, pop, tree = None):
        self.phenotype (p, pop, tree = tree)
        s = 0
        max_hits = np.size (self.Y)
        hits = 0
        for n, Y in enumerate (self.Y):
            self.k = n + 1.0
            v = np.abs (Y - self.eval (self.X))
            hits += (v < 1e-3).sum ()
            s += v.sum ()
        points = self.individuum.n_func + self.individuum.n_terminals
        if self.individuum.p_eval is not None and self.args.min_change:
            chg = abs (self.individuum.p_eval - s)
            chg /= self.individuum.p_eval
            if chg < self.args.min_change:
                if self.args.multiobjective:
                    return 1e5, 1e5
                return 1e5
        if self.args.multiobjective:
            return s, max (points, 20)
        return s
    # end def evaluate

    def phenotype (self, p, pop, tree = None):
        if tree is None:
            tree = self.get_gene (p, pop)
        self.individuum = tree
        self.k = 1.0
        if len (self.prange) >= 1:
            p0 = np.ones (len (self.prange))
            try:
                popt, pcov = curve_fit (self.eval, self.X, self.Y [0], p0 = p0)
            except RuntimeError:
                popt = p0
            assert len (popt) == len (self.prange)
            for n, t in enumerate (self.prange):
                self.param [t.name] = popt [n]
    # end def phenotype

    def print_string (self, file, p, pop, tree = None):
        print ("%d iterations" % self.GA_iter, file = file)
        self.phenotype (p, pop, tree)
        fmts = []
        args = []
        for t in self.prange:
            fmt = '%s=%%.11g' % t.name
            fmts.append (fmt)
            args.append (self.param [t.name])
        if args:
            print ('; '.join (fmts) % tuple (args))
        file.flush ()
        if not tree:
            ev = self.get_evaluation (p, pop)
            if np.isscalar (ev):
                ev = [ev]
            print \
                ( 'evaluation: %s' % ', '.join ('%.8g' % x for x in ev)
                , file = file
                )
            file.flush ()
            super ().print_string (file, p, pop)
            print (file = file)
        if self.args.verbose:
            print (self.individuum.as_dot (), file = file)
        file.flush ()
    # end def print_string

    def stop_cond (self):
        best_eval = self.get_best_report (pga.PGA_OLDPOP, 0)
        if best_eval <= 1e-4:
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Find_Integral

def main ():
    cmd = ArgumentParser () 
    cmd.add_argument \
        ( '-r', '--random-seed'
        , help    = "Seed for random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    # Other possible functions:
    #    Original from [1]: np.cos (x) + 2 * x + 1
    #    half-circle sinus-like function using fourier series:
    #        np.sqrt (np.pi * x - x * x) * np.sin (k * x)
    cmd.add_argument \
        ( '-f', '--function'
        , help    = "Function to integrate, default=%(default)s"
        , default = "np.sin (k * x)"
        )
    cmd.add_argument \
        ( '-g', '--generations'
        , help    = "Maximum number of generations, default=%(default)s"
        , type    = int
        , default = 50
        )
    cmd.add_argument \
        ( '--min-change'
        , help    = "Minimum change compared to parent, default=%(default)s"
        , type    = float
        , default = 0
        )
    cmd.add_argument \
        ( '-M', '--multiobjective'
        , help    = "Use multi-objective optimization"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-p', '--popsize'
        , help    = "Population size, default=%(default)s"
        , type    = int
        , default = 2000
        )
    cmd.add_argument \
        ( '-t', '--n-terminals'
        , help    = "Number of terminals for fitting, default=%(default)s"
        , type    = int
        , default = 3
        )
    cmd.add_argument \
        ( '-T', '--tournament-size'
        , help    = "Number of individuals in tournament, default=%(default)s"
        , type    = int
        , default = 7
        )
    cmd.add_argument \
        ( '-v', '--verbose'
        , help    = "Print graphviz (.dot) output when printing a string"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-1', '--use-constant-1'
        , help    = "Use the constant 1 as terminal"
        , action  = 'store_true'
        )
    args = cmd.parse_args ()
    if args.n_terminals > 10:
        print \
            ( "Max number of terminals for fitting = 10, got %s"
            % args.n_terminals
            )
        
    ga = Find_Integral (args)
    #ga.eval_solution ()
    #print (ga.mpi_rank, sys.getrefcount (ga))
    ga.run ()
    #print (ga.mpi_rank, sys.getrefcount (ga))
    #ga = None
# end def main

if __name__ == '__main__':
    main ()
