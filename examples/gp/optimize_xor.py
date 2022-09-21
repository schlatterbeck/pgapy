#!/usr/bin/python3

import pga
import sys
from gp import F_nand, Terminal, Genetic_Programming, Function
from random import Random
from argparse import ArgumentParser

class PGA_Random (Random) :

    def __init__ (self, pga_instance) :
        self.pga_instance = pga_instance
        super ().__init__ (1)
    # end def __init__

    def getstate (self, *args, **kw) :
        raise NotImplementedError ("Getting the state is not possible")
    # end def getstate

    def setstate (self, *args, **kw) :
        raise NotImplementedError ("Setting the state is not possible")
    # end def setstate

    def random (self) :
        return self.pga_instance.random01 ()
    # end def random

# end class PGA_Random


class Find_XOR (pga.PGA, Genetic_Programming):

    def __init__ (self, args):
        self.args    = args
        self.randpop = []
        self.popsize = 500
        self.random  = PGA_Random (self)
        terms = [Terminal ('a'), Terminal ('b')]
        Genetic_Programming.__init__ (self, [F_nand], terms)
        pga.PGA.__init__ \
            ( self, Function, 10
            , maximize      = False
            , pop_size      = self.popsize
            , num_replace   = self.popsize - self.popsize // 10
            , mutation_prob = 0.0
            , max_GA_iter   = 50
            , random_seed   = args.random_seed
            , print_options = [pga.PGA_REPORT_STRING]
            )
    # end def __init__

    def mutation (self, p, pop, pm):
        """ No mutation """
        return 0
    # end def mutation

    def crossover (self, p1_i, p2_i, ppop, c1_i, c2_i, cpop):
        p1 = self.get_gene (p1_i, ppop)
        p2 = self.get_gene (p2_i, ppop)
        c1, c2 = p1.crossover (p2, self.random)
        self.set_gene (c1_i, cpop, c1)
        self.set_gene (c2_i, cpop, c2)
    # end def crossover

    def initstring (self, p, pop):
        if not self.randpop:
            self.randpop = self.ramped_half_and_half (self.popsize, 6)
        self.set_gene (p, pop, self.randpop.pop ())
    # end def initstring

    def evaluate (self, p, pop):
        tree = self.get_gene (p, pop)
        r = 0
        for a in (False, True):
            for b in (False, True):
                self.terminals [0].value = a
                self.terminals [1].value = b
                v = tree.eval ()
                if v != ((a and not b) or (not a and  b)):
                    r += 1
        return r
    # end def evaluate

    def print_string (self, file, p, pop):
        print ("%d iterations" % self.GA_iter, file = file)
        file.flush ()
        super ().print_string (file, p, pop)
    # end def print_string

    def stop_cond (self):
        best_eval = self.get_best_report (pga.PGA_OLDPOP, 0)
        if best_eval <= 0:
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Find_XOR

if __name__ == '__main__':
    cmd = ArgumentParser () 
    cmd.add_argument \
        ( '-r', '--random-seed'
        , help    = "Seed for random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    args = cmd.parse_args ()
    xor = Find_XOR (args)
    xor.run ()

