#!/usr/bin/python3

import pga
import sys
import weakref
from gp import F_nand, Terminal, Genetic_Programming, Function
from random import Random
from argparse import ArgumentParser

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


class Find_Parity_3 (pga.PGA, Genetic_Programming):

    def __init__ (self, args):
        self.args    = args
        self.randpop = []
        self.popsize = 500
        self.random  = PGA_Random (self)
        terms = [Terminal ('a'), Terminal ('b'), Terminal ('c')]
        Genetic_Programming.__init__ (self, [F_nand], terms)
        d = dict \
            ( maximize        = False
            , pop_size        = self.popsize
            , num_replace     = self.popsize - self.popsize // 10
            , tournament_size = 7
            , mutation_prob   = 0.0
            , max_GA_iter     = 50
            , random_seed     = args.random_seed
            , print_options   = [pga.PGA_REPORT_STRING]
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        pga.PGA.__init__ (self, Function, 10, **d)
    # end def __init__

    truthtable = set \
        (( (True,  False, False)
        ,  (False, True,  False)
        ,  (False, False, True)
        ,  (True,  True,  True)
        ))

    def evaluate (self, p, pop):
        tree = self.get_gene (p, pop)
        r = 0
        for a in (False, True):
            for b in (False, True):
                for c in (False, True):
                    self.terminals [0].value = a
                    self.terminals [1].value = b
                    self.terminals [2].value = c
                    v = tree.eval ()
                    rv = 0
                    if (a, b, c) in self.truthtable:
                        rv = 1
                    if v != rv:
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

# end class Find_Parity_3

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( '-r', '-R', '--random-seed'
        , help    = "Seed for random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    args = cmd.parse_args (argv)
    ga = Find_Parity_3 (args)
    #print (ga.mpi_rank, sys.getrefcount (ga))
    ga.run ()
    #print (ga.mpi_rank, sys.getrefcount (ga))
    #ga = None
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
