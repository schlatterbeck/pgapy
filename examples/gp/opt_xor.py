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
        self.popsize = 100
        self.random  = PGA_Random (self)
        terms = [Terminal ('a'), Terminal ('b')]
        Genetic_Programming.__init__ (self, [F_nand], terms, debug = args.debug)
        d = dict \
            ( maximize      = False
            , pop_size      = self.popsize
            , num_replace   = self.popsize - self.popsize // 10
            , mutation_prob = 0.0
            , max_GA_iter   = 50
            , random_seed   = args.random_seed
            , print_options = [pga.PGA_REPORT_STRING]
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        pga.PGA.__init__ (self, Function, 10, **d)
    # end def __init__

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
        if self.args.verbose:
            print (self.get_gene (p, pop).as_dot (), file = file)
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

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-D", "--debug"
        , help    = "Turn on debugging assertions"
        , action  = 'store_true'
        )
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
    cmd.add_argument \
        ( "-v", "--verbose"
        , help    = "Print graphviz output"
        , action  = 'store_true'
        )
    args = cmd.parse_args (argv)
    xor = Find_XOR (args)
    xor.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
