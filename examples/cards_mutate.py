#!/usr/bin/python3

from operator   import mul
from argparse   import ArgumentParser
from functools  import reduce
import pga
import sys

class Cards (pga.PGA):

    def __init__ (self, args):
        self.args = args
        d = dict \
            ( maximize      = False
            , pop_size      = 20
            , num_replace   = 19
            , mutation_only = True
            , mutation_prob = 0.2
            , max_GA_iter   = 1000
            , print_options = [pga.PGA_REPORT_STRING]
            , random_seed   = self.args.random_seed
            )
        if self.args.output_file:
            d ['output_file'] = self.args.output_file
        super (self.__class__, self).__init__ (int, 10, **d)
    # end def __init__

    def build_pheno (self, p, pop):
        g = []
        for i in range (len (self)):
            g.append (self.get_allele (p, pop, i) + 1)
        return g
    # end def build_pheno

    def evaluate (self, p, pop):
        g = self.build_pheno (p, pop)
        return \
            ( abs (sum (g [0:5]) - 36) * 10
            + abs (reduce (mul, g [5:10], 1) - 360)
            )
    # end def evaluate

    def print_string (self, file, p, pop):
        g = self.build_pheno (p, pop)
        r = [sum (g [0:5]), reduce (mul, g [5:10], 1)]
        s = []
        for n, rr in enumerate (r):
            s1 = ', '.join (str (x) for x in g [5*n: 5*(n+1)])
            s.append ('%s: %s' % (s1, r [n]))
        print (' -- '.join (s), file = file)
        print ("%d iterations" % self.GA_iter, file = file)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval <= 0:
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Cards

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( "-R", "--random-seed"
        , help    = "Seed random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    args = cmd.parse_args (argv)
    pg = Cards (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
