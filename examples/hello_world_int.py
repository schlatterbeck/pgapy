#!/usr/bin/python3

from __future__ import print_function
from argparse   import ArgumentParser
import pga
import sys

class Hello_World (pga.PGA):

    alphabet  = "abcdefghijklmnopqrstuvwxyz"
    alphabet += alphabet.upper ()
    alphabet += "0123456789!,. "

    def __init__ (self, args):
        self.args   = args
        self.target = target = args.target
        d = dict \
            ( maximize      = True
            , pop_size      = 100
            , num_replace   = 90
            , mutation_prob = 0.1
            , init          = [(0, len (self.alphabet) - 1)] * len (target)
            , print_options = [pga.PGA_REPORT_STRING]
            , random_seed   = args.random_seed
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super (self.__class__, self).__init__ (type (2), len (target), **d)
    # end def __init__

    def evaluate (self, p, pop):
        r = 0.0
        for i in range (len (self)):
            idx = self.get_allele (p, pop, i)
            if self.alphabet [idx] == self.target [i]:
                r += 1.0
        return r
    # end def evaluate

    def print_string (self, file, p, pop):
        s = []
        for i in range (len (self)):
            idx = self.get_allele (p, pop, i)
            s.append (self.alphabet [idx])
        print (''.join (s), file = file)
        super ().print_string (file, p, pop)
        print ("Evaluations: %d" % self.eval_count, file = file)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval == len (self):
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Hello_World


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
    cmd.add_argument \
        ( "-t", "--target"
        , help    = "Target string to search for, default=%(default)s"
        , default = "Hello World!"
        )
    args = cmd.parse_args (argv)
    pg = Hello_World (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
