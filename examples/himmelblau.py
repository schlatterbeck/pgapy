#!/usr/bin/python3

from __future__ import print_function
from argparse         import ArgumentParser
import sys
import pga

class Himmelblau (pga.PGA):
    """ Constrained Himmelblau's function
        Example 4.3.1 p.157 from Kalyanmoy Deb. Optimization for
        Engineering Design – Algorithms and Examples. PHI Learning,
        Delhi, India, second edition, 2012.
        The optimum according to Deb is at (0.829, 2.933). The result we
        find is about the same (due to rounding). The function value we
        get at the optimum point is 60.37361
    """

    def __init__ (self, args):
        self.args = args
        minmax = ((0, 30), (0, 30))
        d = dict \
            ( maximize             = False
            , random_seed          = args.random_seed
            , pop_size             = 6
            , num_replace          = 6
            , select_type          = pga.PGA_SELECT_LINEAR
            , pop_replace_type     = pga.PGA_POPREPL_PAIRWISE_BEST
            , mutation_only        = True
            , mutation_type        = pga.PGA_MUTATION_DE
            , mutation_bounce_back = True
            , DE_crossover_prob    = 0.8
            , DE_crossover_type    = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant           = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor      = 0.85
            , init                 = minmax
            , max_GA_iter          = 500
            , num_eval             = 2
            , num_constraint       = 1
            , print_options        = [pga.PGA_REPORT_STRING]
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (float, 2, **d)
    # end def __init__

    def evaluate (self, p, pop):
        x1 = self.get_allele (p, pop, 0)
        x2 = self.get_allele (p, pop, 1)
        f = (x1 ** 2 + x2 - 11) ** 2 + (x1 + x2 ** 2 - 7) ** 2
        g = 26 - (x1 - 5) ** 2 - x2 ** 2
        return f, g
    # end def evaluate

    def print_string (self, file, p, pop):
        f, g = self.evaluate (p, pop)
        #print (file.fileno ())
        print ("f:%e g:%e" % (f, g), file = file)
        print ("Evaluations: %d" % self.eval_count, file = file)
        super ().print_string (file, p, pop)
    # end def print_string

# end class Himmelblau

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
    pg = Himmelblau (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
