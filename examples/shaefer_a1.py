#!/usr/bin/python3

from argparse import ArgumentParser
import sys
import pga
import numpy as np

class Shaefer_A1 (pga.PGA):
    """ Example function to be optimized in
        Craig G. Shaefer. The ARGOT strategy: Adaptive representation
        genetic optimizer technique. In John J. Grefenstette, editor,
        Proceedings of the Second International Conference on Genetic
        Algorithms and Their Applications (ICGA), page 50–58. Lawrence
        Erlbaum Associates, Cambridge, MA, July 1987.
        That algorithm is not very interesting today, as can be seen in
        this implementation using Differential Evolution finds a
        solution quite fast. Unfortunately the original paper doesn't
        give the population size used.
    """

    def __init__ (self, args):
        self.args = args
        minmax = ((-25, 25), (-25, 25), (-25, 25), (-25, 25))
        d = dict \
            ( maximize             = False
            , random_seed          = args.random_seed
            , pop_size             = args.pop_size
            , num_replace          = args.pop_size
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
            , max_GA_iter          = args.generations
            , num_eval             = 1
            , num_constraint       = 0
            , print_options        = [pga.PGA_REPORT_STRING]
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (float, len (minmax), **d)
    # end def __init__

    @staticmethod
    def func (x5, x6, x10, x13, x_range = 10, n_data = 50):
        """
        >>> f = Shaefer_A1.func
        >>> print ("%.6f" % f (0.5, 1, 0.75, 0))
        0.000000
        >>> print ("%.6f" % f (-0.5, -1, 0.75, 0))
        0.000000
        >>> print ("%.4f" % f (.49449, .99511, .75091, 6.0997e-3))
        1164.1268
        >>> print ("%.4f" % f (-1.879062, -0.005072817, -24.99929, 0.05054976))
        220125.9556
        """
        s = 0.0
        r = x_range
        step = 2 * r / n_data
        a = np.arange (-r, r + step, step)
        s = np.sum \
            ( ( x13 ** 2 * a ** 4
              + 2 * x6 * x13 * a ** 3
              + (x6 ** 2 + 2 * x5 * x13) * a ** 2
              + 2 * x5 * x6 * a
              + x5 ** 2
              + x10
              - a ** 2 - a - 1
              ) ** 2
            )
        return s
    # end def func

    def evaluate (self, p, pop):
        x5  = self.get_allele (p, pop, 0)
        x6  = self.get_allele (p, pop, 1)
        x10 = self.get_allele (p, pop, 2)
        x13 = self.get_allele (p, pop, 3)
        r = self.args.x_range
        return self.func (x5, x6, x10, x13, r, self.args.n_data)
    # end def evaluate

    def print_string (self, file, p, pop):
        s = self.evaluate (p, pop)
        print ("s:%e" % s, file = file)
        print ("Evaluations: %d" % self.eval_count, file = file)
        super ().print_string (file, p, pop)
    # end def print_string

# end class Shaefer_A1

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( "-g", "--generations"
        , help    = "Number of generations, default=%(default)s"
        , type    = int
        , default = 500
        )
    cmd.add_argument \
        ( "-k", "--n-data"
        , help    = "Number of data points, default=%(default)s"
        , type    = int
        , default = 50
        )
    cmd.add_argument \
        ( "-p", "--pop-size"
        , help    = "Population size, default=%(default)s"
        , type    = int
        , default = 50
        )
    cmd.add_argument \
        ( "-R", "--random-seed"
        , help    = "Seed random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    cmd.add_argument \
        ( "-x", "--x-range"
        , help    = "X-range -x..x, default=%(default)s"
        , type    = float
        , default = 10
        )
    args = cmd.parse_args (argv)
    pg = Shaefer_A1 (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
