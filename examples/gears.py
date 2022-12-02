#!/usr/bin/python3

from __future__       import print_function, division
from argparse         import ArgumentParser
from math             import gcd
import pga
import sys

class Gears (pga.PGA):
    """ Example from presentation by Deb 2008
    """

    def __init__ (self, args):
        self.args   = args
        self.factor = self.args.numerator / self.args.denominator
        minmax = [(self.args.min_tooth, self.args.max_tooth)] * 4
        d = dict \
            ( maximize             = False
            , num_eval             = 2
            , num_constraint       = 1
            , pop_size             = 60
            , num_replace          = 60
            , select_type          = pga.PGA_SELECT_LINEAR
            , pop_replace_type     = pga.PGA_POPREPL_PAIRWISE_BEST
            , mutation_only        = True
            , mutation_type        = pga.PGA_MUTATION_DE
            , DE_crossover_prob    = 0.8
            , DE_crossover_type    = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant           = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor      = 0.85
            , init                 = minmax
            , max_GA_iter          = 1000
            , print_options        = [pga.PGA_REPORT_STRING]
            , mutation_bounce_back = True
            )
        if args.random_seed:
            d ['random_seed'] = args.random_seed
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (float, 4, **d)
    # end def __init__

    def err (self, x1, x2, x3, x4):
        return abs (self.factor - (x3 * x4) / (x1 * x2)) / self.factor * 100
    # end def err

    def evaluate (self, p, pop):
        x = []
        for i in range (4):
            x.append (round (self.get_allele (p, pop, i)))
        gc1 = gcd (x [0], x [2]) + gcd (x [1], x [3])
        gc2 = gcd (x [0], x [3]) + gcd (x [1], x [2])
        f   = (1 / self.factor - (x [0] * x [1]) / (x [2] * x [3])) ** 2
        return f, min (gc1, gc2) - 2
    # end def evaluate

    def print_string (self, file, p, pop):
        x = []
        for i in range (4):
            x.append (round (self.get_allele (p, pop, i)))
        print (x, file = file)
        print ("Gear Error: %12.9f%%" % self.err (*x), file = file)
        print ("Random seed: %d" % self.random_seed, file = file)
        super ().print_string (file, p, pop)
    # end def print_string

# end class Gears

def main (argv):
    """ A nice problem is
        -l 17 -u 90 -n 950 -d 150
        With a solution with GCD violation but perfect match:
        17,18,38,51
        And several good solutions without GCD in factors, the first is
        with constraints on the tooth lower bound to 12, upper bound 60
        20,17,49,44 0.123839009%
        23,49,83,86 0.004670060%
        31,34,75,89 0.004993508%
        37,28,81,81 0.005080526%
        43,23,72,87 0.005080268%
        29,23,66,64 0.007890791%
        31,31,78,78 0.038337258%
        The last is nice from a manufacturing perspective because there
        are only two different gears. Found by hand not by program ;-)
        Generally if we only want two gears we can use:
        >>> from fractions import Fraction
        >>> from math import sqrt
        >>> x = sqrt (950/150)
        >>> Fraction (x).limit_denominator (50)
        Fraction(78, 31)
    """
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( '-c', '--check'
        , help    = "A comma-separated list of 4 integers for gears to check"
        )
    cmd.add_argument \
        ( '-l', '--min-tooth'
        , type    = int
        , default = 12
        )
    cmd.add_argument \
        ( '-u', '--max-tooth'
        , type    = int
        , default = 60
        )
    cmd.add_argument \
        ( '-d', '--denominator'
        , type    = float
        , default = 1.0
        )
    cmd.add_argument \
        ( '-n', '--numerator'
        , type    = float
        , default = 6.931
        )
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( '-r', '-R', '--random-seed'
        , type    = int
        )
    args = cmd.parse_args (argv)
    pg = Gears (args)
    if args.check:
        x = [int (i) for i in args.check.split (',')]
        print ("Factor: %f" % pg.factor)
        print ("Gear Error: %12.9f%%" % pg.err (*x))
    else:
        pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
