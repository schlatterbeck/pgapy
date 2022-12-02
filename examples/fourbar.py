#!/usr/bin/python3

from __future__       import print_function, division
from argparse         import ArgumentParser
from math             import sqrt
import pga
import sys

class Four_Bar (pga.PGA):
    """ Example from [1]
    [1] Tapabrata Ray, Kang Tai, and Kin Chye Seow.  Multiobjective
        design optimization by an evolutionary algorithm. Engineering
        Optimization, 33(4):399–424, 2001.
    """

    def __init__ (self, args):
        self.args  = args
        self.F     = 10.
        self.E     = 2e5
        self.L     = 200.
        self.sigma = 10.
        q2         = self.q2 = sqrt (2)
        b = self.F / self.sigma
        minmax = [(b, 3*b), (q2 * b, 3 * b), (q2 * b, 3 * b), (b, 3 * b)]
        d = dict \
            ( maximize             = False
            , num_eval             = 2
            , num_constraint       = 0
            , pop_size             = 60
            , num_replace          = 60
            , select_type          = pga.PGA_SELECT_LINEAR
            , pop_replace_type     = pga.PGA_POPREPL_NSGA_II
            , mutation_only        = True
            , mutation_type        = pga.PGA_MUTATION_DE
            , DE_crossover_prob    = 0.8
            , DE_crossover_type    = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant           = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor      = 0.85
            , init                 = minmax
            , max_GA_iter          = 100
            , print_options        = [pga.PGA_REPORT_STRING]
            , mutation_bounce_back = True
            )
        if args.random_seed:
            d ['random_seed'] = args.random_seed
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (float, 4, **d)
    # end def __init__

    def evaluate (self, p, pop):
        x = []
        for i in range (len (self)):
            x.append (self.get_allele (p, pop, i))
        q2  = self.q2
        f1  = self.L * (2 * x [0] + q2 * x [1] + q2 * x [2] + x [3])
        f2  = self.F * self.L / self.E \
            * (2 / x [0] + 2 * q2 / x [1] - 2 * q2 / x [2] + 2 / x [3])
        return f1, f2
    # end def evaluate

# end class Four_Bar

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( '-r', '-R', '--random-seed'
        , type    = int
        )
    args = cmd.parse_args (argv)
    pg = Four_Bar (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
