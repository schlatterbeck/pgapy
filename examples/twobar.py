#!/usr/bin/python3

from argparse import ArgumentParser
from math import sqrt
import pga
import sys

class Two_Bar (pga.PGA):
    """ Example from [1]
    [1] Tapabrata Ray, Kang Tai, and Kin Chye Seow.  Multiobjective
        design optimization by an evolutionary algorithm. Engineering
        Optimization, 33(4):399–424, 2001.
    """

    def __init__ (self, args):
        self.args   = args
        minmax = [(0.1, 2.25), (0.5, 2.5)]
        d = dict \
            ( maximize             = False
            , num_eval             = 4
            , num_constraint       = 2
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
        super ().__init__ (float, 2, **d)
    # end def __init__

    def evaluate (self, p, pop):
        rho     = 0.283
        h       = 100
        P       = 1e4
        E       = 3e7
        sigma_0 = 2e4
        A_min   = 1
        x = []
        for i in range (len (self)):
            x.append (self.get_allele (p, pop, i))
        f1  = 2 * rho * h * x [1] * sqrt (1 + x [0] ** 2)
        f2  = P * h * (1 + x [0] ** 2) ** 1.5 * (1 + x [0] ** 4) ** 0.5 \
            / (2 * sqrt (2) * E * x [0] ** 2 * x [1])
        g1  = P * (1 + x [0]) * (1 + x [0] ** 2) ** 0.5 \
            / (2 * sqrt (2) * x [0] * x [1]) \
            - sigma_0
        g2  = P * (1 - x [0]) * (1 + x [0] ** 2) ** 0.5 \
            / (2 * sqrt (2) * x [0] * x [1]) \
            - sigma_0
        return f1, f2, g1, g2
    # end def evaluate

# end class Two_Bar

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
    pg = Two_Bar (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
