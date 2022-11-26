#!/usr/bin/python3

from __future__ import print_function, division
from argparse import ArgumentParser
from math import pi
import pga
import sys

class Vibr (pga.PGA):
    """ Example from [1]
    [1] Tapabrata Ray, Kang Tai, and Kin Chye Seow.  Multiobjective
        design optimization by an evolutionary algorithm. Engineering
        Optimization, 33(4):399–424, 2001.
    """
    material = \
        [ (2770, 70e9, 1500), (100, 1.6e9, 500), (7780, 200e9, 800) ]

    def __init__ (self, args):
        self.args = args
        minmax = [ (3, 6), (.35, .5), (.01, .6), (.01, .6), (.01, .6)
                 , (0, 2.99), (0, 2.99), (0, 2.99)
                 ]

        d = dict \
            ( maximize             = False
            , num_eval             = 7
            , num_constraint       = 5
            , pop_size             = 40
            , num_replace          = 40
            , select_type          = pga.PGA_SELECT_LINEAR
            , pop_replace_type     = pga.PGA_POPREPL_NSGA_II
            , mutation_only        = True
            , mutation_type        = pga.PGA_MUTATION_DE
            , DE_crossover_prob    = 0.9
            , DE_crossover_type    = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant           = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor      = 0.85
            , DE_jitter            = 0.001
            , DE_dither            = 0.4
            , init                 = minmax
            , max_GA_iter          = 200
            , print_options        = [pga.PGA_REPORT_STRING]
            , mutation_bounce_back = True
            )
        if args.random_seed:
            d ['random_seed'] = args.random_seed
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (float, 8, **d)
    # end def __init__

    def evaluate (self, p, pop):
        x = []
        for i in range (len (self)):
            x.append (self.get_allele (p, pop, i))
        L, b, d1, d2, d3, m1, m2, m3 = x
        m1    = int (m1)
        m2    = int (m2)
        m3    = int (m3)
        rhom1 = self.material [m1][0]
        rhom2 = self.material [m2][0]
        rhom3 = self.material [m3][0]
        E1    = self.material [m1][1]
        E2    = self.material [m2][1]
        E3    = self.material [m3][1]
        c1    = self.material [m1][2]
        c2    = self.material [m2][2]
        c3    = self.material [m3][2]

        mu  = 2 * b * (rhom1 * d1 + rhom2 * (d2 - d1) + rhom3 * (d3 - d2))
        E   = 2 * b / 3 \
            * ( E1 * d1 ** 3
              + E2 * (d2 ** 3 - d1 ** 3)
              + E3 * (d3 ** 3 - d2 ** 3)
              )
        g1  = mu * L - 2800
        g2  = 1e-5 - d2 + d1
        g3  = 1e-5 - d3 + d2
        g4  = d2 - d1 - 0.01
        g5  = d3 - d3 - 0.01
        if g1 > 0 or g2 > 0 or g3 > 0 or g4 > 0 or g5 > 0:
            f1 = f2 = 0.0
        else:
            f1  = (pi / (2 * L ** 2)) * (E / mu) ** 0.5
            f2  = 2 * b * (c1 * d1 + c2 * (d2 - d1) + c3 * (d3 - d2))
        return -f1, f2, g1, g2, g3, g4, g5
    # end def evaluate

# end class Vibr

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
    pg = Vibr (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
