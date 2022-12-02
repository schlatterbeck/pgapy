#!/usr/bin/python3

from __future__       import print_function, division
from argparse         import ArgumentParser
from math             import gcd
import numpy as np
import pga
import sys

class DTLZ2 (pga.PGA):

    def __init__ (self, args):
        self.args = args
        self.nobj = args.nobjective
        self.dim  = args.dimension
        refpoints = pga.das_dennis (self.nobj, self.args.das_dennis_partitions)
        l         = len (refpoints)
        l         = (l + 3) // 4 * 4
        self.k    = self.dim - self.nobj + 1
        minmax    = [[0.0, 1.0]] * self.dim
        d = dict \
            ( maximize             = False
            , num_eval             = self.nobj
            , num_constraint       = 0
            , num_replace          = l
            , select_type          = pga.PGA_SELECT_LINEAR
            , pop_replace_type     = pga.PGA_POPREPL_NSGA_III
            , mutation_only        = True
            , mutation_type        = pga.PGA_MUTATION_DE
            , DE_crossover_prob    = 0.0
            , DE_crossover_type    = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant           = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor      = 0.40
            , DE_jitter            = 0.30
            , init                 = minmax
            , max_GA_iter          = 250
            , print_options        = [pga.PGA_REPORT_STRING]
            , mutation_bounce_back = True
            , no_duplicates        = True
            #, print_frequency      = 1
            )
        if args.reference_directions:
            l = 96
            d.update \
                ( refdir_scale         = 0.05
                , refdir_partitions    = 8
                , reference_directions = [[1, 1, 1], [1, 1, 1]]
                , num_replace          = l
                )
        else:
            d.update (reference_points = refpoints)
        if args.random_seed:
            d ['random_seed'] = args.random_seed
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (float, self.dim, **d)
        assert l == self.pop_size
    # end def __init__

    def evaluate (self, p, pop):
        x = []
        y = []
        for i in range (self.dim):
            x.append (self.get_allele (p, pop, i))
        x = np.array (x)
        g = sum ((x [-self.k:] - 0.5) ** 2)
        for i in range (self.nobj):
            p = 1.0
            for j in range (self.nobj - i - 1):
                p *= np.cos (x [j] * np.pi / 2)
            if i > 0:
                p *= np.sin (x [self.nobj - i - 1] * np.pi / 2)
            scale = 1
            if self.args.scale:
                scale = self.args.scale ** i
            y.append (p * (1 + g) * scale)
        return y
    # end def evaluate

# end class DTLZ2

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( '-d', '--dimension'
        , help    = "Dimension of problem, default=%(default)s"
        , type    = int
        , default = 12
        )
    cmd.add_argument \
        ( '-m', '--nobjective'
        , help    = "Number of objectives, default=%(default)s"
        , type    = int
        , default = 3
        )
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( '-p', '--das-dennis-partitions'
        , help    = "Number of das dennis partitions, default=%(default)s"
        , type    = int
        , default = 12
        )
    cmd.add_argument \
        ( '-r', '--reference-directions'
        , help   = 'Use reference directions instead of reference points'
        , action = 'store_true'
        )
    cmd.add_argument \
        ( '-s', '--scale'
        , help    = "Scaling of problem, default=%(default)s"
        , type    = float
        , default = 1
        )
    cmd.add_argument \
        ( '-R', '--random-seed'
        , type    = int
        )
    args = cmd.parse_args (argv)
    pg = DTLZ2 (args)
    pg.run ()

if __name__ == '__main__':
    main (sys.argv [1:])
