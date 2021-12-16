#!/usr/bin/python3

from __future__ import print_function
from rsclib.autosuper import autosuper
from math import cos, atan, pi
import pga
import sys

def f1 (x1, x2) :
    return x1
# end def f

def f2 (x1, x2) :
    return x2
# end def f

def g1 (x1, x2) :
    return - x1 ** 2 - x2 ** 2 + 1 + 0.1 * cos (16 * atan (x1 / x2))
# end def g1

def g2 (x1, x2) :
    return (x1 - 0.5) ** 2 + (x2 - 0.5) ** 2 - 0.5
# end def g2

class Multi_Objective (pga.PGA, autosuper) :
    """ An example of multi-objective optimization, from NSGA-II paper
        "TNK", see also pgapack nsga-ii directory tnk.c
    """

    def __init__ (self) :
        minmax = ((0, pi), (0, pi))
        self.__super.__init__ \
            ( float, 2
            , maximize             = False
            , random_seed          = 1
            , pop_size             = 100
            , num_replace          = 100
            , select_type          = pga.PGA_SELECT_LINEAR
            , pop_replace_type     = pga.PGA_POPREPL_NSGA_II
            , mutation_only        = True
            , mutation_type        = pga.PGA_MUTATION_DE
            , DE_crossover_prob    = 0.8
            , DE_crossover_type    = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant           = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor      = 0.85
            , init                 = minmax
            , max_GA_iter          = 250
            , num_eval             = 4
            , num_constraint       = 2
            , no_duplicates        = True
            , mutation_bounce_back = True
            , print_options        = [pga.PGA_REPORT_STRING]
            )
    # end def __init__

    def evaluate (self, p, pop) :
        x1 = self.get_allele (p, pop, 0)
        x2 = self.get_allele (p, pop, 1)
        return (f1 (x1, x2), f2 (x1, x2), g1 (x1, x2), g2 (x1, x2))
    # end def evaluate

    def print_string (self, file, p, pop) :
        x1 = self.get_allele (p, pop, 0)
        x2 = self.get_allele (p, pop, 1)
        print \
            ( "f1:%e f2:%e g1:%e g2:%e"
            % (f1 (x1, x2), f2 (x1, x2), g1 (x1, x2), g2 (x1, x2))
            , file = file
            )
        self.__super.print_string (file, p, pop)
    # end def print_string

# end class Multi_Objective

if __name__ == '__main__' :
    print ("Example: Tanaka et. al. (TNK) in Python")
    pg = Multi_Objective ()
    pg.run ()

