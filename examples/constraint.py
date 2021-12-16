#!/usr/bin/python3

from __future__ import print_function
from rsclib.autosuper import autosuper
import pga
import sys

def f (x1, x2) :
    return (x1 - 0.8) ** 2 + (x2 - 0.3) ** 2
# end def f

def g1 (x1, x2) :
    return ((x1 - 0.2) ** 2 + (x2 - 0.5) ** 2) / 0.16 - 1
# end def g1

def g2 (x1, x2) :
    return 1 - ((x1 + 0.5) ** 2 + (x2 - 0.5) ** 2) / 0.81
# end def g2

class Constrained (pga.PGA, autosuper) :
    """ First example from Deb 2000, two constraints, single objective
        function.
    """

    def __init__ (self) :
        minmax = ((-5, 5), (-5, 5))
        self.__super.__init__ \
            ( float, 2
            , maximize          = False
            , random_seed       = 1
            , pop_size          = 60
            , num_replace       = 60
            , select_type       = pga.PGA_SELECT_LINEAR
            , pop_replace_type  = pga.PGA_POPREPL_PAIRWISE_BEST
            , mutation_only     = True
            , mutation_type     = pga.PGA_MUTATION_DE
            , DE_crossover_prob = 0.8
            , DE_crossover_type = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant        = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor   = 0.85
            , init              = minmax
            , max_GA_iter       = 100
            , num_eval          = 3
            , print_options     = [pga.PGA_REPORT_STRING]
            )
    # end def __init__

    def evaluate (self, p, pop) :
        x1 = self.get_allele (p, pop, 0)
        x2 = self.get_allele (p, pop, 1)
        return (f (x1, x2), g1 (x1, x2), g2 (x1, x2))
    # end def evaluate

    def print_string (self, file, p, pop) :
        x1 = self.get_allele (p, pop, 0)
        x2 = self.get_allele (p, pop, 1)
        print \
            ( "f:%e g1:%e g2:%e"
            % (f (x1, x2), g1 (x1, x2), g2 (x1, x2))
            , file = file
            )
        self.__super.print_string (file, p, pop)
    # end def print_string

# end class Constrained

if __name__ == '__main__' :
    pg = Constrained ()
    pg.run ()

