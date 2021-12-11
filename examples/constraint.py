#!/usr/bin/python

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

    def __init__ (self) :
        self.__super.__init__ \
            ( float, 2
            , maximize          = False
            , num_eval          = 3
            , pop_size          = 20
            , num_replace       = 20
            , mutation_only     = True
            , max_GA_iter       = 200
            , pop_replace_type  = pga.PGA_POPREPL_PAIRWISE_BEST
            , mutation_type     = pga.PGA_MUTATION_DE
            , DE_crossover_prob = 0.8
            , DE_crossover_type = pga.PGA_DE_CROSSOVER_BIN
            , DE_variant        = pga.PGA_DE_VARIANT_RAND
            , DE_scale_factor   = 0.85
            , random_seed       = 1
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

