#!/usr/bin/python3

from __future__ import print_function
from rsclib.autosuper import autosuper
import pga

class Himmelblau (pga.PGA, autosuper) :
    """ Constrained Himmelblau's function
        Example 4.3.1 p.157 from Kalyanmoy Deb. Optimization for
        Engineering Design – Algorithms and Examples. PHI Learning,
        Delhi, India, second edition, 2012.
    """

    def __init__ (self) :
        minmax = ((0, 30), (0, 30))
        self.__super.__init__ \
            ( float, 2
            , maximize             = False
            , random_seed          = 1
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
    # end def __init__

    def evaluate (self, p, pop) :
        x1 = self.get_allele (p, pop, 0)
        x2 = self.get_allele (p, pop, 1)
        f = (x1 ** 2 + x2 - 11) ** 2 + (x1 + x2 ** 2 - 7) ** 2
        g = 26 - (x1 - 5) ** 2 - x2 ** 2
        return f, g
    # end def evaluate

    def print_string (self, file, p, pop) :
        f, g = self.evaluate (p, pop)
        print ("f:%e g:%e" % (f, g), file = file)
        print ("Evaluations: %d" % self.eval_count)
        self.__super.print_string (file, p, pop)
    # end def print_string

# end class Himmelblau

if __name__ == '__main__' :
    pg = Himmelblau ()
    pg.run ()

