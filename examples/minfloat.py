#!/usr/bin/python

from __future__ import print_function
from rsclib.autosuper import autosuper
import pga
import sys

class Minfloat (pga.PGA, autosuper) :
    """ This demonstrates the use of the new Differential Evolution
        strategy. Note that there is no selection pressure in the linear
        selection method. Instead selection pressure is employed with
        the population replacement method PGA_POPREPL_PAIRWISE_BEST.
        If a pop_replace_type without selection pressure *and* a
        selection method without selection pressure is used, the
        search degenerates to a random walk.
    """

    def __init__ (self) :
        self.gene_distance = self.euclidian_distance
        popsize = 30
        super (self.__class__, self).__init__ \
            ( float, 9
            , maximize                    = False
            , pop_size                    = popsize
            , num_replace                 = popsize
            , print_options               = [pga.PGA_REPORT_STRING]
            , init                        = [(-10.0, 10.0)] * 9
            , select_type                 = pga.PGA_SELECT_LINEAR
            , pop_replace_type            = pga.PGA_POPREPL_PAIRWISE_BEST
            # Either use mutation_bounded, which sets a value that
            # violates the bounds to the bound or mutation_bounce_back
            # which sets a value that violates the bound to a random
            # value between the old value and the bound. The latter
            # method is said to retain better population variance.
            #, mutation_bounded            = True
            , mutation_bounce_back        = True
            , mutation_only               = True
            , mutation_type               = pga.PGA_MUTATION_DE
            , tournament_with_replacement = False
            , random_seed                 = 42
            , DE_variant                  = pga.PGA_DE_VARIANT_BEST
            , DE_crossover_prob           = 0.2
            , DE_jitter                   = 0.001
            , DE_scale_factor             = 0.85 - (popsize * 0.0005)
            )
        self.from_stop_cond = False
        self.num_evals = 0
    # end def __init__

    def evaluate (self, p, pop) :
        assert self.from_stop_cond or not self.get_evaluation_up_to_date (p, pop)
        self.num_evals += 1
        return sum (self.get_allele (p, pop, k) for k in range (len (self)))
    # end def evaluate

    def stop_cond (self) :
        """ Stop when the evaluation has reached 0.99 * 10 * len (self)
        """
        best = self.get_best_index (pga.PGA_OLDPOP)
        self.from_stop_cond = True
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        self.from_stop_cond = False
        if eval <= -0.99999 * 10 * len (self) :
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

    def print_string (self, file, p, pop) :
        print ("evals: %s" % self.num_evals, file = file)
        print ("best index: %d" % self.get_best_index (pop), file = file);
        file.flush ()
        self.__super.print_string (file, p, pop)
    # end def print_string

# end class Minfloat

if __name__ == '__main__' :
    pg = Minfloat ()
    pg.run ()
