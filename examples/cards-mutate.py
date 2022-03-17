#!/usr/bin/python3

from __future__ import print_function
from operator import mul
import pga
import sys
try :
    from functools import reduce
except ImportError:
    pass

class Cards (pga.PGA) :

    def __init__ (self) :
        super (self.__class__, self).__init__ \
            ( int, 10
            , maximize      = False
            , pop_size      = 20
            , num_replace   = 19
            , mutation_only = True
            , mutation_prob = 0.2
            , max_GA_iter   = 1000
            #, random_seed   = 1
            , print_options = [pga.PGA_REPORT_STRING]
            )
    # end def __init__

    def build_pheno (self, p, pop) :
        g = []
        for i in range (len (self)) :
            g.append (self.get_allele (p, pop, i) + 1)
        return g
    # end def build_pheno

    def evaluate (self, p, pop) :
        g = self.build_pheno (p, pop)
        return \
            ( abs (sum (g [0:5]) - 36) * 10
            + abs (reduce (mul, g [5:10], 1) - 360)
            )
    # end def evaluate

    def print_string (self, file, p, pop) :
        g = self.build_pheno (p, pop)
        r = [sum (g [0:5]), reduce (mul, g [5:10], 1)]
        s = []
        for n, rr in enumerate (r) :
            s1 = ', '.join (str (x) for x in g [5*n: 5*(n+1)])
            s.append ('%s: %s' % (s1, r [n]))
        print (' -- '.join (s), file = file)
        print ("%d iterations" % self.GA_iter)
    # end def print_string

    def stop_cond (self) :
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval <= 0 :
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Cards

if __name__ == '__main__' :
    pg = Cards ()
    pg.run ()

