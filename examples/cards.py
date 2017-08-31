#!/usr/bin/python

from __future__ import print_function
from pga import PGA, PGA_REPORT_STRING
from operator import mul
import sys
try :
    from functools import reduce
except ImportError:
    pass

class Cards (PGA) :

    def __init__ (self) :
        super (self.__class__, self).__init__ \
            ( bool, 10
            , maximize      = False
            , pop_size      = 30
            , num_replace   = 28
            , print_options = [PGA_REPORT_STRING]
            )
    # end def __init__

    def build_pheno (self, p, pop) :
        g = ([], [])
        for i in range (len (self)) :
            g [self.get_allele (p, pop, i)].append (i + 1)
        return g
    # end def build_pheno

    def evaluate (self, p, pop) :
        g = self.build_pheno (p, pop)
        return \
            ( (1 + 1000 * abs (len (g [0]) - 5))
            * (1 +   10 * abs (sum (g [0]) - 36))
            * (1 +   10 * abs (reduce (mul, g [1], 1) - 360))
            )
    # end def evaluate

    def print_string (self, file, p, pop) :
        g = self.build_pheno (p, pop)
        r = [sum (g [0]), reduce (mul, g [1], 1)]
        s = []
        for n, gg in enumerate (g) :
            s.append (', '.join (str (x) for x in gg) + ': %s' % r [n])
        print (' -- '.join (s))
    # end def print_string

# end class Cards

if __name__ == '__main__' :
    pg = Cards ()
    pg.run ()

