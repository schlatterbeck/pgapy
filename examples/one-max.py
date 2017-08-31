#!/usr/bin/python

from __future__ import print_function
from pga import PGA, PGA_REPORT_STRING \
            , PGA_REPORT_ONLINE        \
            , PGA_REPORT_OFFLINE       \
            , PGA_REPORT_WORST         \
            , PGA_REPORT_AVERAGE       \
            , PGA_REPORT_HAMMING
import pga
from operator import add

class One_Max (PGA) :

    def __init__ (self, length) :
        super (self.__class__, self).__init__ \
            ( bool, length
            , maximize      = True
            , print_options = \
                [ PGA_REPORT_STRING
                , PGA_REPORT_WORST
                , PGA_REPORT_AVERAGE
                ]
            )
    def evaluate (self, p, pop) :
        return sum \
            (self.get_allele (p, pop, i) for i in range (len (self)))
    # end def evaluate

    def print_string (self, file, p, pop) :
        r = []
        for k in range (10) :
            r.append (str (int (self.get_allele (p, pop, k))))
        r.append ('..')
        for k in range (90, 100) :
            r.append (str (int (self.get_allele (p, pop, k))))
        print (''.join (r))
    # end def print_string

# end class One_Max

if __name__ == '__main__' :
    pg = One_Max (100)
    pg.run ()

