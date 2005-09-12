#!/usr/bin/python

from pga import PGA, PGA_REPORT_STRING \
            , PGA_REPORT_ONLINE        \
            , PGA_REPORT_OFFLINE       \
            , PGA_REPORT_WORST         \
            , PGA_REPORT_AVERAGE       \
            , PGA_REPORT_HAMMING
import pga
from operator import add

class My_PGA (PGA) :
    def evaluate (self, p, pop) :
        return reduce \
            (add, [self.get_allele (p, pop, i) for i in range (len (self))])
    # end def evaluate

    def print_string (self, file, p, pop) :
        file.write ("HUHU\n");
        PGA.print_string (self, file, p, pop)
    # end def print_string

# end class My_PGA

if __name__ == '__main__' :
    pg = My_PGA \
        ( type (2), 100
        , maximize      = True
        , init          = [(0,5)] * 100
        , print_options = \
            [ PGA_REPORT_STRING
            , PGA_REPORT_ONLINE
            , PGA_REPORT_OFFLINE
            , PGA_REPORT_WORST
            , PGA_REPORT_AVERAGE
            ]
        )
    pg.run ()

