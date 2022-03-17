#!/usr/bin/python3

from __future__ import print_function
from pga import PGA, PGA_REPORT_STRING
from argparse import ArgumentParser
import sys

class Sorter (PGA) :

    def __init__ (self, length) :
        super (self.__class__, self).__init__ \
            ( int, length
            , maximize      = False
            , pop_size      = 30
            , num_replace   = 28
            , mutation_prob = 0.3
            , init          = [(0, 10 * length - 1)] * length
            , print_options = [PGA_REPORT_STRING]
            )
    # end def __init__

    def evaluate (self, p, pop) :
        r = 0.0
        last = None
        violation = 1.0
        distance  = 0.0
        for i in range (1, len (self)) :
            last = self.get_allele (p, pop, i - 1)
            this = self.get_allele (p, pop, i)
            if last >= this :
                violation *= 100
                distance  += (last + 1 - this)
            else :
                distance  += this - last
        return distance * violation
    # end def evaluate

    def print_string (self, file, p, pop) :
        print \
            (', '.join \
                (str (self.get_allele (p, pop, x)) for x in range (len (self)))
            , file = file
            )
    # end def print_string

# end class Sorter

if __name__ == '__main__' :
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-l", "--length"
        , help    = "Length of string to sort"
        , type    = int
        , default = 10
        )
    args = cmd.parse_args ()
    pg = Sorter (args.length)
    pg.run ()

