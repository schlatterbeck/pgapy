#!/usr/bin/python

from __future__ import print_function
from pga import PGA, PGA_REPORT_STRING
import sys

class Hello_World (PGA) :

    alphabet  = "abcdefghijklmnopqrstuvwxyz"
    alphabet += alphabet.upper ()
    alphabet += "0123456789!,. "

    def __init__ (self, target) :
        self.target = target
        super (self.__class__, self).__init__ \
            ( type (2), len (target)
            , maximize      = True
            , pop_size      = 100
            , num_replace   = 90
            , mutation_prob = 0.1
            , init          = [(0, len (self.alphabet) - 1)] * len (target)
            , print_options = [PGA_REPORT_STRING]
            )
    # end def __init__

    def evaluate (self, p, pop) :
        r = 0.0
        for i in range (len (self)) :
            idx = self.get_allele (p, pop, i)
            if self.alphabet [idx] == self.target [i] :
                r += 1.0
        return r
    # end def evaluate

    def print_string (self, file, p, pop) :
        s = []
        for i in range (len (self)) :
            idx = self.get_allele (p, pop, i)
            s.append (self.alphabet [idx])
        print (''.join (s), file = file)
        #PGA.print_string (self, file, p, pop)
    # end def print_string

# end class Hello_World

if __name__ == '__main__' :
    target = "Hello World!"
    if len (sys.argv) > 1 :
        target = ' '.join (sys.argv [1:])
    pg = Hello_World (target)
    pg.run ()

