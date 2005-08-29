#!/usr/bin/python

from pga import PGA
import pga

class My_PGA (PGA) :
    def evaluate (self, pop, p) :
        count = 0
        for i in range (len (self)) :
            if self.get_allele (pop, p, i) :
                count += 1
        return count
# end class My_PGA

if __name__ == '__main__' :
    pg = My_PGA ("Hallo", type (1), 100, True)
    print pg
    print pg.context
    print pg.__dict__
    pg.run ()

