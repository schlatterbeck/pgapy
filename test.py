#!/usr/bin/python

from pga import PGA
import pga

class My_PGA (PGA) :
    def evaluate (self, pop, p) :
        count = 0
        for i in range (len (self)) :
            if self.get_allele (pop, p) :
                count += 1
        print "Count:", count
        return count
# end class My_PGA

pg = My_PGA (str ("Hallo"), type (1), 100)
pg.run ()

