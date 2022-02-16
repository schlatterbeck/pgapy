#!/usr/bin/python3

import matplotlib.pyplot as plt
import sys
import os
from tsplib95 import load as tspload

class TSP :
    def __init__ (self, filename) :
        self.tsp  = tspload (filename)
        bn, ext   = os.path.splitext (filename)
        self.tour = tspload (bn + '.tour')
        self.off = 0
        try :
            e = self.tsp.get_weight (0, 0)
        except (IndexError, KeyError) :
            self.off = 1
    # end def __init__

    def plot (self, seq, seqoff = 1, close = 1) :
        X = []
        Y = []
        for i in range (len (seq)) :
            a = seq [i]
            node = self.tsp.node_coords [a - seqoff + 1]
            #print (node)
            X.append (node [1])
            Y.append (node [0])
        if close :
            X.append (X [0])
            Y.append (Y [0])
        plt.plot (X, Y, 'bo-')
        plt.show ()
    # end def plot

    def seqlen (self, seq, seqoff = 1) :
        s = 0
        for i in range (len (seq)) :
            i1 = seq [i] - seqoff
            i2 = seq [(i+1) % self.tsp.dimension] - seqoff
            s += self.tsp.get_weight (i1 + self.off, i2 + self.off)
        return s
    # end def seqlen

    def plotall (self) :
        for t in self.tour.tours :
            print (self.seqlen (t))
            self.plot (t)
    # end def plotall
# end class TSP

if __name__ == '__main__' :
    tsp = TSP (sys.argv [1])
    tsp.plotall ()
