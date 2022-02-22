#!/usr/bin/python3

import matplotlib.pyplot as plt
import sys
import os
from tsplib95 import load as tspload
from argparse import ArgumentParser

class TSP :
    def __init__ (self, args) :
        self.args = args
        self.tsp  = tspload (args.tsplibfile)
        if self.tsp.tours :
            self.tour = self.tsp
        else :
            if args.tourfile is None :
                bn, ext = os.path.splitext (args.tsplibfile)
                args.tourfile = bn + '.tour'
        if args.tourfile :
            self.tour = tspload (args.tourfile)
        self.off = 0
        try :
            e = self.tsp.get_weight (0, 0)
        except (IndexError, KeyError) :
            self.off = 1
    # end def __init__

    def plot (self, seq, seqoff = 1) :
        X = []
        Y = []
        for i in range (len (seq)) :
            a = seq [i]
            x, y = self.tsp.node_coords [a - seqoff + 1]
            if self.args.exchange_x_y :
                x, y = y, x
            X.append (x)
            Y.append (y)
        if not self.args.open :
            X.append (X [0])
            Y.append (Y [0])
        style = '-'
        if self.args.dots :
            style = 'o-'
        plt.plot (X, Y, style)
        slen = self.seqlen (seq, seqoff = seqoff)
        plt.title ("%s (%s)" % (self.tsp.name, slen))
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
            #print (self.seqlen (t))
            self.plot (t)
    # end def plotall
# end class TSP

if __name__ == '__main__' :
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( 'tsplibfile'
        , help    = 'TSP-Lib compatible file to read as problem'
        )
    cmd.add_argument \
        ( '-d', '--dots'
        , help    = 'Plot dots at nodes'
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-x', '--exchange-x-y'
        , help    = 'Exchange X and Y coordinates (e.g. lin318)'
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-t', '--tourfile'
        , help    = 'TSP-Lib compatible file to read as problem'
        )
    cmd.add_argument \
        ( '-o', '--open'
        , help    = 'Do not plot last (closing) edge'
        , action  = 'store_true'
        )
    args    = cmd.parse_args ()
    tsp     = TSP (args)
    tsp.plotall ()

