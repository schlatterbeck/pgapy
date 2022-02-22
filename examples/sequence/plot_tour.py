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

    def edge_weight (self, edge) :
        return self.tsp.get_weight (edge [0] + self.off, edge [1] + self.off)
    # end def edge_weight

    def plot (self, seq, seqoff = 1) :
        l = len (seq)
#        edges = set ()
#        for i in range (l) :
#            j = (i + 1) % l
#            edges.add ((i, j))
#        edges = list (sorted (edges, key = self.edge_weight))
#        tenpc = l // 10
#        edges = set ((seq [e [0]], seq [e [1]]) for e in edges [:tenpc])
        style = '-'
        if self.args.dots :
            style = 'o-'
        X = []
        Y = []
#        fmt = []
#        last_i = None
        for i in range (l) :
#            color = ''
#            if (last_i, i) in edges :
#                color = 'r'
            a = seq [i]
            try :
                x, y = self.tsp.node_coords [a - seqoff + 1]
            except KeyError :
                x, y = self.tsp.display_data [a - seqoff + 1]
            if self.args.exchange_x_y :
                x, y = y, x
            X.append (x)
            Y.append (y)
#            fmt.append (color + style)
            last_i = i
        if not self.args.open :
            X.append (X [0])
            Y.append (Y [0])
#            fmt.append (fmt [0])
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
            s += self.edge_weight ((i1, i2))
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

