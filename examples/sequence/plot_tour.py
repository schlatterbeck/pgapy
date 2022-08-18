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
        if self.tour.name == 'linhp318.tour' :
            self.fixed_edge_len = self.edge_weight ((1, 214))
    # end def __init__

    def edge_weight (self, edge) :
        return self.tsp.get_weight (edge [0] + self.off, edge [1] + self.off)
    # end def edge_weight

    def is_valid_tour (self, seq, seqoff = 1) :
        l = self.tsp.dimension
        d = {}
        if len (seq) != l :
            print ("Invalid length: %s" % len (seq))
            print (seq)
            return False
        for i in range (l) :
            a = seq [i]
            if a <= 0 or a > l :
                print ("Invalid: %s" % (a))
                print (seq)
                return False
            if a in d :
                print ("Dupe: %s" % (a))
                print (seq)
                return False
            d [a] = 1
        return True
    # end def is_valid_tour

    def plot (self, seq, seqoff = 1) :
        l = len (seq)
        # k% Longest edges
        edges = []
        if self.args.percent_longest :
            m = l
            if self.args.open :
                m = l - 1
            for i in range (m) :
                j = (i + 1) % l
                edges.append ((seq [i] - seqoff, seq [j] - seqoff))
            edges = list (sorted (edges, key = lambda e: -self.edge_weight (e)))
            toppc = round (l * (self.args.percent_longest / 100))
            edges = edges [:toppc]
        style = '-'
        if self.args.dots :
            style = 'o-'
        X = []
        Y = []
        for i in range (l) :
            a = seq [i]
            try :
                x, y = self.tsp.node_coords [a - seqoff + 1]
            except KeyError :
                x, y = self.tsp.display_data [a - seqoff + 1]
            if self.args.exchange_x_y :
                x, y = y, x
            X.append (x)
            Y.append (y)
            last_i = i
        if not self.args.open :
            X.append (X [0])
            Y.append (Y [0])
        plt.plot (X, Y, style)
        for e in edges :
            x1, y1 = self.tsp.node_coords [e [0] + 1]
            x2, y2 = self.tsp.node_coords [e [1] + 1]
            if self.args.exchange_x_y :
                x1, y1 = y1, x1
                x2, y2 = y2, x2
            plt.plot ([x1, x2], [y1, y2], 'r' + style)
        slen = self.seqlen (seq, seqoff = seqoff)
        if self.tour.name == 'linhp318.tour' :
            plt.title ( "%s (%s %s)"
                      % (self.tsp.name, slen, slen - self.fixed_edge_len)
                      )
        else :
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
            self.is_valid_tour (t)
            #print (self.seqlen (t))
            self.plot (t)
    # end def plotall

    def seqall (self) :
        for t in self.tour.tours :
            self.is_valid_tour (t)
            yield self.seqlen (t)
    # end def seqall
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
        ( '-l', '--percent-longest'
        , help    = 'Color percent longest edges, default=%(default)s'
        , type    = float
        , default = 0.0
        )
    cmd.add_argument \
        ( '-n', '--noplot'
        , help    = 'Do not plot, only print weights'
        , action  = 'store_false'
        , default = True
        , dest    = 'do_plot'
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
    if args.do_plot :
        tsp.plotall ()
    else :
        for s in tsp.seqall () :
            print ("Len: %s" % s)

