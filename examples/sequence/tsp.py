#!/usr/bin/python3

import sys
import pga
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from tsplib95 import load as tspload
from random   import Random
from copy     import copy, deepcopy
from bisect   import bisect_left, bisect_right
from collections import OrderedDict

class PGA_Random (Random) :

    def __init__ (self, pga_instance) :
        self.pga_instance = pga_instance
        super (self.__class__, self).__init__ (1)
    # end def __init__

    def getstate (self, *args, **kw) :
        raise NotImplementedError ("Getting the state is not possible")
    # end def getstate

    def setstate (self, *args, **kw) :
        raise NotImplementedError ("Setting the state is not possible")
    # end def setstate

    def random (self) :
        return self.pga_instance.random01 ()
    # end def random

# end class PGA_Random

class Long_Edge_Iter :
    def __init__ (self, parent, allele) :
        self.parent   = parent
        self.edges    = parent.long_edges (allele)
        self.allele   = allele
        self.edge_idx = list (range (len (self.edges)))
        self.idx      = 0
        self.parent.random.shuffle (self.edge_idx)
    # end def __init__

    @property
    def n (self) :
        return (self.e [1][0] - self.e [0][0]) % len (self.parent)
    # end def n

    def _next (self) :
        if self.idx > len (self.edge_idx) - 3 :
            self.parent.random.shuffle (self.edge_idx)
            self.idx = 0
        r = self.edge_idx [self.idx : self.idx + 3]
        self.idx += 3
        self.e = [self.edges [i] for i in r]
        return self.e
    # end def _next

    def is_valid (self) :
        return self.parent.valid_index \
            (self.e [2][0], self.allele, self.e [0][1], self.n)
    # end def is_valid

    def next (self) :
        l = len (self.parent)
        self._next ()
        while not self.is_valid () :
            self._next ()
        return self.e
    # end def next

# end class Long_Edge_Iter

class Segment_Pointer :
    """ Pointer (currently index) to next/prev segment
        At some point we may want to modify this to store a ref to the
        next segment.
    """

    def __init__ (self, tail, dir) :
        self.tail = tail
        self.dir  = dir
        self.idx  = None
        self.ptr  = None
    # end def __init__

    def set_ptr (self, ptr, idx) :
        assert self.ptr is None
        self.ptr = ptr
        self.idx = idx
    # end def set_ptr

    def __str__ (self) :
        p = '' if self.ptr is None else self.ptr
        if self.dir > 0 :
            return "%d>%s" % (self.tail, p)
        return "%s<%d" % (p, self.tail)
    __repr__ = __str__

# end class Segment_Pointer

class Graph_Segment :

    def __init__ (self, start, end, dir) :
        self.start = Segment_Pointer (start, -1)
        self.end   = Segment_Pointer (end, 1)
        self.dir   = dir
        self.pred  = None
        self.succ  = None
    # end def __init__

    def split (self, node) :
        """ Always split *against* direction dir
            Return the edge resulting from split operation
        """
        assert self.start.tail <= node <= self.end.tail
        if node == self.start.tail and self.dir == 1 :
            if not self.pred :
                return None
            return (self.start.tail, self.pred.end.tail)
        if node == self.end.tail and self.dir == -1 :
            if not self.succ :
                return None
            return (self.end.tail, self.succ.start.tail)
        return (node, (node - self.dir))
    # end def split

    @property
    def enext (self) :
        return self.end.ptr
    # end def enext

    @property
    def eprev (self) :
        return self.start.ptr
    # end def eprev

    @property
    def t_start (self) :
        """ Transitive version of start
        """
        if self.pred :
            return self.pred.start
        return self.start
    # end def transitive_start

    @property
    def t_end (self) :
        """ Transitive version of end
        """
        if self.succ :
            return self.succ.end
        return self.end
    # end def transitive_start

    @classmethod
    def key (cls, item) :
        if isinstance (item, cls) :
            return item.start
        return item
    # end def key

    def __str__ (self) :
        sg = '+' if self.dir > 0 else '-'
        return "%s(%s)%s" % (self.start, sg, self.end)
    __repr__ = __str__

    def __eq__ (self, other) :
        if isinstance (other, self.__class__) :
            return self.start.tail == other.start.tail
        return self.start == other
    # end def __eq__

    def __lt__ (self, other) :
        if isinstance (other, self.__class__) :
            return self.start.tail < other.start.tail
        return self.start.tail < other
    # end def __lt__

# end class Graph_Segment

class Segmented_Graph :
    """ Graph initialized with a given broken edge.
        This always keeps the pieces sorted in the direction of the end
        node given in the constructor.
    """

    def __init__ (self, dim, start, end) :
        self.dim = dim
        self.t1  = start
        self.tn  = end
        self.i   = 0
        if start == 0 and end == dim - 1 or end == 0 and start == dim - 1 :
            dir = -1
            if end == 0 :
                dir = 1
                start, end = end, start
            self.segments = [Graph_Segment (start, end, dir)]
        elif start < end :
            g1 = Graph_Segment (0, start, 1)
            g2 = Graph_Segment (end, dim - 1, 1)
            g2.succ = g1
            g1.pred = g2
            self.segments = [g1, g2]
        else :
            g1 = Graph_Segment (start, dim - 1, -1)
            g2 = Graph_Segment (0, end, -1)
            g1.succ = g2
            g2.pred = g1
            self.segments = [g2, g1]
    # end def __init__

    def get_segment_index (self, node) :
        # FIXME: Would be nice to use Graph_Segment.key as key
        # parameter to bisect, only supported in python 3.10
        #idx = bisect_left (self.segments, node, key = Graph_Segment.key)
        idx = bisect_left (self.segments, node)
        if idx == len (self.segments) or self.segments [idx].start.tail > node :
            idx -= 1
        assert 0 <= idx < len (self.segments)
        return idx
    # end def get_segment_index
    
    def split_edge (self, node) :
        """ Return the edge resulting from splitting at node.
            May return None if splitting is not possible.
        """
        idx = self.get_segment_index (node)
        return self.segments [idx].split (node)
    # end def split_edge

    def split (self, node, edge_next) :
        """ Perform the split at node
            Splits the given segment at node and inserts the two new
            segments in its place. Joins node to edge.
            Invert (change direction) of all segments in the direction
            of the split (Which is the counter-direction of the original
            segment).
        """
        assert edge_next == self.tn
        idx1    = self.get_segment_index (node)
        edge    = self.split_edge (node)
        idx2    = self.get_segment_index (edge [1])
        self.tn = edge [1]
        self.i += 1
        if idx1 == idx2 :
            seg = self.segments [idx1]
            if edge [0] < edge [1] :
                dir = 1
                st  = idx1 + 2
                end = len (self.segments) + 1
                f1 = Graph_Segment (seg.start.tail, edge [0], seg.dir)
                f2 = Graph_Segment (edge [1], seg.end.tail, seg.dir * -1)
                self.segments [idx1:idx1+1] = [f1, f2]
                assert edge [0] == f1.end.tail
                f1.end.set_ptr (edge_next, self.i)
            else :
                dir = -1
                st  = idx1 - 1
                end = -1
                f1 = Graph_Segment (seg.start.tail, edge [1], seg.dir * -1)
                f2 = Graph_Segment (edge [0], seg.end.tail, seg.dir)
                self.segments [idx1:idx1+1] = [f1, f2]
                assert edge [0] == f2.start.tail
                f2.start.set_ptr (edge_next, self.i)
            if seg.succ :
                f2.succ = seg.succ
                f2.succ.pred = f2
                f2.succ.dir  = f2.dir
            if seg.pred :
                f1.pred = seg.pred
                f1.pred.succ = f1
                f1.pred.dir  = f1.dir
            if seg.eprev is not None :
                f1.start.set_ptr (seg.eprev, seg.start.idx)
            if seg.enext is not None :
                f2.end.set_ptr (seg.enext, seg.end.idx)
        else :
            if edge [0] < edge [1] :
                one, two = idx1, idx2
                end = -1
                assert self.segments [idx1].start.tail == edge [0]
                assert self.segments [idx2].end.tail == edge [1]
                self.segments [idx1].start.set_ptr (edge_next, self.i)
            else :
                two, one = idx1, idx2
                end = len (self.segments)
                assert self.segments [idx1].end.tail == edge [0]
                assert self.segments [idx2].start.tail == edge [1]
                self.segments [idx1].end.set_ptr (edge_next, self.i)
            assert self.segments [two].succ == self.segments [one]
            assert self.segments [one].pred == self.segments [two]
            self.segments [two].succ = None
            self.segments [one].pred = None
        eidx = self.get_segment_index (edge_next)
        es   = self.segments [eidx]
        assert edge_next == es.t_start.tail or edge_next == es.t_end.tail
        if es.t_start.tail == es.t_end.tail :
            assert es.t_start.ptr is not None or es.t_end.ptr is not None
            if es.t_start.ptr is None :
                es.t_start.set_ptr (edge [0], self.i)
            else :
                es.t_end.set_ptr (edge [0], self.i)
        elif edge_next == es.t_start.tail :
            es.t_start.set_ptr (edge [0], self.i)
        else :
            es.t_end.set_ptr (edge [0], self.i)
        self.fix_directions (edge [1])
    # end def split

    def fix_directions (self, n) :
        """ Proceed through all new edges and turn segments that are in
            the wrong direction until we hit a segment that is ok or
            we're at the end.
        """
        se = None
        f  = self.segments [self.get_segment_index (n)]
        while True :
            assert se == f.t_end.ptr or se == f.t_start.ptr
            # One of those two will not be identical
            if  (  (se == f.t_start.ptr  and f.t_start.ptr  != f.t_end.ptr)
                or (n  == f.t_start.tail and f.t_start.tail != f.t_end.tail)
                ) :
                assert n == f.t_start.tail
                f.dir = 1
                nn = f.t_end.ptr
                se = f.t_end.tail
            else :
                assert n == f.t_end.tail
                f.dir = -1
                nn = f.t_start.ptr
                se = f.t_start.tail
            if f.succ :
                f.succ.dir = f.dir
            if f.pred :
                f.pred.dir = f.dir
            if nn is None :
                assert se == self.t1
                return
            n = nn
            f = self.segments [self.get_segment_index (n)]
    # end def fix_directions

    def walk (self, last_idx) :
        idx  = self.t1
        prev = None
        c    = 0
        dir  = None
        while c < self.dim :
            f = self.segments [self.get_segment_index (idx)]
            if dir is None :
                if prev is None :
                    dir = -f.dir
                else :
                    if c == self.dim - 1 :
                        dir = 1
                    elif f.start.tail == f.end.tail :
                        assert prev == f.start.ptr or prev == f.end.ptr
                        if f.start.ptr == prev :
                            dir = 1
                        else :
                            dir = -1
                    else :
                        dir = 1 if f.start.tail == idx else -1
            e = f.end
            if dir < 0 :
                e = f.start
            for i in range (idx, e.tail + dir, dir) :
                yield i
                c += 1
            if e.idx is not None and e.idx <= last_idx :
                idx  = e.ptr
                prev = e.tail
                dir  = None
            else :
                idx = (i + dir) % self.dim
                # keep dir
    # end def walk

# end class Segmented_Graph

class TSP (pga.PGA) :

    # eil75 with rand-seed 2 and 0.8 or-op (max 4) yields 535
    # eil50 with rand-seed 1 and 0.8 or-op (max 4) yields 426
    minvals = \
        { 'oliver30.tsp'  : 420 # 421 according to WSF89
        #, 'eil50.tsp'     : 426 # 428 according to WSF89
        , 'eil51.tsp'     : 426
        , 'eil75.tsp'     : 535 # 545 according to WSF89
        , 'eil76.tsp'     : 538
        , 'eil101.tsp'    : 629
        , 'grid.tsp'      :  72
        , 'croes.tsp'     : 246
        , 'dantzig42.tsp' : 699
        , 'gr24.tsp'      : 1272
        , 'gr48.tsp'      : 5046
        , 'gr120.tsp'     : 6942
        }

    def __init__ (self, args) :
        self.args = args
        self.tsp  = tspload (args.tsplibfile)
        self.tsp_offset = 0
        try :
            self.tsp.get_weight (0, 0)
        except (IndexError, KeyError) :
            self.tsp_offset = 1
        self.minval = self.minvals.get (self.args.tsplibfile, None)
        d = dict \
            ( random_seed      = args.random_seed
            , maximize         = False
            , max_GA_iter      = self.args.max_iter
            , pop_size         = args.popsize
            , num_replace      = round (args.popsize / 10)
            #, num_replace      = round (args.popsize / 20)
            #, num_replace      = round (args.popsize / 5)
            #, num_replace      = round (args.popsize * 2 / 3)
            #, num_replace      = args.popsize
            #, select_type      = pga.PGA_SELECT_LINEAR
            , select_type      = pga.PGA_SELECT_TOURNAMENT
            #, tournament_size  = 1.7
            #, tournament_size  = 4
            , tournament_size  = 2.5
            #, tournament_size  = 2.0
            , mixing_type      = pga.PGA_MIX_TRADITIONAL
            , pop_replace_type = pga.PGA_POPREPL_BEST
            #, pop_replace_type = pga.PGA_POPREPL_RTR
            ##, rtr_window_size  = 20
            , no_duplicates    = True
            , crossover_type   = pga.PGA_CROSSOVER_EDGE
            #, crossover_prob   = 1.0
            , print_options    = [pga.PGA_REPORT_STRING]
            #, tournament_with_replacement = False
            )
        self.fixed_edges = set ()
        if self.tsp.fixed_edges :
            fe = np.array (self.tsp.fixed_edges) - 1
            d.update (fixed_edges = fe)
            for t in fe :
                self.fixed_edges.add (tuple (t))
                self.fixed_edges.add (tuple (reversed (t)))
        super (self.__class__, self).__init__ (int, self.tsp.dimension, **d)
        self.random  = PGA_Random (self)
        self.shuffle = list (range (len (self)))
        self.random.shuffle (self.shuffle)
        self.sidx    = 0
        self.v_idx1  = list (range (len (self)))
        self.n_longest = round (len (self) * self.args.long_edge_ratio)
        if self.n_longest < 4 :
            self.n_longest = 4
        self.long_edge_probability = self.args.long_edge_probability
        if not self.long_edge_probability :
            self.long_edge_probability = self.n_longest / len (self)
        # Statistics
        self.tries             = 0
        self.or_op_tries       = 0
        self.two_op_tries      = 0
        self.or_op_success     = 0
        self.two_op_success    = 0
        self.try_harder_fail   = 0
        self.hard_tries        = 0
        self.hard_success      = 0
        self.normal_fail       = 0
        self.long_edge_tries   = 0
        self.long_edge_success = 0
        # For brute-force search note edges that failed completely
        self.good_edges        = {}
        self.prune_edges       = {}
        self.good_edge_hits    = 0
        # Parameters of brute-force search
        self.prune_value       = 10
        self.n_generation      = 300
        self.n_remove          = 0
    # end def __init__

    def edge_weight (self, i, j) :
        return self.tsp.get_weight (j + self.tsp_offset, i + self.tsp_offset)
    # end def edge_weight

    def evaluate (self, p, pop) :
        s = 0
        for i in range (self.tsp.dimension) :
            j = (i + 1) % self.tsp.dimension
            a1 = self.get_allele (p, pop, i)
            a2 = self.get_allele (p, pop, j)
            s += self.edge_weight (a1, a2)
        return s
    # end def evaluate

    def plot (self, idx = None, pop = None) :
        if pop is None :
            pop = pga.PGA_OLDPOP
        if idx is None :
            idx = self.get_best_index (pop)
        X = []
        Y = []
        for i in range (self.tsp.dimension) :
            a = self.get_allele (idx, pop, i)
            if not self.tsp.node_coords and not self.tsp.display_data :
                print ("No coordinates to plot", file = sys.stderr)
                return
            try :
                node = self.tsp.node_coords [a+1]
            except KeyError :
                node = self.tsp.display_data [a+1]
            X.append (node [0])
            Y.append (node [1])
        X.append (X [0])
        Y.append (Y [0])
        plt.plot (X, Y)
        plt.title ('%s: %s' % (self.tsp.name, self.get_best_report (pop, 0)))
        plt.show ()
    # end def plot

    def stop_cond (self) :
        if self.minval is not None :
            best = self.get_best_report (pga.PGA_OLDPOP, 0)
            if best <= self.minval :
                return True
        if self.check_stopping_conditions () :
            return True
        return False
    # end stop_cond

    def print_string (self, file, p, pop) :
        print ("evals: %s" % self.eval_count, file = file)
        print ("best index: %d" % self.get_best_index (pop), file = file);
        print ( "Tries: %d 2op: %d/%d Or-op: %d/%d Fail: %d"
              % ( self.tries
                , self.two_op_success, self.two_op_tries
                , self.or_op_success, self.or_op_tries
                , self.normal_fail
                )
              )
        print ( "Long edges: %d/%d"
              % (self.long_edge_success, self.long_edge_tries)
              )
        if self.args.optimize_harder :
            print ( "Hard-tries: %d/%d hardfail: %d"
                  % (self.hard_success, self.hard_tries, self.try_harder_fail)
                  )
            print ( "Good edges: %d prune: %d hits: %d"
                  % ( len (self.good_edges)
                    , len (self.prune_edges)
                    , self.good_edge_hits
                    )
                  )
        print ("")
        l      = len (self)
        allele = [self.get_allele (p, pop, i) for i in range (l)]
        fe = 0
        # Find first fixed edge if any
        if self.fixed_edges :
            for i in range (l) :
                j = (i + 1) % l
                if (allele [i], allele [j]) in self.fixed_edges :
                    fe = j
                    break
        print (" ".join (str (allele [(i + fe) % l] + 1) for i in range (l)))
        file.flush ()
        #return super (self.__class__, self).print_string (file, p, pop)
    # end def print_string

    def two_op (self, allele, idx1, idx2) :
        """ Exchange two edges if result is better
            This operation was first conceived by Croes 1958 who calls
            it an "inversion" operation.
            Input is the allele and the two nodes which start the edge
            (in the direction set by the gene)
        """
        assert idx1 != idx2
        l  = len (self)
        i1 = allele [idx1]
        j1 = allele [(idx1 + 1) % l]
        i2 = allele [idx2]
        j2 = allele [(idx2 + 1) % l]
        eo = self.edge_weight (i1, j1) + self.edge_weight (i2, j2)
        en = self.edge_weight (i1, i2) + self.edge_weight (j1, j2)
        if en < eo :
            # Invert half of the cycle
            for i in range (l) :
                ap1  = allele [(idx1 + i + 1) % l]
                ap2  = allele [(idx2 - i) % l]
                if ap1 == ap2 :
                    break
                allele [(idx1 + i + 1) % l] = ap2
                allele [(idx2 - i) % l]     = ap1
                if (idx1 + i + 2) % l == (idx2 - i) % l :
                    break
            return eo - en
        return 0
    # end def two_op

    def or_op (self, allele, idx1, idx2, idx3, inv, force = False) :
        """ Operation by Ilhan Or, 1976
            Move n consecutive nodes to another place, optionally
            inverting the sequence on re-insertion. We specify three
            nodes after which the edge is split. We move n nodes *after*
            idx1 to the point *after* idx3. The idx2 defines the length
            of the sequence to move. So idx2 must be > idx1.
        """
        l   = len (self)
        n   = (idx2 - idx1) % l
        st  = (idx1 + 1) % l
        assert idx1 != idx3
        assert n > 0
        eo  = 0
        seq = []
        lv  = None
        for i in range (n) :
            ix  = (st + i) % l
            assert (ix != idx3)
            v = allele [ix]
            if lv is not None :
                eo += self.edge_weight (lv, v)
            lv = v
            seq.append (v)
        eb  = eo
        bef = allele [(st - 1) % l]
        aft = allele [(st + n) % l]
        eo += self.edge_weight (bef, seq [0])
        eo += self.edge_weight (seq [-1], aft)
        en  = self.edge_weight (bef, aft)
        y1  = allele [idx3]
        y2  = allele [(idx3 + 1) % l]
        eo += self.edge_weight (y1, y2)
        if inv :
            en += self.edge_weight (y1, seq [-1])
            en += self.edge_weight (seq [0], y2)
            for i in range (n - 1) :
                en += self.edge_weight (seq [i], seq [i + 1])
        else :
            en += eb
            en += self.edge_weight (y1, seq [0])
            en += self.edge_weight (seq [-1], y2)
        if (en < eo or force) :
            for i in range (l) :
                ix1 = (st + i) % l
                ix2 = (st + i + n) % l
                allele [ix1] = allele [ix2]
                if ix2 == idx3 :
                    break
            iter = seq
            if inv :
                iter = reversed (seq)
            for i, v in enumerate (iter) :
                ix = (ix1 + i + 1) % l
                allele [ix] = v
            return eo - en
        return 0
    # end def or_op

    def lk_candidates (self, allele, t2, ewo) :
        l = len (self)
        candidates = []
        for idx in range (l) :
            if idx == t2 or idx == self.t1 :
                continue
            edge = self.lk_graph.split_edge (idx)
            if edge is None :
                continue
            idx2 = edge [1]
            if idx2 == self.t1 :
                continue
            if (idx, idx2) in self.lk_joined :
                continue
            if (idx, t2) in self.lk_joined :
                continue
            if (idx, idx2) in self.lk_broken :
                continue
            if (t2, idx) in self.lk_broken :
                continue
            if (allele [idx], allele [idx2]) in self.fixed_edges :
                continue
            # This would yield the same edge as the one being split:
            if idx2 == t2 :
                continue
            # Lookahead: FIXME
            #if (idx2, (idx2 - dir) % l) in self.lk_joined :
            #    continue
            # Gain condition
            ewn = self.edge_weight (allele [t2], allele [idx])
            if self.lk_gain + ewo - ewn <= self.lk_best_g :
                continue
            candidates.append ((edge, ewo - ewn))
        candidates.sort (key = lambda c: -c [1])
        return candidates
    # end def lk_candidates

    def is_valid_tour (self, allele) :
        l = len (self)
        d = {}
        for i in range (l) :
            a = allele [i]
            if a < 0 or a >= l :
                print (allele)
                return False
            if a in d :
                print (allele)
                return False
            d [a] = 1
        return True
    # end def is_valid_tour

    def lk_next (self, allele, t1, t2) :
        self.lk_i += 1
        l = len (self)
        self.lk_broken [(t1, t2)] = 1
        self.lk_broken [(t2, t1)] = 1
        ewo = self.edge_weight (allele [t1], allele [t2])
        candidates = self.lk_candidates (allele, t2, ewo)
        for (tk1, tk2), g in candidates :
            assert (t2, tk1)  not in self.lk_joined
            assert (tk1, tk2) not in self.lk_broken
            ewo = self.edge_weight (allele [tk1], allele [tk2])
            ewc = self.edge_weight (allele [tk2], allele [self.t1])
            gn  = g - ewc + ewo
            if gn + self.lk_gain > self.lk_best_g :
                self.lk_best_g = gn + self.lk_gain
                self.lk_best_i = self.lk_i
                self.lk_best_t = tk2
            self.lk_gain += g
            if self.args.debug :
                print ( "g: %5d gn: %5d lkg: %5d best_g: %5d"
                      % (g, gn, self.lk_gain, self.lk_best_g)
                      )
            self.lk_joined [(t2, tk1)] = -1
            self.lk_joined [(tk1, t2)] = -1
            self.lk_edges.append ((t2, tk1))
            oldgraph = deepcopy (self.lk_graph)
            self.lk_graph.split (tk1, t2)
            if self.args.debug >= 2 :
                a = []
                for i in self.lk_graph.walk (self.lk_best_i) :
                    a.append (allele [i])
                print (np.array (a) + 1)
            # Debug:
            if self.args.debug >= 3 :
                a = self.lk_take_tour (allele)
                assert self.is_valid_tour (a)

            self.lk_next (allele, tk1, tk2)
            if self.lk_best_g :
                return
            del self.lk_joined [(t2, tk1)]
            del self.lk_joined [(tk1, t2)]
            del self.lk_broken [(tk1, tk2)]
            del self.lk_broken [(tk2, tk1)]
            self.lk_edges.pop ()
            self.lk_gain -= g
            self.lk_graph = oldgraph
        self.lk_i -= 1
    # end def lk_next

    def lk_take_tour (self, allele) :
        an = []
        for i in self.lk_graph.walk (self.lk_best_i) :
            an.append (allele [i])
        #print (an)
        return an
    # end def lk_take_tour

    def lk_op (self, allele, t1) :
        """ Lin-Kernighan Heuristics (1973)
            This starts with an edge from node t1 and tries to find
            multiple additional edges to replace.
        """
        self.lk_broken = {}
        self.lk_joined = {}
        self.lk_edges  = []
        self.t1        = t1
        l = len (self)
        self.lk_gain   = 0
        self.lk_best_g = 0
        self.lk_best_i = -1
        # The very first iteration tries to break before and after t1
        self.lk_i = 0
        for dir in -1, 1 :
            t2 = (t1 + dir) % l
            self.lk_graph = Segmented_Graph (l, t1, t2)
            self.lk_next (allele, t1, t2)
            if self.lk_best_g :
                n_allele = self.lk_take_tour (allele)
                if self.args.debug :
                    print (np.array (n_allele) + 1)
                return self.lk_best_g, n_allele
    # end def lk_op

    def lk_optimize (self, p = 0, pop = pga.PGA_OLDPOP) :
        """ Optimize given individuum
        """
        l = len (self)
        shuffle = [i for i in range (l)]
        allele = [self.get_allele (p, pop, i) for i in range (l)]
        #print (allele)
        #print ("Eval: %s" % self.evaluate (p, pop))
        while True :
            self.random.shuffle (shuffle)
            for idx in shuffle :
                #print ("Try: %s" % idx)
                r = self.lk_op (allele, idx)
                if r is not None :
                    gain, n_allele = r
                    #print ("gain: %s" % gain)
                    for i in range (l) :
                        self.set_allele (p, pop, i, n_allele [i])
                    allele = n_allele
                    #print (allele)
                    #print ("Eval: %s" % self.evaluate (p, pop))
                    break
            else :
                break
        print ("Eval: %s" % self.evaluate (p, pop))
        return allele
    # end def lk_optimize

    def long_edges (self, allele) :
        l = len (self)
        edges = []
        for i in range (l) :
            j = (i + 1) % l
            edge = (i, j)
            if (allele [i], allele [j]) in self.fixed_edges :
                continue
            edges.append (edge)
        edges.sort \
            (key = lambda e: -self.edge_weight (allele [e [0]], allele [e [1]]))
        return edges [:self.n_longest]
    # end def long_edges

    def valid_index (self, new_index, allele, idx, n) :
        # Yech modulo in window that may overlap
        l = len (self)
        if idx <= new_index <= idx + n :
            return False
        if (idx + n) % l < idx and 0 <= new_index <= (idx + n) % l :
            return False
        if new_index == (idx - 1) % l :
            return False
        if  ( (allele [new_index], allele [(new_index + 1) % l])
              in self.fixed_edges
            ) :
            return False
        return True
    # end def valid_index

    def next_shuffle_index (self, allele, idx, n = 0) :
        l = len (self)
        while not self.valid_index (self.shuffle [self.sidx], allele, idx, n) :
            self.sidx += 1
            if self.sidx >= len (self) :
                self.random.shuffle (self.shuffle)
                self.sidx = 0
        r = self.shuffle [self.sidx]
        self.sidx += 1
        if self.sidx >= len (self) :
            self.random.shuffle (self.shuffle)
            self.sidx = 0
        return r
    # end def next_shuffle_index

    def update_eval (self, p, pop, v) :
        ev = self.get_evaluation (p, pop)
        self.set_evaluation (p, pop, ev - v)
        if self.evaluate (p, pop) != self.get_evaluation (p, pop) :
            import pdb; pdb.set_trace ()
    # end def update_eval

    def try_or_op_two_op (self, allele) :
        long_edge_iter = Long_Edge_Iter (self, allele)
        l = len (self)
        self.random.shuffle (self.v_idx1)
        if self.args.or_op_for_gene :
            orop = self.random_flip (self.args.or_op_probability)
        eval = 0
        for idx in self.v_idx1 :
            if self.random_flip (1 - self.args.end_of_gene_probability) :
                continue
            self.tries += 1
            if not self.args.or_op_for_gene :
                orop = self.random_flip (self.args.or_op_probability)
            if orop :
                self.or_op_tries += 1
                il = (idx - 1) % l
                if (allele [idx], allele [il]) in self.fixed_edges :
                    continue
                n    = self.random_interval (1, self.args.or_op_max)
                if  ( (allele [(idx + n - 1) % l], allele [(idx + n) % l])
                      in self.fixed_edges
                    ) :
                    continue
                inv  = self.random_flip (0.5)
                idx2 = self.next_shuffle_index (allele, idx, n)
                eval = self.or_op (allele, il, (il + n) % l, idx2, inv)
                if eval :
                    self.or_op_success += 1
                    break
            else :
                self.two_op_tries += 1
                if (allele [idx], allele [(idx + 1) % l]) in self.fixed_edges :
                    continue
                idx2 = self.next_shuffle_index (allele, idx)
                eval = self.two_op (allele, idx, idx2)
                if eval :
                    self.two_op_success += 1
                    break
            if self.random_flip (self.long_edge_probability) :
                e = long_edge_iter.next ()
                for inv in False, True :
                    self.long_edge_tries += 1
                    eval = self.or_op \
                        (allele, e [0][0], e [1][0], e [2][0], inv)
                    if eval :
                        self.long_edge_success += 1
                        return eval
        return eval
    # end def try_or_op_two_op

    def dont_use_edge (self, edge) :
        if edge in self.fixed_edges :
            return True
        ge = self.good_edges.get (edge)
        if ge :
            # Too long ago ?
            if self.GA_iter - ge < self.n_generation :
                self.good_edge_hits += 1
                return True
            del self.good_edges [edge]
        pe = self.prune_edges.get (edge)
        if pe and pe > self.prune_value :
            del self.prune_edges [edge]
            self.good_edges [edge] = self.GA_iter
            self.good_edges [tuple (reversed (edge))] = self.GA_iter
            return True
        return False
    # end def dont_use_edge

    def try_or_op_harder (self, allele) :
        """ Try to apply or_op with minimum length and successive
            reduction of length.
        """
        self.random.shuffle (self.shuffle)
        l = len (self)
        for idx in self.v_idx1 :
            idxedge = (allele [idx], allele [(idx - 1) % l])
            if self.dont_use_edge (idxedge) :
                continue
            for idx2 in self.shuffle :
                idx2edge = (idx2, (idx2 + 1) % l)
                if self.dont_use_edge (idxedge) :
                    continue
                dif = (idx - idx2) % l
                if dif < 2 :
                    continue
                for inv in False, True :
                    for n in range (dif, 0, -1) :
                        if  ( ( allele [(idx + n - 1) % l]
                              , allele [(idx + n) % l]
                              ) in self.fixed_edges
                            ) :
                            continue
                        if not self.valid_index (idx2, allele, idx, n) :
                            continue
                        self.tries       += 1
                        self.or_op_tries += 1
                        self.hard_tries  += 1
                        il = (idx - 1) % l
                        eval = self.or_op (allele, il, (il + n) % l, idx2, inv)
                        if eval :
                            self.or_op_success += 1
                            self.hard_success  += 1
                            return eval
                if idxedge not in self.prune_edges :
                    self.prune_edges [idxedge] = 0
                    self.prune_edges [tuple (reversed (idxedge))] = 0
                self.prune_edges [idxedge] += 1
                self.prune_edges [tuple (reversed (idxedge))] += 1
            # breaking this edge did not yield improvement so it must be good
            self.good_edges [idxedge] = self.GA_iter
            self.good_edges [tuple (reversed (idxedge))] = self.GA_iter
        self.try_harder_fail += 1
        # Remove n_remove oldest edges from good edges
        # If n_remove is 0, remove *all* and double the prune_value
        if not self.n_remove :
            self.good_edges  = {}
            self.prune_edges = {}
            self.prune_value *= 2
        else :
            for n, ge in enumerate (sorted
                ( list (self.good_edges)
                , key = lambda e : -self.good_edges [e]
                )) :
                del self.good_edges [ge]
                if n >= (self.n_remove - 1) :
                    break
        return 0
    # end def try_or_op_harder

    def endofgen (self) :
        # Generate permutation
        l = len (self)
        pop    = pga.PGA_NEWPOP
        lopt   = dict \
            ( worst = self.get_worst_index (pop)
            , best  = self.get_best_index  (pop)
            , rand  = self.random_interval (0, self.pop_size - 1)
            )
        local_opt = lopt.get (self.args.optimize_harder, -1)
        for p in range (self.pop_size) :
            assert self.get_evaluation_up_to_date (p, pop)
            allele_a = [self.get_allele (p, pop, i) for i in range (l)]
            # allele_b is later updated in-place
            allele_b = copy (allele_a)
            if p == local_opt :
                eval = self.try_or_op_harder (allele_b)
            else :
                eval = self.try_or_op_two_op (allele_b)
            if eval :
                for i, v in enumerate (allele_b) :
                    self.set_allele (p, pop, i, v)
                self.update_eval (p, pop, eval)
                if self.fixed_edges :
                    for i in range (l) :
                        j = (i + 1) % l
                        if (allele_b [i], allele_b [j]) in self.fixed_edges :
                            break
                    else :
                        print (allele_a)
                        print ('')
                        print (allele_b)
                        import pdb; pdb.set_trace ()
            else :
                self.normal_fail += 1
    # end def endofgen

    def gene_difference (self, p1, pop1, p2, pop2) :
        """ Used when RTR population replacement is used
            We count the number of common edges.
        """
        e1 = set ()
        e2 = set ()
        l  = len (self)
        o1 = self.get_allele (p1, pop1, 0)
        o2 = self.get_allele (p2, pop2, 0)
        for i in range (l) :
            j = (i + 1) % l
            a1 = self.get_allele (p1, pop1, j)
            a2 = self.get_allele (p2, pop2, j)
            if o1 > a1 :
                e1.add ((a1, o1))
            else :
                e1.add ((o1, a1))
            if o2 > a2 :
                e2.add ((a2, o2))
            else :
                e2.add ((o2, a2))
            assert o1 != a1
            assert o2 != a2
            o1 = a1
            o2 = a2
        return len (e1 & e2)
    # end def gene_difference

# end class TSP

if __name__ == '__main__' :
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( 'tsplibfile'
        , help = 'TSP-Lib compatible file to read as problem'
        )
    cmd.add_argument \
        ( '--optimize-harder'
        , help    = 'Use local optimizer for single individuum, possible '
                    'values: worst, best, rand'
        )
    cmd.add_argument \
        ( '-D', '--debug'
        , help    = 'Enable debug output and checks'
        , action  = 'count'
        , default = 0
        )
    cmd.add_argument \
        ( '-e', '--end-of-gene-probability'
        , help    = 'End-of-gene optimization probability, default=%(default)s'
        , type    = float
        , default = 0.2
        )
    cmd.add_argument \
        ( '-G', '--or-op-for-gene'
        , help    = 'Decide Or-op vs. Two-op once per gene'
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-l', '--long-edge-ratio'
        , help    = 'Ratio of long edges to consider, default=%(default)s'
        , type    = float
        , default = 0.05
        )
    cmd.add_argument \
        ( '-L', '--long-edge-probability'
        , help    = 'Probability of long-edge heuristics, default is '
                    'long_edges/popsize'
        , type    = float
        )
    cmd.add_argument \
        ( '--lin-kernighan'
        , help    = 'Use Lin-Kernighan heuristics'
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-m', '--max-iter'
        , help    = 'Maximum generations, default=%(default)s'
        , type    = int
        , default = 1000
        )
    cmd.add_argument \
        ( '-o', '--or-op-max'
        , help    = 'Maximum length or Or-Op, default=%(default)s'
        , type    = int
        , default = 4
        )
    cmd.add_argument \
        ( '-O', '--or-op-probability'
        , help    = 'Probability Or-op vs Two-op, default=%(default)s'
        , type    = float
        , default = 0.2
        )
    cmd.add_argument \
        ( '--plot'
        , help    = 'Plot result'
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-p', '--popsize'
        , help    = 'Population size, default=%(default)s'
        , type    = int
        , default = 50
        )
    cmd.add_argument \
        ( '-R', '--random-seed'
        , help    = 'Random seed to for initializing random number generator'
        , type    = int
        , default = 42
        )
    args    = cmd.parse_args ()
    tsp     = TSP (args)
    if args.lin_kernighan :
        allele = np.array (tsp.lk_optimize ()) + 1
        print (allele)
    else :
        tsp.run ()
    if args.plot :
        tsp.plot ()

