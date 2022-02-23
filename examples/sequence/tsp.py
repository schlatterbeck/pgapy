#!/usr/bin/python3

import pga
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from tsplib95 import load as tspload
from random   import Random
from copy     import copy

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

    def plot (self, pop = None) :
        if pop is None :
            pop = pga.PGA_OLDPOP
        X = []
        Y = []
        idx = self.get_best_index (pop)
        for i in range (self.tsp.dimension) :
            a = self.get_allele (idx, pop, i)
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
    cmd.add_argument \
        ( '--test'
        , help    = 'Perform some regression tests'
        , action  = 'store_true'
        )
    args    = cmd.parse_args ()
    tsp     = TSP (args)
    if args.test :
        allele = [
            0, 213, 239, 221, 227, 217, 215, 109, 216, 211, 210, 108,
            220, 230, 238, 224, 219, 107, 105, 116, 113, 112, 124, 120,
            40, 122, 121, 111, 13, 126, 134, 132, 137, 135, 136, 125,
            12, 119, 118, 117, 218, 214, 222, 235, 234, 225, 223, 226,
            212, 233, 236, 237, 232, 249, 259, 256, 264, 268, 165, 288,
            286, 276, 280, 292, 191, 295, 198, 300, 311, 204, 310, 302,
            290, 284, 279, 272, 273, 285, 289, 283, 269, 270, 303, 296,
            297, 309, 304, 274, 317, 275, 298, 299, 278, 308, 307, 301,
            305, 306, 192, 185, 180, 177, 178, 173, 179, 174, 163, 168,
            176, 183, 205, 206, 89, 69, 74, 73, 65, 64, 172, 175, 87,
            86, 186, 60, 45, 129, 57, 55, 54, 44, 51, 42, 46, 103, 36,
            152, 130, 140, 146, 164, 291, 169, 166, 167, 162, 158, 161,
            148, 147, 145, 151, 156, 143, 287, 277, 281, 282, 293, 294,
            202, 199, 197, 203, 200, 99, 196, 201, 93, 80, 195, 171,
            184, 154, 149, 144, 159, 160, 209, 181, 188, 182, 187, 189,
            190, 193, 194, 170, 314, 265, 263, 261, 255, 252, 267, 266,
            271, 251, 250, 242, 228, 229, 138, 139, 258, 254, 257, 260,
            262, 253, 246, 245, 244, 247, 243, 248, 313, 142, 241, 240,
            312, 128, 131, 127, 123, 141, 231, 155, 157, 150, 208, 153,
            33, 133, 207, 34, 38, 37, 59, 41, 62, 53, 58, 63, 61, 52,
            50, 56, 104, 67, 78, 72, 68, 88, 90, 83, 79, 76, 81, 70, 85,
            77, 82, 84, 94, 97, 98, 92, 100, 95, 96, 101, 91, 66, 71,
            75, 48, 316, 315, 49, 43, 47, 39, 35, 24, 25, 23, 26, 15,
            115, 114, 110, 106, 3, 4, 7, 2, 10, 9, 6, 5, 102, 11, 32,
            31, 29, 30, 27, 19, 16, 17, 18, 21, 20, 28, 14, 22, 8, 1]
        eval = tsp.or_op (allele, 0, 20, 125, 0, force = True)
        print (eval)
        print (allele)

    else :
        tsp.run ()
    if args.plot :
        tsp.plot ()

