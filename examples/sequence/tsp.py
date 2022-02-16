#!/usr/bin/python3

import pga
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from tsplib95 import load as tspload
from random   import Random

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

class TSP (pga.PGA) :

    # eil75 with rand-seed 5 and two-op yields 543
    # eil75 with rand-seed 42 and 0.8 or-op (max 4) yields 538
    # eil50 with rand-seed 3 and 0.8 or-op (max 4) yields 426
    minvals = \
        { 'oliver30.tsp'  : 420 # 421 according to WSF89
        , 'eil50.tsp'     : 426 # 428 according to WSF89
        , 'eil51.tsp'     : 426
        , 'eil75.tsp'     : 538 # 545 according to WSF89
        , 'eil76.tsp'     : 538
        , 'eil101.tsp'    : 629
        , 'grid.tsp'      :  72
        , 'croes.tsp'     : 246
        , 'dantzig42.tsp' : 699
        }

    def __init__ (self, tsp, args) :
        self.tsp        = tsp
        self.tsp_offset = 0
        try :
            self.tsp.get_weight (0, 0)
        except (IndexError, KeyError) :
            self.tsp_offset = 1
        self.args   = args
        self.minval = self.minvals.get (self.args.tsplibfile, None)
        popsize = 50
        super (self.__class__, self).__init__ \
            ( int, tsp.dimension
            , random_seed      = args.random_seed
            , maximize         = False
            , max_GA_iter      = self.args.max_iter
            , pop_size         = popsize
            , num_replace      = round (popsize / 10)
            #, num_replace      = round (popsize / 20)
            #, num_replace      = round (popsize / 5)
            #, num_replace      = round (popsize * 2 / 3)
            #, num_replace      = popsize
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
        self.random  = PGA_Random (self)
        self.shuffle = None
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
        file.flush ()
        return super (self.__class__, self).print_string (file, p, pop)
    # end def print_string

    def two_op (self, p, pop, idx1, idx2) :
        """ Exchange two edges if result is better
            This operation was first conceived by Croes 1958 who calls
            it an "inversion" operation.
            Input are the two nodes which start the edge (in the
            direction set by the gene)
        """
        assert idx1 != idx2
        l  = len (self)
        i1 = self.get_allele (p, pop, idx1)
        j1 = self.get_allele (p, pop, (idx1 + 1) % l)
        i2 = self.get_allele (p, pop, idx2)
        j2 = self.get_allele (p, pop, (idx2 + 1) % l)
        eo = self.edge_weight (i1, j1) + self.edge_weight (i2, j2)
        en = self.edge_weight (i1, i2) + self.edge_weight (j1, j2)
        if en < eo :
            # Invert half of the cycle
            for i in range (l) :
                ap1  = self.get_allele (p, pop, (idx1 + i + 1) % l)
                ap2  = self.get_allele (p, pop, (idx2 - i) % l)
                if ap1 == ap2 :
                    break
                self.set_allele (p, pop, (idx1 + i + 1) % l, ap2)
                self.set_allele (p, pop, (idx2 - i) % l, ap1)
                if (idx1 + i + 2) % l == (idx2 - i) % l :
                    break
            ev = self.get_evaluation (p, pop)
            self.set_evaluation (p, pop, ev - eo + en)
            if self.evaluate (p, pop) != ev - eo + en :
                import pdb; pdb.set_trace ()
            return True
        return False
    # end def two_op

    def or_op (self, p, pop, idx1, idx2, n, inv) :
        """ Operation by Ilhan Or, 1976
            Move n consecutive nodes to another place, optionally
            inverting the sequence on re-insertion
        """
        l   = len (self)
        assert (idx1 - 1) % l != idx2
        eo  = 0
        seq = []
        la  = None
        for i in range (n) :
            ix  = (idx1 + i) % l
            assert (ix != idx2)
            a = self.get_allele (p, pop, ix)
            if la is not None :
                eo += self.edge_weight (la, a)
            la = a
            seq.append (a)
        eb  = eo
        bef = self.get_allele (p, pop, (idx1 - 1) % l)
        aft = self.get_allele (p, pop, (idx1 + n) % l)
        eo += self.edge_weight (bef, seq [0])
        eo += self.edge_weight (seq [-1], aft)
        en  = self.edge_weight (bef, aft)
        y1  = self.get_allele (p, pop, idx2)
        y2  = self.get_allele (p, pop, (idx2 + 1) % l)
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
        if (en < eo) :
            #for i in range (l) :
            #    print (self.get_allele (p, pop, i), end = ' ')
            #print ('')
            #import pdb; pdb.set_trace ()
            for i in range (l) :
                ix1 = (idx1 + i) % l
                ix2 = (idx1 + i + n) % l
                a = self.get_allele (p, pop, ix2)
                self.set_allele (p, pop, ix1, a)
                if ix2 == idx2 :
                    break
            iter = seq
            if inv :
                iter = reversed (seq)
            for i, a in enumerate (iter) :
                ix = (ix1 + i + 1) % l
                self.set_allele (p, pop, ix, a)
            #for i in range (l) :
            #    print (self.get_allele (p, pop, i), end = ' ')
            #print ('')
            #import pdb; pdb.set_trace ()
            ev = self.get_evaluation (p, pop)
            self.set_evaluation (p, pop, ev - eo + en)
            if self.evaluate (p, pop) != self.get_evaluation (p, pop) :
                import pdb; pdb.set_trace ()
            return True
        return False
    # end def or_op

    def next_shuffle_index (self, i, n = 0) :
        if not self.shuffle :
            self.shuffle = list (range (len (self)))
            self.random.shuffle (self.shuffle)
            self.sidx    = 0
        # Yech modulo in window that may overlap
        l = len (self)
        while (  i <= self.shuffle [self.sidx] <= i + n
              or (   (i + n) % l < i
                 and 0 <= self.shuffle [self.sidx] <= (i + n) % l
                 )
              or self.shuffle [self.sidx] == (i - 1) % l
              ) :
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

    def endofgen (self) :
        # Generate permutation
        l = len (self)
        v_idx1 = list (range (l))
        pop = pga.PGA_NEWPOP
        for p in range (self.pop_size) :
            assert self.get_evaluation_up_to_date (p, pop)
            self.random.shuffle (v_idx1)
            for idx in v_idx1 :
                orop = self.random_flip (0.2)
                if orop :
                    inv  = self.random_flip (0.5)
                    n    = self.random_interval (1, 4)
                    idx2 = self.next_shuffle_index (idx, n)
                    if self.or_op (p, pop, idx, idx2, n, inv) :
                        break
                else :
                    idx2 = self.next_shuffle_index (idx)
                    if self.two_op (p, pop, idx, idx2) :
                        break
    # end def endofgen

    def gene_difference (self, p1, pop1, p2, pop2) :
        """ Used when RTR population replacement is used
            We count the number of common edges.
        """
        import pdb; pdb.set_trace ()
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
        ( '-m', '--max-iter'
        , help    = 'Maximum generations, default=%(default)s'
        , type    = int
        , default = 1000
        )
    cmd.add_argument \
        ( '-R', '--random-seed'
        , help    = 'Random seed to for initializing random number generator'
        , type    = int
        , default = 42
        )
    cmd.add_argument \
        ( '-p', '--plot'
        , help    = 'Plot result'
        , action  = 'store_true'
        )
    args    = cmd.parse_args ()
    problem = tspload (args.tsplibfile)
    tsp     = TSP (problem, args)
    tsp.run ()
    if args.plot :
        tsp.plot ()

