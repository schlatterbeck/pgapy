#!/usr/bin/python3

from __future__ import print_function
from argparse import ArgumentParser
from copy import copy
import numpy as np
import pga
import sys

class Tile:

    def __init__ (self, parent, row, col):
        b = parent.magic
        self.parent = parent
        self.col    = col
        self.row    = row
        self.value  = parent.square [row, col]
        self.rerr   = parent.rsum [row] - b
        self.cerr   = parent.csum [col] - b
        self.is_d1  = col == row
        self.is_d2  = col == parent.n - row - 1
        self.d1err  = parent.d1sum if self.is_d1 else 0
        self.d2err  = parent.d2sum if self.is_d2 else 0
        self.err    = abs (self.rerr) + abs (self.cerr) \
                    + abs (self.d1err) + abs (self.d2err)
    # end def __init__

    def new_error (self, value):
        diff  = value - self.value
        rerr  = self.rerr + diff
        cerr  = self.cerr + diff
        d1err = self.d1err + diff if self.is_d1 else 0
        d2err = self.d2err + diff if self.is_d2 else 0
        return abs (rerr) + abs (cerr) + abs (d1err) + abs (d2err)
    # end def new_error

# end class Tile

class Magic_Square (pga.PGA):

    def __init__ (self, args):
        self.args   = args
        self.n      = args.length
        self.ev     = None
        nsq         = self.n**2
        self.magic  = self.n * (nsq + 1) // 2
        self.nr     = np.arange (self.n, dtype = np.int32)
        self.shape  = (self.n, self.n)
        self.rnr    = np.flip (self.nr)
        self.r      = range (nsq)
        self.square = np.zeros ((self.n, self.n), dtype = np.int32)
        if args.use_euclidian_gene_distance:
            self.gene_distance = self.euclidian_distance
        # The default mutation is exchange of two alleles, thats the
        # best for this problem from the builtin operations.
        # Crossover operator needs to be one of the permutation
        # preserving ops.
        # all numbers from 0 to n**2 - 1
        crossover_type = getattr \
            (pga, 'PGA_CROSSOVER_' + args.crossover_type.upper ())
        mutation_type = getattr \
            (pga, 'PGA_MUTATION_' + args.mutation_type.upper ())
        replace_type = getattr \
            (pga, 'PGA_POPREPL_' + args.population_replacement.upper ())
        p = dict \
            ( maximize              = False
            , random_seed           = args.random_seed
            , pop_size              = args.population_size
            , num_replace           = int (args.population_size * 0.9)
            , max_GA_iter           = 1000
            , max_no_change         = 400
            , print_options         = [pga.PGA_REPORT_STRING]
            , print_frequency       = args.print_frequency
            , select_type           = pga.PGA_SELECT_TRUNCATION
            , pop_replace_type      = replace_type
            , mutation_type         = mutation_type
            , mutation_scramble_max = 5
            , crossover_type        = crossover_type
            , no_duplicates         = args.no_duplicates
            , truncation_proportion = 0.5
            , stopping_rule_types   =
                [pga.PGA_STOP_NOCHANGE, pga.PGA_STOP_MAXITER]
            )
        if args.mutation_rate:
            p ['mutation_prob'] = args.mutation_rate
        if self.args.output_file:
            p ['output_file'] = args.output_file
        self.cache = {}
        self.cache_hits = 0
        if args.hillclimb:
            self.hillclimb = self.hillclimb_
            p ['random_deterministic'] = True
        if args.cache:
            self.pre_eval = self.pre_eval_
            self.endofgen = self.endofgen_
        super (self.__class__, self).__init__ (int, nsq, **p)
        self.random = pga.PGA_Random (self)
        self.diag_funcs = \
            [ self.try_diag1_row, self.try_diag1_col
            , self.try_diag2_row, self.try_diag2_col
            ]
        self.tile_funcs = \
            [ self.try_tile_both, self.try_tile_row, self.try_tile_col
            , self.try_tile_max, self.try_tile_best, self.try_tile_max_two
            , self.try_tile_same_multi
            ]
    # end def __init__

    @property
    def corners (self):
        n  = self.n
        return self.square \
            [np.array ([0, n-1, n-1, 0]), np.array ([0, 0, n-1, n-1])]
    # end def corners

    def endofgen_ (self):
        pop = pga.PGA_NEWPOP
        for p in range (self.pop_size):
            assert self.get_evaluation_up_to_date (p, pop)
            self.pheno (p, pop)
            self.eval_from_pheno ()
            assert self.ev == self.get_evaluation (p, pop)
            a = self.get_ind (p, pop)
            if a not in self.cache:
                self.cache [a] = self.get_evaluation (p, pop)
    # end def endofgen_

    def eval_from_pheno (self):
        self.rsum  = np.sum (self.square, axis = 1)
        self.csum  = np.sum (self.square, axis = 0)
        self.d1sum = np.sum (np.diagonal (self.square))
        self.d2sum = np.sum (self.square [self.nr, self.rnr])
        self.sr    = sorted ([(v, i) for i, v in enumerate (self.rsum)])
        self.sc    = sorted ([(v, i) for i, v in enumerate (self.csum)])
        self.ev    = \
            ( np.sum (np.abs (self.rsum - self.magic))
            + np.sum (np.abs (self.csum - self.magic))
            + abs (self.d1sum - self.magic)
            + abs (self.d2sum - self.magic)
            )
        b = self.magic
        self.errs = (self.rsum - b) [:, None] + (self.csum - b) [None, :]
        self.errs [np.diag_indices (self.n)] += self.d1sum - b
        self.errs [self.nr, self.rnr]        += self.d2sum - b
        self.max_err = np.unravel_index (self.errs.argmax (), self.shape)
        self.min_err = np.unravel_index (self.errs.argmin (), self.shape)
        self.abs_err = np.unravel_index \
            (np.abs (self.errs).argmax (), self.shape)
        return self.ev
    # end def eval_from_pheno

    def evaluate (self, p, pop):
        self.pheno (p, pop)
        return self.eval_from_pheno ()
    # end def evaluate

    def get_ind (self, p, pop):
        a = []
        for i in range (len (self)):
            a.append (self.get_allele (p, pop, i))
        return tuple (a)
    # end def get_ind

    def hillclimb_ (self, p, pop):
        self.pheno (p, pop)
        self.eval_from_pheno ()
        dirty = False
        if self.args.normalize:
            dirty = self.normalize ()
        assert self.ev is not None
        ev = self.ev
        self.eval_from_pheno ()
        assert ev == self.ev
        sr = self.sr
        sc = self.sc
        # Check if only diagonals remain
        found = False
        if sr [0][0] == sr [-1][0] and sc [0][0] == sc [-1][0]:
            s = self.random.sample (self.diag_funcs, len (self.diag_funcs))
            for f in s:
                if f ():
                    dirty = True
                    break
        else:
            s = self.random.sample (self.tile_funcs, len (self.tile_funcs))
            for f in s:
                if f ():
                    dirty = True
                    break
        if dirty:
            self.to_gene (p, pop)
        else:
            self.set_evaluation (p, pop, self.ev)
    # end def hillclimb_

    def normalize (self):
        """ A normalized magic square has the lowest number of all
            corners in the upper left corner *and* the second item in
            the first row is smaller than the second item in the first
            *column*.
        >>> a = np.array ([[1,2,3],[4,5,6],[7,8,9]])
        >>> class X:
        ...     corners   = Magic_Square.corners
        ...     normalize = Magic_Square.normalize
        ...     n = 3
        ...     square = a
        >>> x = X ()
        >>> assert (x.corners == np.array ([1, 7, 9, 3])).all ()
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        >>> x = X ()
        >>> x.square = np.rot90 (a, 1)
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        >>> x = X ()
        >>> x.square = np.rot90 (a, 2)
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        >>> x = X ()
        >>> x.square = np.rot90 (a, 3)
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        >>> x = X ()
        >>> x.square = a.T
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        >>> x = X ()
        >>> x.square = np.rot90 (a.T, 1)
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        >>> x = X ()
        >>> x.square = np.rot90 (a.T, 2)
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        >>> x = X ()
        >>> x.square = np.rot90 (a.T, 3)
        >>> x.normalize ()
        >>> assert (x.square == a).all ()
        """
        n = self.n
        corner = np.argmin (self.corners)
        dirty = False
        if corner == 0:
            if self.square [0][1] > self.square [1][0]:
                self.square = self.square.T
                dirty = True
        elif corner == 1:
            if self.square [n-2][0] > self.square [n-1][1]:
                self.square = np.flipud (self.square)
            else:
                self.square = np.rot90 (self.square, -1)
            dirty = True
        elif corner == 2:
            if self.square [n-1][n-2] > self.square [n-2][n-1]:
                self.square = np.rot90 (self.square.T, 2)
            else:
                self.square = np.rot90 (self.square, 2)
            dirty = True
        else:
            if self.square [1][n-1] > self.square [0][n-2]:
                self.square = np.fliplr (self.square)
            else:
                self.square = np.rot90 (self.square, 1)
            dirty = True
        if dirty:
            self.eval_from_pheno ()
        return dirty
    # end def normalize

    def pheno (self, p, pop):
        for i in self.nr:
            for j in self.nr:
                a = self.get_allele (p, pop, i * self.n + j) + 1
                self.square [i, j] = a
        self.dirty = False
    # end def pheno

    def pre_eval_ (self, pop):
        for p in range (self.pop_size):
            self.pheno (p, pop)
            self.eval_from_pheno ()
            if self.get_evaluation_up_to_date (p, pop):
                assert self.ev == self.get_evaluation (p, pop)
                continue
            a = self.get_ind (p, pop)
            if a in self.cache:
                assert self.ev == self.cache [a]
                self.set_evaluation (p, pop, self.cache [a])
                assert self.get_evaluation_up_to_date (p, pop)
                self.cache_hits += 1
    # end def pre_eval_

    def print_string (self, file, p, pop):
        sl = len ('%s' % int (self.magic))
        f  = '%%%dd' % sl
        fr = (f + '  ') * self.n
        self.pheno (p, pop)
        self.eval_from_pheno ()
        for r, row in enumerate (self.square):
            rt = tuple (row)
            print (' ' * (sl + 2), fr % rt + f % self.rsum [r], file = file)
        cs = tuple (self.csum)
        print ('', f % self.d2sum, '', fr % cs + f % self.d1sum, file = file)
        print ("Best idx:", self.get_best_index (pop), file = file)
        print ("Cache hits: %d" % self.cache_hits, file = file)
        file.flush ()
    # end def print_string

    def sorted_err_idx (self):
        """ Sort errors and their indeces, we need the errors only for
            sorting the index.
        """
        e  = self.errs
        return list (sorted (np.ndindex (e.shape), key = lambda x: e [x]))
    # end def sorted_err_idx

    def stop_cond (self):
        """ Stop when the evaluation has reached 0
        """
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.get_evaluation (best, pga.PGA_OLDPOP)
        if eval == 0:
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

    def to_gene (self, p, pop):
        for i, row in enumerate (self.square):
            for j, a in enumerate (row):
                self.set_allele (p, pop, i * self.n + j, a - 1)
        self.set_evaluation (p, pop, self.ev)
    # end def to_gene

    def try_diag1_row (self):
        b  = self.magic
        d1 = self.d1sum - b
        ev = self.ev
        if d1 != 0:
            for r in self.nr:
                row  = self.square [r]
                item = self.square [r, r]
                mask = np.isin (row, item - d1, assume_unique = 1)
                if mask.any ():
                    c2 = self.nr [mask][0]
                    self.square [:, [r, c2]] = self.square [:, [c2, r]]
                    if self.eval_from_pheno () < ev:
                        return True
                    self.square [:, [r, c2]] = self.square [:, [c2, r]]
                    self.eval_from_pheno ()
                    assert ev == self.ev
    # end def try_diag1_row

    def try_diag1_col (self):
        b  = self.magic
        d1 = self.d1sum - b
        ev = self.ev
        if d1 != 0:
            for r in self.nr:
                col  = self.square [:, r]
                item = self.square [r, r]
                mask = np.isin (col, item - d1, assume_unique = 1)
                if mask.any ():
                    r2 = self.nr [mask][0]
                    self.square [[r, r2], :] = self.square [[r2, r], :]
                    if self.eval_from_pheno () < ev:
                        return True
                    self.square [[r, r2], :] = self.square [[r2, r], :]
                    self.eval_from_pheno ()
                    assert ev == self.ev
    # end def try_diag1_col

    def try_diag2_row (self):
        b  = self.magic
        d2 = self.d2sum - b
        ev = self.ev
        if d2 != 0:
            for r in self.nr:
                c = self.n - 1 - r
                row = self.square [r]
                item = self.square [r, c]
                mask = np.isin (row, item - d2, assume_unique = 1)
                if mask.any ():
                    c2 = self.nr [mask][0]
                    self.square [:, [r, c2]] = self.square [:, [c2, r]]
                    if self.eval_from_pheno () < ev:
                        return True
                    self.square [:, [r, c2]] = self.square [:, [c2, r]]
                    self.eval_from_pheno ()
                    assert ev == self.ev
    # end def try_diag1_row

    def try_diag2_col (self):
        b  = self.magic
        d2 = self.d2sum - b
        ev = self.ev
        if d2 != 0:
            for r in self.nr:
                c = self.n - 1 - r
                col = self.square [:, c]
                item = self.square [r, c]
                mask = np.isin (col, item - d2, assume_unique = 1)
                if mask.any ():
                    r2 = self.nr [mask][0]
                    self.square [[r, r2], :] = self.square [[r2, r], :]
                    if self.eval_from_pheno () < ev:
                        return True
                    self.square [[r, r2], :] = self.square [[r2, r], :]
                    self.eval_from_pheno ()
                    assert ev == self.ev
    # end def try_diag1_col

    def try_tile_both (self):
        if self.sr [0][0] < self.magic and self.sr [-1][0] > self.magic:
            if self.sc [0][0] < self.magic and self.sc [-1][0] > self.magic:
                idx1 = [self.sr  [0][1], self.sc  [0][1]]
                idx2 = [self.sr [-1][1], self.sc [-1][1]]
                return self.try_flip (idx1, idx2)
    # end def try_tile_both

    def try_tile_row (self):
        if self.sr [0][0] < self.magic and self.sr [-1][0] > self.magic:
            mn  = np.argmin (self.square [self.sr  [0][1], :])
            mx  = np.argmax (self.square [self.sr [-1][1], :])
            idx1 = [self.sr  [0][1], mn]
            idx2 = [self.sr [-1][1], mx]
            return self.try_flip (idx1, idx2)
    # end def try_tile_row

    def try_tile_col (self):
        if self.sc [0][0] < self.magic and self.sc [-1][0] > self.magic:
            mn  = np.argmin (self.square [:, self.sc  [0][1]])
            mx  = np.argmax (self.square [:, self.sc [-1][1]])
            idx1 = [mn, self.sc  [0][1]]
            idx2 = [mx, self.sc [-1][1]]
            return self.try_flip (idx1, idx2)
    # end def try_tile_col

    def try_tile_max (self):
        return self.try_flip (self.max_err, self.min_err)
    # end def try_tile_max

    def try_tile_best (self):
        t = Tile (self, *self.abs_err)
        bestidx = None
        best = 0
        for i in self.nr:
            for j in self.nr:
                if (i, j) == self.abs_err:
                    continue
                c = Tile (self, i, j)
                errdiff = t.new_error (c.value) - t.err \
                        + c.new_error (t.value) - c.err
                if errdiff < best:
                    best = errdiff
                    bestidx = (i, j)
        if bestidx is not None:
            return self.try_flip (self.abs_err, bestidx)
    # end def try_tile_best

    def try_tile_max_two (self):
        """ Use the two max and the two min errors
            Try all the permutations of exchanging one of the maxima
            with one of the minima
            Try highest with 2nd lowest and vice-versa
        """
        idx = self.sorted_err_idx ()
        if self.try_flip (idx [:2], [idx [-1], idx [-2]]):
            return True
        # Cross over: Reverse only idx2
        return self.try_flip (idx [:2], idx [-2:])
    # end def try_tile_max_two

    def try_tile_same_multi (self):
        """ Same complementary value from both directions, more than 2
        """
        idx = self.sorted_err_idx ()
        dif = 0
        for i in range (len (self) // 2):
            if self.errs [idx [i]] != -self.errs [idx [-i - 1]]:
                break
            if dif == 0:
                dif = -self.errs [idx [i]]
                continue
            if -self.errs [idx [i]] != dif:
                break
        if i <= 2 or dif == 0:
            return
        # Shuffle both ends
        idx [:i]  = self.random.sample (idx [:i],  i)
        idx [-i:] = self.random.sample (idx [-i:], i)
        # And exchange a subset
        #r = self.random_interval (2, i)
        #return self.try_flip (idx [:r], idx [-r:])
        return self.try_flip (idx [:2], idx [-2:])
    # end def try_tile_same_multi

    def try_flip (self, idx1, idx2):
        ev  = self.ev
        assert ev is not None
        if not isinstance (idx1 [0], (tuple, list, np.ndarray)):
            idx1 = np.array ([idx1])
        else:
            idx1 = np.array (idx1)
        if not isinstance (idx2 [0], (tuple, list, np.ndarray)):
            idx2 = np.array ([idx2])
        else:
            idx2 = np.array (idx2)
        ir, ic = np.concatenate ((idx1, idx2)).T
        l = len (ic)
        assert l == len (ir) and l % 2 == 0
        irf, icf = np.concatenate ((idx2, idx1)).T
        self.square [ir, ic] = self.square [irf, icf]
        if self.eval_from_pheno () > ev:
            # Allow bad moves for a certain percentage
            if self.random_flip (0.4):
                return True
            self.square [ir, ic] = self.square [irf, icf]
            self.eval_from_pheno ()
            assert self.ev == ev
            return False
        return True
    # end def try_flip

# end class Magic_Square

def main (argv):
    xtype = \
        ('pmx', 'modified', 'order', 'cycle', 'obx', 'pbx', 'uox', 'aex', 'nox')
    mtype = ('permute', 'scramble', 'position')
    poprepl = ('rtr',  'best', 'pairwise_best', 'random_norep', 'random_rep')
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( '-c', '--crossover-type'
        , choices = xtype
        , help    = "Crossover type, default=%(default)s"
        , default = xtype [0]
        )
    cmd.add_argument \
        ( "--cache"
        , help    = "Use a cache: Prevents hillclimber seeing all new"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-H', '--hillclimb'
        , help    = "Do hillclimbing"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-l', '--length'
        , type    = int
        , default = 3
        , help    = "Side length of the magic square, default=%(default)s"
        )
    cmd.add_argument \
        ( '-m', '--mutation-rate'
        , type    = float
        , help    = "Mutation rate, default is 1/l**2"
        )
    cmd.add_argument \
        ( '--mutation-type'
        , help    = "Mutation type, default=%(default)s"
        , choices = mtype
        , default = mtype [0]
        )
    cmd.add_argument \
        ( "--no-duplicates"
        , help    = "Avoid duplicates"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "--normalize"
        , help    = "Use normalization on each newly created individual"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( '-p', '--population-size'
        , type    = int
        , default = 500
        , help    = "Population size, default=%(default)s"
        )
    cmd.add_argument \
        ( '--population-replacement', '--poprepl'
        , help    = "Population replacement, default=%(default)s"
        , choices = poprepl
        , default = poprepl [0]
        )
    cmd.add_argument \
        ( '-P', '--print-frequency'
        , type    = int
        , default = 10
        , help    = "Print frequency, default=%(default)s"
        )
    cmd.add_argument \
        ( '-r', '-R', '--random-seed'
        , type    = int
        , default = 23
        , help    = "Random seed, default=%(default)s"
        )
    cmd.add_argument \
        ( '--use-euclidian-gene-distance'
        , help    = "Use euclidian gene distance function"
        , action  = 'store_true'
        )
    args = cmd.parse_args (argv)
    pg = Magic_Square (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
