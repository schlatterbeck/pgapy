#!/usr/bin/python3

from __future__ import print_function
from argparse import ArgumentParser
from copy import copy
import numpy as np
import pga
import sys

class Magic_Square (pga.PGA):

    def __init__ (self, args):
        self.args   = args
        self.n      = args.length
        nsq         = self.n**2
        self.best   = self.n * (nsq + 1) / 2
        self.nr     = range (self.n)
        self.r      = range (nsq)
        self.square = np.zeros ((self.n, self.n), dtype = np.int32)
        self.dirty  = False # need to update gene
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
        p = dict \
            ( maximize              = False
            , random_seed           = args.random_seed
            #, random_deterministic  = True
            , pop_size              = args.population_size
            , num_replace           = int (args.population_size * 0.9)
            , max_GA_iter           = 1000
            , max_no_change         = 400
            , print_options         = [pga.PGA_REPORT_STRING]
            , print_frequency       = args.print_frequency
            , select_type           = pga.PGA_SELECT_TRUNCATION
            , pop_replace_type      = pga.PGA_POPREPL_RTR
            , mutation_type         = mutation_type
            , mutation_scramble_max = 5
            , crossover_type        = crossover_type
            , no_duplicates         = True
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
        super (self.__class__, self).__init__ (int, nsq, **p)
    # end def __init__

    @property
    def corners (self):
        n  = self.n
        return self.square \
            [np.array ([0, n-1, n-1, 0]), np.array ([0, 0, n-1, n-1])]
    # end def corners

    def endofgen (self):
        pop = pga.PGA_NEWPOP
        for p in range (self.pop_size):
            assert self.get_evaluation_up_to_date (p, pop)
            a = self.get_ind (p, pop)
            if a not in self.cache:
                self.cache [a] = self.get_evaluation (p, pop)
    # end def endofgen

    def evaluate (self, p, pop):
        self.pheno (p, pop)
        eval = \
            ( sum (abs (x - self.best) for x in self.rsum)
            + sum (abs (x - self.best) for x in self.csum)
            + abs (self.d1sum - self.best)
            + abs (self.d2sum - self.best)
            )

        return eval
    # end def evaluate

    def get_ind (self, p, pop):
        a = []
        for i in range (len (self)):
            a.append (self.get_allele (p, pop, i))
        return tuple (a)
    # end def get_ind

    def pheno (self, p, pop):
        self.rsum  = [0] * self.n
        self.csum  = [0] * self.n
        self.d1sum = 0
        self.d2sum = 0
        for i in self.nr:
            for j in self.nr:
                a = self.get_allele (p, pop, i * self.n + j) + 1
                self.square [i, j] = a
                self.rsum [i] += a
                self.csum [j] += a
                if i == j:
                    self.d1sum += a
                if i == self.n - j - 1:
                    self.d2sum += a
        self.dirty = False
    # end def pheno

    def pre_eval (self, pop):
        for p in range (self.pop_size):
            if self.get_evaluation_up_to_date (p, pop):
                continue
            a = self.get_ind (p, pop)
            if a in self.cache:
                self.set_evaluation (p, pop, self.cache [a])
                assert self.get_evaluation_up_to_date (p, pop)
                self.cache_hits += 1
    # end def pre_eval

    def print_string (self, file, p, pop):
        sl = len ('%s' % int (self.best))
        f  = '%%%dd' % sl
        fr = (f + '  ') * self.n
        self.pheno (p, pop)
        for r, row in enumerate (self.square):
            rt = tuple (row)
            print (' ' * (sl + 2), fr % rt + f % self.rsum [r], file = file)
        cs = tuple (self.csum)
        print ('', f % self.d2sum, '', fr % cs + f % self.d1sum, file = file)
        print ("Best idx:", self.get_best_index (pop), file = file)
        print ("Cache hits: %d" % self.cache_hits, file = file)
        file.flush ()
    # end def print_string

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
    # end def to_gene

# end class Magic_Square

def main (argv):
    xtype = \
        ('pmx', 'modified', 'order', 'cycle', 'obx', 'pbx', 'uox', 'aex', 'nox')
    mtype = ('permute', 'scramble', 'position')
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( '-c', '--crossover-type'
        , choices = xtype
        , help    = "Crossover type, default=%(default)s"
        , default = xtype [0]
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
