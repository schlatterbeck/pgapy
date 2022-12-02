#!/usr/bin/python3

from __future__ import print_function
from argparse import ArgumentParser
from copy import copy
import pga
import sys

class Magic_Square (pga.PGA):

    def __init__ (self, args):
        self.args = args
        self.n    = args.length
        nsq       = self.n**2
        self.best = self.n * (nsq + 1) / 2
        self.nr   = range (self.n)
        self.r    = range (nsq)
        if args.use_custom_mutation:
            self.mutation = self.mutation_
        if args.use_euclidian_gene_distance:
            self.gene_distance = self.euclidian_distance
        # Default initialisation for integer is permutation
        # But we need to explicitly specify the init range
        # Mutation will also be a permutation but we need a special
        # crossover operator to preserve the state where we have
        # all numbers from 1 to n**2
        p = dict \
            ( maximize              = False
            , integer_init_permute  = [1, nsq]
            , pop_size              = args.population_size
            , num_replace           = int (args.population_size * 0.9)
            , max_GA_iter           = 1000
            , max_no_change         = 400
            , print_options         = [pga.PGA_REPORT_STRING]
            , random_seed           = args.random_seed
            , select_type           = pga.PGA_SELECT_TRUNCATION
            , pop_replace_type      = pga.PGA_POPREPL_RTR
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

    def pheno (self, p, pop):
        rsum  = [0] * self.n
        csum  = [0] * self.n
        d1sum = 0
        d2sum = 0
        rows  = []
        for i in self.nr:
            row  = []
            for j in self.nr:
                row.append (self.get_allele (p, pop, i * self.n + j))
                rsum [i] += row [-1]
                csum [j] += row [-1]
                if i == j:
                    d1sum += row [-1]
                if i == self.n - j - 1:
                    d2sum += row [-1]
            rows.append (row)
        return rows, rsum, csum, d1sum, d2sum
    # end def pheno

    def invert (self, seq):
        """ Inversion operator, inspired by
            https://github.com/guisehn/genetic-magic-square-finder
            used for crossover without needs for repair
            Original paper is
            http://www.tandfonline.com/doi/abs/10.1080/10798587.2000.10642829
        >>> class args (object): pass
        >>> args.length        = 3
        >>> args.random_seed   = 23
        >>> args.mutation_rate = None
        >>> ga = Magic_Square (args)
        >>> ga.invert ([2, 7, 6, 9, 5, 1, 4, 3, 8])
        [5, 0, 5, 4, 3, 1, 0, 1, 0]
        >>> ga.invert ([2, 9, 6, 8, 5, 1, 4, 3, 7])
        [5, 0, 5, 4, 3, 1, 2, 1, 0]
        >>> ga.invert ([4, 6, 2, 7, 3, 1, 5])
        [5, 2, 3, 0, 2, 0, 0]
        """
        inv = []
        d = dict ((k, n) for n, k in enumerate (seq))
        for k in range (1, len (seq) + 1):
            i = d [k]
            c = 0
            for j in reversed (range (i)):
                if seq [j] > k:
                    c += 1
            inv.append (c)
        return inv
    # end def invert

    def deinvert (self, inv):
        """ De-Inversion operator, inspired by
            https://github.com/guisehn/genetic-magic-square-finder
            used for crossover without needs for repair
            Original paper is
            http://www.tandfonline.com/doi/abs/10.1080/10798587.2000.10642829
        >>> class args (object): pass
        >>> args.length      = 3
        >>> args.random_seed = 23
        >>> args.mutation_rate = None
        >>> ga = Magic_Square (args)
        >>> ga.deinvert ([5, 0, 5, 4, 3, 1, 2, 1, 0])
        [2, 9, 6, 8, 5, 1, 4, 3, 7]
        >>> ga.deinvert ([5, 0, 5, 4, 3, 1, 0, 1, 0])
        [2, 7, 6, 9, 5, 1, 4, 3, 8]
        """
        pos = copy (inv)
        for k in reversed (range (len (inv))):
            for j in range (k + 1, len (inv)):
                if pos [j] >= pos [k]:
                    pos [j] += 1
        seq = [0] * len (inv)
        for k in range (len (inv)):
            seq [pos [k]] = k + 1
        return seq
    # end def deinvert

    def crossover (self, p1, p2, pop1, c1, c2, pop2):
        p1g = list (self.get_allele (p1, pop1, i) for i in self.r)
        p2g = list (self.get_allele (p2, pop1, i) for i in self.r)
        p1i = self.invert (p1g)
        p2i = self.invert (p2g)
        rnd = self.random_interval (0, len (self) - 1)
        c1i = p1i [:rnd] + p2i [rnd:]
        c2i = p2i [:rnd] + p1i [rnd:]
        c1g = self.deinvert (c1i)
        c2g = self.deinvert (c2i)
        for i in self.r:
            self.set_allele (c1, pop2, i, c1g [i])
            self.set_allele (c2, pop2, i, c2g [i])
    # end def crossover

    def endofgen (self):
        pop = pga.PGA_NEWPOP
        for p in range (self.pop_size):
            assert self.get_evaluation_up_to_date (p, pop)
            a = self.get_ind (p, pop)
            if a not in self.cache:
                self.cache [a] = self.get_evaluation (p, pop)
    # end def endofgen

    def get_ind (self, p, pop):
        a = []
        for i in range (len (self)):
            a.append (self.get_allele (p, pop, i))
        return tuple (a)
    # end def get_ind

    def mutation_ (self, p, pop, r):
        """ Custom mutation operator:
            We take the shape of the square into account and do a
            crossover between rows and columns of the square.
            Crossover point and row/column to crossover with are
            selected at random.
            Inspired by "Magic Squares Construction with Genetic
            Algorithm", Adil Al-Rammahi, in Al-Takany 18(2), 2005
            Note that the paper uses this operator as a *crossover*
            operator but we use it as mutation operator.
        """
        mutations = 0
        # Rows:
        for is_row in (True, False):
            for rc in self.nr:
                if not self.random_flip (r):
                    continue
                other = rc
                while other == rc:
                    other = self.random_interval (0, self.n - 1)
                co    = self.random_interval (0, self.n - 1)
                first = self.random_flip (0.5)
                if first:
                    rng = range (co)
                else:
                    rng = range (co, self.n)
                for idx in rng:
                    if is_row:
                        i  = self.n * rc + idx
                        io = self.n * other + idx
                    else:
                        i  = self.n * idx + rc
                        io = self.n * idx + other
                    a  = self.get_allele (p, pop, i)
                    self.set_allele (p, pop, i, self.get_allele (p, pop, io))
                    self.set_allele (p, pop, io, a)
                    mutations += 1
        return mutations
    # end def mutation

    def pre_eval (self, pop):
        for p in range (self.pop_size):
            if self.get_evaluation_up_to_date (p, pop):
                continue
            a = self.get_ind (p, pop)
            if a in self.cache:
                self.set_evaluation (p, pop, self.cache [a])
                self.set_evaluation_up_to_date (p, pop, True)
                self.cache_hits += 1
    # end def pre_eval

    def evaluate (self, p, pop):
        best = self.best
        rows, rsum, csum, d1sum, d2sum = self.pheno (p, pop)
        eval = \
            ( sum (abs (x - best) for x in rsum)
            + sum (abs (x - best) for x in csum)
            + abs (d1sum - best)
            + abs (d2sum - best)
            )
        return (eval) ** (1/50.)
    # end def evaluate

    def print_string (self, file, p, pop):
        sl = len ('%s' % self.best)
        f  = '%%%ss' % sl
        rows, rsum, csum, d1sum, d2sum = self.pheno (p, pop)
        for r, row in enumerate (rows):
            print (' ' * sl, row, rsum [r], file = file)
        fitness = self.get_fitness (p, pop)
        print (f % d2sum, csum, f % d1sum, file = file)
        print ("Best idx:", self.get_best_index (pop), file = file)
        print ("Cache hits: %d" % self.cache_hits, file = file)
        file.flush ()
    # end def print_string

    def stop_cond (self):
        """ Stop when the evaluation has reached 0
        """
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval == 0:
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Magic_Square

def main (argv):
    cmd = ArgumentParser ()
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
        ( '-r', '-R', '--random-seed'
        , type    = int
        , default = 23
        , help    = "Random seed, default=%(default)s"
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
        ( '--use-euclidian-gene-distance'
        , help    = "Use euclidian gene distance function"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '--use-custom-mutation'
        , help    = "Use custom mutation operator"
        , action  = 'store_true'
        )
    args = cmd.parse_args (argv)
    pg = Magic_Square (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
