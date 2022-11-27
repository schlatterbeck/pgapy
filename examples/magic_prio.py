#!/usr/bin/python3

from argparse import ArgumentParser
import pga
import sys

class Magic_Square (pga.PGA):
    """ This example uses priorities (aka "Random Keys", see James C.
        Bean. Genetic algorithms and random keys for sequencing and
        optimization. ORSA Journal on Computing, 6(2):154-160, 1994.)
        and sorts the numbers in the square according to evolved prios.
        There is a twist, we fill in the numbers in concentric circles
        around the middle, the hope was that this would improve things,
        compared with a naive insertion from top left to bottom right,
        but it didn't.
    """

    def __init__ (self, args):
        self.args = args
        self.n    = args.length
        nsq       = self.n**2
        self.best = self.n * (nsq + 1) / 2
        p = dict \
            ( maximize              = False
            , init                  = [[0, nsq * 100]] * nsq
            , pop_size              = args.population_size
            , num_replace           = int (args.population_size * 0.9)
            , max_GA_iter           = 500
            , max_no_change         = 400
            , print_options         = [pga.PGA_REPORT_STRING]
            , random_seed           = args.random_seed
            , crossover_type        = pga.PGA_CROSSOVER_ONEPT
            , pop_replace_type      = pga.PGA_POPREPL_RTR
            , mutation_type         = pga.PGA_MUTATION_RANGE
            , stopping_rule_types   =
                [pga.PGA_STOP_NOCHANGE, pga.PGA_STOP_MAXITER]
            )
        if args.use_euclidian_gene_distance:
            self.gene_distance = self.euclidian_distance
        if args.mutation_rate:
            p ['mutation_prob'] = args.mutation_rate
        if self.args.output_file:
            p ['output_file'] = args.output_file
        super (self.__class__, self).__init__ (float, nsq, **p)
    # end def __init__

    def circle_iter (self):
        """ Yield coordinates of concentric rings around middle.
            if self.n is odd, we first yield the middle point.
        >>> class args (object):
        ...     pass
        >>> args.random_seed = 23
        >>> args.length = 3
        >>> m = Magic_Square (args)
        >>> for p in m.circle_iter ():
        ...    print (p)
        [1, 1]
        [0, 0]
        [1, 0]
        [2, 0]
        [2, 1]
        [2, 2]
        [1, 2]
        [0, 2]
        [0, 1]

        >>> args.length = 4
        >>> m = Magic_Square (args)
        >>> for p in m.circle_iter ():
        ...    print (p)
        [1, 1]
        [2, 1]
        [2, 2]
        [1, 2]
        [0, 0]
        [1, 0]
        [2, 0]
        [3, 0]
        [3, 1]
        [3, 2]
        [3, 3]
        [2, 3]
        [1, 3]
        [0, 3]
        [0, 2]
        [0, 1]
        
        >>> args.length = 5
        >>> m = Magic_Square (args)
        >>> for p in m.circle_iter ():
        ...    print (p)
        [2, 2]
        [1, 1]
        [2, 1]
        [3, 1]
        [3, 2]
        [3, 3]
        [2, 3]
        [1, 3]
        [1, 2]
        [0, 0]
        [1, 0]
        [2, 0]
        [3, 0]
        [4, 0]
        [4, 1]
        [4, 2]
        [4, 3]
        [4, 4]
        [3, 4]
        [2, 4]
        [1, 4]
        [0, 4]
        [0, 3]
        [0, 2]
        [0, 1]
        """
        m  = (self.n - 1) / 2.
        mi = int (m)
        for k in range (mi, -1, -1):
            # concentric ring around middle
            p = [k, k]
            yield (p)
            while abs (p [0] + 1 - m) <= m - k:
                p [0] += 1
                yield p
            while abs (p [1] + 1 - m) <= m - k:
                p [1] += 1
                yield p
            while abs (p [0] - 1 - m) <= m - k:
                p [0] -= 1
                yield p
            while p [1] - 1 > k:
                p [1] -= 1
                yield p
    # end def circle_iter

    def pheno (self, p, pop):
        rsum  = [0] * self.n
        csum  = [0] * self.n
        d1sum = 0
        d2sum = 0
        rows  = [[0] * self.n for x in range (self.n)]
        g     = list \
            ( sorted
                ( range (1, len (self) + 1)
                , key = lambda x: (self.get_allele (p, pop, x - 1), x)
                )
            )
        for n, (x, y) in enumerate (self.circle_iter ()):
            v = g [n]
            rows [y][x] = v
            rsum [x] += v
            csum [y] += v
            if x == y:
                d1sum += v
            if y == self.n - x - 1:
                d2sum += v
        return rows, rsum, csum, d1sum, d2sum
    # end def pheno

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
        print ("%d iterations" % self.GA_iter, file = file)
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
    args = cmd.parse_args (argv)
    pg = Magic_Square (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
