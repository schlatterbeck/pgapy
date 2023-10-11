#!/usr/bin/python3

from argparse   import ArgumentParser
import pga
import sys

class Royal_Road (pga.PGA):
    """ Royal Road function by John Holland
        See also README.rst.
    """

    def __init__ (self, args):
        self.args = args
        l = 2 ** args.region_size * (args.block_len + args.gap_len)

        d = dict \
            ( maximize                    = True
            , pop_size                    = args.pop_size
            , num_replace                 = int (0.9 * args.pop_size)
            , print_options               = [pga.PGA_REPORT_STRING]
            , random_seed                 = args.random_seed
            , print_frequency             = 1
            , crossover_type              = args.crossover_type
            , max_GA_iter                 = args.generations
            , pop_replace_type            = args.pop_replace
            , mixing_type                 = args.mixing_type
            , crossover_prob              = args.crossover_prob
            , no_duplicates               = args.no_duplicates
            , select_type                 = args.select_type
            , tournament_with_replacement = args.with_replacement
            , nam_window_size             = args.nam_window_size
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super (self.__class__, self).__init__ (bool, l, **d)
        mxvec = [[1] * self.args.block_len] * (2 ** self.args.region_size)
        self.maxeval = self.reward (mxvec)
    # end def __init__

    def evaluate (self, p, pop):
        vec = self.phenotype (p, pop)
        return self.reward (vec)
    # end def evaluate

    def phenotype (self, p, pop):
        vec = []
        l   = self.args.block_len + self.args.gap_len
        for i in range (2 ** self.args.region_size):
            bits = []
            for bi in range (self.args.block_len):
                bits.append (self.get_allele (p, pop, i * l + bi))
            vec.append (bits)
        return vec
    # end def phenotype

    def reward (self, vec):
        """ This only gets the relevant bits (non-introns)
        >>> args = get_args ([])
        >>> rr = Royal_Road (args)
        >>> ones = [1] * 8
        >>> rr.reward ([ones] * 16)
        12.8
        >>> maxpart = [1, 0, 1, 0, 1, 0, 1, 0]
        >>> rr.reward ([maxpart] * 16)
        1.28
        >>> minpart = [1, 1, 1, 1, 1, 1, 1, 0]
        >>> print ('%.3f' % rr.reward ([minpart] * 16))
        -0.960
        """
        r = 0.0
        blocks = [False] * (2 ** self.args.region_size)
        for i, bits in enumerate (vec):
            v    = 0
            ones = 0
            for bi, bit in enumerate (bits):
                v |= bit << bi
                if bit:
                    ones += 1
            if v == 2 ** self.args.block_len - 1:
                # Mark true for BONUS calculation
                blocks [i] = True
            else:
                # PART calculation
                if ones <= self.args.part_bits:
                    r += ones * self.args.part_value
                elif ones < self.args.block_len:
                    r -= (ones - self.args.part_bits) * self.args.part_value
                else:
                    assert 0
        # BONUS calculation
        for j in range (self.args.region_size + 1):
            u = self.args.bonus_u_star
            for idx in range (2 ** self.args.region_size // (2 ** j)):
                if blocks [idx]:
                    r += u
                    u = self.args.bonus_u
            itr    = iter (blocks)
            blocks = [a and b for a, b in zip (itr, itr)]
        return r
    # end def reward

    def print_string (self, file, p, pop):
        super ().print_string (file, p, pop)
        vec = self.phenotype (p, pop)
        for n, bits in enumerate (vec):
            e = ' '
            if (n + 1) % 8 == 0:
                e = '\n'
            print (''.join (str (b) for b in bits), file = file, end = e)
        print ("Evaluations: %d" % self.eval_count, file = file)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval == self.maxeval:
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Royal_Road

def get_args (argv = sys.argv [1:]):
    crosstypes = ['onept', 'twopt', 'uniform']
    repltypes  = ['best', 'nsga_ii', 'rtr']
    mixtypes   = ['mutate_or_cross', 'mutate_and_cross', 'traditional']
    seltypes   = ['tournament', 'truncation']
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-b", "--block-len"
        , help    = "Block length of optimized bits, default=%(default)s"
        , type    = int
        , default = 8
        )
    cmd.add_argument \
        ( "-c", "--crossover-type"
        , help    = "Crossover type, default=%%(default)s, one of %s"
                  % ', '.join (crosstypes)
        , default = 'twopt'
        )
    cmd.add_argument \
        ( "-C", "--crossover-prob"
        , help    = "Crossover probability, default=%(default)s"
        , type    = float
        , default = 0.85
        )
    cmd.add_argument \
        ( "-G", "--generations"
        , help    = "Max number of generations, default=%(default)s"
        , type    = int
        , default = 100
        )
    cmd.add_argument \
        ( "-g", "--gap-len"
        , help    = "Length of gap between blocks, default=%(default)s"
        , type    = int
        , default = 7
        )
    cmd.add_argument \
        ( "-k", "--region-size"
        , help    = "Base-2 logarithm of region size, default=%(default)s"
        , type    = int
        , default = 4
        )
    cmd.add_argument \
        ( "-m", "--part-bits"
        , help    = "Number of bits m* in PART calculation, default=%(default)s"
        , type    = int
        , default = 4
        )
    cmd.add_argument \
        ( "-M", "--mutation-prob"
        , help    = "Probability of mutation, default is reciprocal "
                    "of string length"
        , type    = float
        )
    cmd.add_argument \
        ( "-n", "--no-duplicates"
        , help    = "Avoid duplicates"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "-N", "--nam-window-size"
        , help    = "Negative assortative mating window size, "
                    "default=%(default)s"
        , type    = int
        , default = 1
        )
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( "-p", "--pop-size"
        , help    = "Population size, default=%(default)s"
        , type    = int
        , default = 512
        )
    cmd.add_argument \
        ( "-P", "--pop-replace"
        , help    = "Population replacement, default=%%(default)s, one of %s"
                  % ', '.join (repltypes)
        , default = 'best'
        )
    cmd.add_argument \
        ( "-R", "--random-seed"
        , help    = "Seed random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    cmd.add_argument \
        ( "-s", "--select-type"
        , help    = "Selection type, default=%%(default)s, one of %s"
                  % ', '.join (seltypes)
        , default = 'tournament'
        )
    cmd.add_argument \
        ( "-u", "--bonus-u"
        , help    = "Reward u of next block at level j, default=%(default)s"
        , type    = float
        , default = 0.3
        )
    cmd.add_argument \
        ( "-U", "--bonus-u-star"
        , help    = "Reward u* of first block at level j, default=%(default)s"
        , type    = float
        , default = 1.0
        )
    cmd.add_argument \
        ( "-v", "--part-value"
        , help    = "Value in PART computation, default=%(default)s"
        , type    = float
        , default = 0.02
        )
    cmd.add_argument \
        ( "-w", "--without-replacement"
        , dest    = 'with_replacement'
        , help    = "For tournament selection use variant without replacement"
        , action  = 'store_false'
        , default = True
        )
    cmd.add_argument \
        ( "-x", "--mixing-type"
        , help    = "Mixing type, default=%%(default)s, one of %s"
                  % ', '.join (mixtypes)
        , default = 'mutate_or_cross'
        )
    args = cmd.parse_args (argv)
    try:
        args.crossover_type = getattr \
            (pga, 'PGA_CROSSOVER_%s' % args.crossover_type.upper ())
    except AttributeError as err:
        cmd.usage ()
        sys.exit (err)
    try:
        args.pop_replace = getattr \
            (pga, 'PGA_POPREPL_%s' % args.pop_replace.upper ())
    except AttributeError as err:
        cmd.usage ()
        sys.exit (err)
    try:
        args.mixing_type = getattr \
            (pga, 'PGA_MIX_%s' % args.mixing_type.upper ())
    except AttributeError as err:
        cmd.usage ()
        sys.exit (err)
    try:
        args.select_type = getattr \
            (pga, 'PGA_SELECT_%s' % args.select_type.upper ())
    except AttributeError as err:
        cmd.usage ()
        sys.exit (err)
    return args
# end def get_args

def main (argv = sys.argv [1:]):
    pg = Royal_Road (get_args (argv))
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
