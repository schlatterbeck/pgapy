#!/usr/bin/python3

from __future__ import print_function
from argparse   import ArgumentParser
import pga
import sys

class Namefull (pga.PGA):
    """ Inspired by pgapack examples/c/namefull
        But not using a custom crossover
        Note that the duplicate checking is actually cheating: Two
        strings are equivalent if they both match the target string or
        do not match the target string at the same positions. This uses
        information from the target string. Note that this string
        matching drastically increases the chance that during init two
        strings are the "same" (both have no matches). Therefore we're
        cheating on init give it a 20% chance that the correct char is
        chosen (instead of 1/94). The "duplicate" checking would do the
        same but also cheating during init makes this faster.
    """

    def __init__ (self, args):
        self.args   = args
        self.target = target = args.target.encode ('ascii')
        d = dict \
            ( maximize       = True
            , pop_size       = 100
            , num_replace    = 90
            , mutation_prob  = 0.1
            , print_options  = [pga.PGA_REPORT_STRING]
            , random_seed    = args.random_seed
            , char_init_type = pga.PGA_CINIT_MIXED
            , no_duplicates  = True
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (bytes, len (target), **d)
    # end def __init__

    def evaluate (self, p, pop):
        r = 0.0
        for i in range (len (self)):
            c = self.get_allele (p, pop, i)
            if c == self.target [i:i+1]:
                r += 1.0
        return r
    # end def evaluate

    def check_duplicate (self, p1, pop1, p2, pop2):
        for i in range (len (self.target)):
            a = self.get_allele (p1, pop1, i)
            b = self.get_allele (p2, pop2, i)
            c = self.target [i:i+1]
            if a == c and b != c or a != c and b == c:
                return False
        return True
    # end def check_duplicate

    def hash (self, p, pop):
        s = []
        c = 0
        for i in range (len (self.target)):
            s.append (self.get_allele (p, pop, i))
            if s [-1] != self.target [i]:
                s [-1] = bytes ([13])
            else:
                c += 1
        return hash (b''.join (s))
    # end def hash

    def initstring (self, p, pop):
        for i in range (len (self.target)):
            r = self.random_flip (0.2)
            if r:
                # Cheating
                v = self.target [i]
            else:
                # All printable characters
                v = self.random_interval (32, 126)
            self.set_allele (p, pop, i, bytes ([v]))
    # end def initstring

    def mutation (self, p, pop, r):
        mutations = 0
        for i in range (len (self)):
            if self.random_flip (r):
                mutations += 1
                v = self.random_interval (32, 126)
                self.set_allele (p, pop, i, bytes ([v]))
        return mutations
    # end def mutation

    def print_string (self, file, p, pop):
        b = []
        for i in range (len (self)):
            b.append (self.get_allele (p, pop, i))
        print (b''.join (b).decode ('utf-8'), file = file)
        super ().print_string (file, p, pop)
        print ("Evaluations: %d" % self.eval_count, file = file)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval == len (self):
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Namefull

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( "-R", "--random-seed"
        , help    = "Seed random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    minlen = 30
    cmd.add_argument \
        ( "-t", "--target"
        , help    = "Target string to search for, default=%%(default)s"
                    ", minimum length: %d" % minlen
        , default = "This is a simple string match: Hello World!"
        )
    args = cmd.parse_args (argv)
    if len (args.target) < minlen:
        print \
            ( "Target string must be at least %d characters long"
            % minlen
            , file = sys.stderr
            )
        sys.exit (1)
    pg = Namefull (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
