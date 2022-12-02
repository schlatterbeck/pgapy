#!/usr/bin/python3

from __future__ import print_function
from argparse   import ArgumentParser
import pga
import sys

class Hello_World (pga.PGA):

    alphabet  = "abcdefghijklmnopqrstuvwxyz"
    alphabet += alphabet.upper ()
    alphabet += "0123456789!,. "
    alphabet  = alphabet.encode ('ascii')

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
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        if self.args.use_custom_init:
            self.initstring = self.initstring_
        if self.args.use_custom_mutation:
            self.mutation = self.mutation_
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

    def initstring_ (self, p, pop):
        for i in range (len (self.target)):
            r = self.random_interval (0, len (self.alphabet) - 1)
            self.set_allele (p, pop, i, self.alphabet [r:r+1])
    # end def initstring_

    def mutation_ (self, p, pop, r):
        mutations = 0
        for i in range (len (self)):
            if self.random_flip (r):
                mutations += 1
                idx = self.random_interval (0, len (self.alphabet) - 1)
                self.set_allele (p, pop, i, self.alphabet [idx:idx+1])
        return mutations
    # end def mutation_

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

# end class Hello_World

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-i", "--use-custom-init"
        , help    = "Use custom init function not only mixed upper/lowercase"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "-m", "--use-custom-mutation"
        , help    = "Use custom mutation function"
        , action  = 'store_true'
        )
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
    cmd.add_argument \
        ( "-t", "--target"
        , help    = "Target string to search for, default=%(default)s"
        , default = "Hello World!"
        )
    args = cmd.parse_args (argv)
    pg = Hello_World (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
