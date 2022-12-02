#!/usr/bin/python3

from argparse import ArgumentParser
import pga
import sys

class Sorter (pga.PGA):

    def __init__ (self, args):
        self.args = args
        d = dict \
            ( maximize      = False
            , pop_size      = 30
            , num_replace   = 28
            , mutation_prob = 0.3
            , init          = [(0, 10 * args.length - 1)] * args.length
            , random_seed   = args.random_seed
            , print_options = [pga.PGA_REPORT_STRING]
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super (self.__class__, self).__init__ (int, args.length, **d)
    # end def __init__

    def evaluate (self, p, pop):
        r = 0.0
        last = None
        violation = 1.0
        distance  = 0.0
        for i in range (1, len (self)):
            last = self.get_allele (p, pop, i - 1)
            this = self.get_allele (p, pop, i)
            if last >= this:
                violation *= 100
                distance  += (last + 1 - this)
            else:
                distance  += this - last
        return distance * violation
    # end def evaluate

    def print_string (self, file, p, pop):
        print \
            (', '.join \
                (str (self.get_allele (p, pop, x)) for x in range (len (self)))
            , file = file
            )
    # end def print_string

# end class Sorter

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( "-l", "--length"
        , help    = "Length of string to sort"
        , type    = int
        , default = 10
        )
    cmd.add_argument \
        ( "-R", "--random-seed"
        , help    = "Seed random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    args = cmd.parse_args (argv)
    pg = Sorter (args)
    pg.run ()

if __name__ == '__main__':
    main (sys.argv [1:])
