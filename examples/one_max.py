#!/usr/bin/python3

from operator import add
from argparse import ArgumentParser
import pga
import sys

class One_Max (pga.PGA):

    def __init__ (self, args):
        self.args = args
        d = dict \
            ( maximize      = True
            , random_seed   = args.random_seed
            , print_options = \
                [ pga.PGA_REPORT_STRING
                , pga.PGA_REPORT_WORST
                , pga.PGA_REPORT_AVERAGE
                ]
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super (self.__class__, self).__init__ (bool, args.length, **d)
    # end def __init__

    def evaluate (self, p, pop):
        return sum \
            (self.get_allele (p, pop, i) for i in range (len (self)))
    # end def evaluate

    def print_string (self, file, p, pop):
        r = []
        for k in range (10):
            r.append (str (int (self.get_allele (p, pop, k))))
        r.append ('..')
        for k in range (90, 100):
            r.append (str (int (self.get_allele (p, pop, k))))
        print (''.join (r), file = file)
        if self.args.verbose:
            super ().print_string (file, p, pop)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval >= len (self):
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class One_Max

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-l", "--length"
        , help    = "Length of bitfield, default=%(default)s"
        , type    = int
        , default = 100
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
        ( "-v", "--verbose"
        , help    = "Verbose output"
        , action  = 'store_true'
        )
    args = cmd.parse_args (argv)
    pg = One_Max (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])

