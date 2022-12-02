#!/usr/bin/python3

from argparse import ArgumentParser
import pga
import sys

class Minfloat (pga.PGA):
    """ This demonstrates the use of the new Differential Evolution
        strategy. Note that there is no selection pressure in the linear
        selection method. Instead selection pressure is employed with
        the population replacement method PGA_POPREPL_PAIRWISE_BEST.
        If a pop_replace_type without selection pressure *and* a
        selection method without selection pressure is used, the
        search degenerates to a random walk.
    """

    def __init__ (self, args):
        self.args = args
        self.gene_distance = self.euclidian_distance
        popsize = 30
        d = dict \
            ( maximize                    = False
            , pop_size                    = popsize
            , num_replace                 = popsize
            , print_options               = [pga.PGA_REPORT_STRING]
            , init                        = [(-10.0, 10.0)] * 9
            , select_type                 = pga.PGA_SELECT_LINEAR
            , pop_replace_type            = pga.PGA_POPREPL_PAIRWISE_BEST
            # Either use mutation_bounded, which sets a value that
            # violates the bounds to the bound or mutation_bounce_back
            # which sets a value that violates the bound to a random
            # value between the old value and the bound. The latter
            # method is said to retain better population variance.
            #, mutation_bounded            = True
            , mutation_bounce_back        = True
            , mutation_only               = True
            , mutation_type               = pga.PGA_MUTATION_DE
            , tournament_with_replacement = False
            , random_seed                 = self.args.random_seed
            , DE_variant                  = pga.PGA_DE_VARIANT_BEST
            , DE_crossover_prob           = 0.2
            , DE_jitter                   = 0.001
            , DE_scale_factor             = 0.85 - (popsize * 0.0005)
            )
        if self.args.output_file:
            d ['output_file'] = args.output_file
        super ().__init__ (float, 9, **d)
    # end def __init__

    def evaluate (self, p, pop):
        assert not self.get_evaluation_up_to_date (p, pop)
        return sum (self.get_allele (p, pop, k) for k in range (len (self)))
    # end def evaluate

    def stop_cond (self):
        """ Stop when the evaluation has reached 0.99 * 10 * len (self)
        """
        best = self.get_best_report_index (pga.PGA_OLDPOP, 0)
        eval = self.get_evaluation (best, pga.PGA_OLDPOP)
        if eval <= -0.99999 * 10 * len (self):
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

    def print_string (self, file, p, pop):
        print ("evals: %s" % self.eval_count, file = file)
        print ("best index: %d" % self.get_best_index (pop), file = file);
        file.flush ()
        super ().print_string (file, p, pop)
    # end def print_string

# end class Minfloat

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
    args = cmd.parse_args (argv)
    pg = Minfloat (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
