#!/usr/bin/python3

# This needs the FANN neural network library.
# http://leenissen.dk/fann/wp/
# Also packaged in Debian as python-fann2, python3-fann2
# We train a simple neural network to recognize the XOR function.

from fann2    import libfann
from argparse import ArgumentParser
import pga
import sys

class Xor (pga.PGA):

    def __init__ (self, args):
        self.args       = args
        self.float_bits = 5
        self.bit_gene   = False
        d = dict \
            ( maximize      = False
            , pop_size      = 30
            , num_replace   = 28
            , print_options = [pga.PGA_REPORT_STRING]
            , random_seed   = args.random_seed
            , max_GA_iter   = args.max_generations
            )
        if self.args.gray_code or self.args.binary:
            self.bit_gene = True
        if not self.bit_gene:
            d.update (init = [(-10.0, 10.0)] * 9)
        if self.args.output_file:
            d ['output_file'] = args.output_file
        datatype = float
        length   = 9
        if self.bit_gene:
            datatype = bool
            length   = 9 * self.float_bits
        super (self.__class__, self).__init__ (datatype, length, **d)
    # end def __init__

    def get_float (self, p, pop, idx):
        if self.bit_gene:
            if self.args.gray_code:
                conv = self.get_real_from_gray_code
            else:
                conv = self.get_real_from_binary
            idx *= self.float_bits
            v = conv (p, pop, idx, idx + self.float_bits - 1, -10.0, 10.0)
        else:
            v = self.get_allele (p, pop, idx)
        return v
    # end def get_float

    def build_pheno (self, p, pop):
        ann = libfann.neural_net ()
        #ann.set_activation_function_output (libfann.SIGMOID_SYMMETRIC_STEPWISE)
        #ann.set_activation_function_hidden (libfann.SIGMOID_SYMMETRIC_STEPWISE)
        # connection_rate, (inputs, hidden, outputs)
        ann.create_sparse_array (1, (2, 2, 1))
        # get_num_layers get_num_input get_num_output
        # get_total_connections get_total_neurons
        # init_weights (self, data) print_connections, print_parameters
        # set_weight (self, from_neuron, to_neuron, weight)
        for frm in range (3):
            for to in range (2):
                v = self.get_float (p, pop, frm * 2 + to)
                ann.set_weight (frm, to + 3, v)
                #print (frm, to + 3, ann.get_weight (frm, to + 3))
        for frm in range (3):
            ann.set_weight (frm + 3, 6, self.get_float (p, pop, frm + 6))
            #print (frm + 3, 6, ann.get_weight (frm + 3, 6))
        return ann
    # end def build_pheno

    def evaluate (self, p, pop):
        ann = self.build_pheno (p, pop)
        s = 0
        for i1 in range (2):
            for i2 in range (2):
                v = (i1 ^ i2)
                r = ann.run ([i1 * 2 - 1, i2 * 2 - 1]) [0]
                #print (i1, i2, r)
                s += abs (v - r)
        return s
    # end def evaluate

    def print_string (self, file, p, pop):
        ann = self.build_pheno (p, pop)
        for i1 in range (2):
            for i2 in range (2):
                v = ann.run ([i1 * 2 - 1, i2 * 2 - 1])
                print ('%s %s: %s' % (i1, i2, v [0]), file = file)
        #print (ann.get_network_type (), file = file)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
        if eval <= 0:
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Xor

def main (argv):
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-b", "--binary"
        , help    = "Use binary code for gene (implies bit-gene)"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "-g", "--gray-code"
        , help    = "Use gray code for gene (implies bit-gene)"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "-m", "--max-generations"
        , help    = "Maximum number of generation, default=%(default)s"
        , type    = int
        , default = 1000
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
    args = cmd.parse_args (argv)
    if args.binary and args.gray_code:
        print ("Either specify binary or gray-code, not both", file=sys.stderr)
        sys.exit (23)
    pg = Xor (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
