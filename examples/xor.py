#!/usr/bin/python

# This needs the FANN neural network library.
# http://leenissen.dk/fann/wp/
# Also packaged in Debian as python-fann2, python3-fann2
# We train a simple neural network to recognize the XOR function.

from __future__ import print_function
from pga import PGA, PGA_REPORT_STRING
from fann2 import libfann
import sys

class Xor (PGA) :

    def __init__ (self) :
        super (self.__class__, self).__init__ \
            ( float, 9
            , maximize      = False
            , pop_size      = 30
            , num_replace   = 28
            , print_options = [PGA_REPORT_STRING]
            , init          = [(-10.0, 10.0)] * 9
            )
    # end def __init__

    def build_pheno (self, p, pop) :
        ann = libfann.neural_net ()
        #ann.set_activation_function_output (libfann.SIGMOID_SYMMETRIC_STEPWISE)
        #ann.set_activation_function_hidden (libfann.SIGMOID_SYMMETRIC_STEPWISE)
        # connection_rate, (inputs, hidden, outputs)
        ann.create_sparse_array (1, (2, 2, 1))
        # get_num_layers get_num_input get_num_output
        # get_total_connections get_total_neurons
        # init_weights (self, data) print_connections, print_parameters
        # set_weight (self, from_neuron, to_neuron, weight)
        for frm in range (3) :
            for to in range (2) :
                v = self.get_allele (p, pop, frm * 2 + to)
                ann.set_weight (frm, to + 3, v)
                #print (frm, to + 3, ann.get_weight (frm, to + 3))
        for frm in range (3) :
            ann.set_weight (frm + 3, 6, self.get_allele (p, pop, frm + 6))
            #print (frm + 3, 6, ann.get_weight (frm + 3, 6))
        return ann
    # end def build_pheno

    def evaluate (self, p, pop) :
        ann = self.build_pheno (p, pop)
        s = 0
        for i1 in range (2) :
            for i2 in range (2) :
                v = (i1 ^ i2)
                r = ann.run ([i1 * 2 - 1, i2 * 2 - 1]) [0]
                #print (i1, i2, r)
                s += abs (v - r)
        return s
    # end def evaluate

    def print_string (self, file, p, pop) :
        ann = self.build_pheno (p, pop)
        for i1 in range (2) :
            for i2 in range (2) :
                v = ann.run ([i1 * 2 - 1, i2 * 2 - 1])
                print ('%s %s: %s' % (i1, i2, v [0]), file = file)
    # end def print_string

# end class Xor

if __name__ == '__main__' :
    pg = Xor ()
    pg.run ()
