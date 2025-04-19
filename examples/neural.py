#!/usr/bin/python3

# This implements several classic neural network training problems.
# It uses the scikit.learn multilayer perceptron implementation.
# We train a simple neural network to recognize the XOR function.
# Training, of course, uses a genetic algorithm, not backpropagation.

try:
    from sklearn.neural_network import MLPRegressor
    from sklearn.exceptions import ConvergenceWarning
except ImportError:
    MLPRegressor = None
from argparse import ArgumentParser
import pga
import sys
import os
import warnings
import numpy as np
# Silence tensorflow warnings
os.environ ['TF_CPP_MIN_LOG_LEVEL'] = '2'

try:
    import tensorflow as tf
    ly = tf.keras.layers
    Layer = ly.Layer
except ImportError:
    tf = None
    ly = None
    class Layer:
        pass

if MLPRegressor is not None:
    warnings.filterwarnings('ignore', category = ConvergenceWarning)

class Scikitlearn_Mixin:
    """ Use neural network from scikit.learn
    """

    def __init__ (self, *args, **kw):
        self.nn = nn = MLPRegressor \
            ( hidden_layer_sizes = (self.n_hidden,)
            , activation         = 'tanh'
            , max_iter           = 1
            )
        # Seems we can set the size only with fit:
        x = np.array ([[0] * self.n_input])
        y = np.array ([0] * self.n_output)
        # Special case for only one dimension in output
        if len (y) > 1:
            y = [y]
        nn.fit (x, y)
        super ().__init__ (*args, **kw)
    # end def __init__

    def set_coefficients (self, cf1, cf2, b1, b2):
        np.copyto (self.nn.intercepts_ [0], b1)
        np.copyto (self.nn.intercepts_ [1], b2)
        np.copyto (self.nn.coefs_ [0],      cf1)
        np.copyto (self.nn.coefs_ [1],      cf2)
    # end def set_coefficients

    def predict (self, *v):
        p = self.nn.predict (*v)
        if self.n_output == 1:
            return np.array ([p]).T
        return np.array (p)
    # end def predict

# end class Scikitlearn_Mixin

class Keras_Predict:

    def predict (self, *v):
        p = self.nn (*v)
        return p.numpy ()
    # end def predict

# end class Keras_Predict

class Keras_Mixin (Keras_Predict):

    def __init__ (self, *args, **kw):
        dt = tf.float64
        input  = ly.Input (shape = (self.n_input,), dtype = dt)
        hidden = ly.Dense \
            (self.n_hidden, activation='tanh',   dtype = dt) (input)
        output = ly.Dense \
            (self.n_output, activation='linear', dtype = dt) (hidden)
        self.nn = tf.keras.Model (input, output)
        super ().__init__ (*args, **kw)
    # end def __init__

    def set_coefficients (self, cf1, cf2, b1, b2):
        b = [b1,  b2]
        c = [cf1, cf2]
        for n, layer in enumerate (self.nn.layers [1:]):
            layer.set_weights ([c [n], b [n]])
    # end def set_coefficients

# end class Keras_Mixin

class Select_Layer (Layer):

    def __init__ (self, indeces, dtype = None):
        super ().__init__ (dtype = dtype)
        self.units   = len (indeces)
        self.indeces = indeces
        self.offsets = []
        s = 0
        for idx in self.indeces:
            self.offsets.append ([s, s + len (idx)])
            s += len (idx)
    # end def __init__

    def build (self, input_shape):
        #w_init = tf.random_normal_initializer ()
        z_init = tf.zeros_initializer ()
        #print ('input_shape:', input_shape)
        self.w = tf.Variable \
            ( initial_value = z_init
                  (shape = (input_shape [-1], self.units), dtype = self.dtype)
            , trainable = True
            )
        self.b = tf.Variable \
            ( initial_value = z_init (shape = (self.units,), dtype = self.dtype)
            , trainable = True
            )
    # end def build

    def call (self, inputs):
        r = []
        for n, idx in enumerate (self.indeces):
            r.append \
                ( tf.matmul
                      ( tf.gather (inputs, idx, axis = 1)
                      , tf.gather (self.w, idx, axis = 0) [:, n:n+1]
                      )
                + self.b [n]
                )
        return tf.concat (r, 1)
# end class Select_Layer

class Neural_Net_Generic (pga.PGA):
    """ This generalizes a neural network with one hidden layer
    """

    def __init__ (self, args):
        self.args       = args
        self.float_bits = 5
        self.bit_gene   = False
        self.do_stop    = False
        d = dict \
            ( maximize        = False
            , pop_size        = self.args.pop_size
            , num_replace     = self.args.pop_size - 2
            , print_options   = [pga.PGA_REPORT_STRING]
            , random_seed     = args.random_seed
            , max_GA_iter     = args.max_generations
            , print_frequency = args.print_rate
            )
        if self.args.gray_code or self.args.binary:
            self.bit_gene = True
        if getattr (self, 'length', None):
            length = self.length
        else:
            length = \
                ( (self.n_input  + 1) * self.n_hidden
                + (self.n_hidden + 1) * self.n_output
                )
        if not self.bit_gene:
            mw = self.args.max_weight
            d.update (init = [(-mw, mw)] * length)
            d.update (self.de_params ())
        if self.args.output_file:
            d ['output_file'] = args.output_file
        datatype = float
        if self.bit_gene:
            datatype = bool
            length   = length * self.float_bits
        super ().__init__ (datatype, length, **d)
        self.expected  = np.array \
            ([self.function (i) for i in self.input_iter ()])
        self.scaled_in = np.array \
            ([i for i in self.input_iter ()]) * 2 - 1
    # end def __init__

    def de_params (self):
        """ Parameters for differential evolution
        """
        if not self.args.differential_evolution:
            return {}
        variant = self.args.de_variant.upper ()
        variant = getattr (pga, 'PGA_DE_VARIANT_' + variant)
        d = dict \
            ( mutation_type        = pga.PGA_MUTATION_DE
            , mutation_only        = True
            , mutation_bounce_back = True
            , select_type          = pga.PGA_SELECT_LINEAR
            , pop_replace_type     = pga.PGA_POPREPL_PAIRWISE_BEST
            , num_replace          = self.args.pop_size
            , DE_crossover_prob    = self.args.de_crossover_prob
            , DE_jitter            = self.args.de_jitter
            , DE_dither            = self.args.de_dither
            , DE_scale_factor      = self.args.de_scale_factor
            , DE_variant           = variant
            )
        return d
    # end def de_params

    def get_float (self, p, pop, idx):
        if self.bit_gene:
            if self.args.gray_code:
                conv = self.get_real_from_gray_code
            else:
                conv = self.get_real_from_binary
            idx *= self.float_bits
            mw = self.args.max_weight
            v = conv (p, pop, idx, idx + self.float_bits - 1, -mw, mw)
        else:
            v = self.get_allele (p, pop, idx)
        return v
    # end def get_float

    def build_pheno (self, p, pop):
        offset = 0
        cf1 = np.array \
            ([self.get_float (p, pop, i + offset)
              for i in range (self.n_input * self.n_hidden)
            ]).reshape (self.n_input, self.n_hidden)
        offset = self.n_input * self.n_hidden
        b1  = np.array \
            ([self.get_float (p, pop, i + offset)
              for i in range (self.n_hidden)
            ])
        offset += self.n_hidden
        cf2 = np.array \
            ([self.get_float (p, pop, i + offset)
              for i in range (self.n_hidden * self.n_output)
            ]).reshape (self.n_hidden, self.n_output)
        offset += self.n_hidden * self.n_output
        b2  = np.array \
            ([self.get_float (p, pop, i + offset)
              for i in range (self.n_output)
            ])
        self.set_coefficients (cf1, cf2, b1, b2)
    # end def build_pheno

    def evaluate (self, p, pop):
        self.build_pheno (p, pop)
        ev   = self.expected
        av   = self.predict (self.scaled_in)
        avn  = 1 / (1 + np.exp (-av))
        minr = np.logical_and (av >= -0.917, av <= 0.917)
        ev1l = np.where (np.logical_and (ev != 0, av <  -0.917))
        ev0l = np.where (np.logical_and (ev == 0, av >   0.917))
        ev1s = np.logical_and (ev != 0, av >= -0.917)
        ev0s = np.logical_and (ev == 0, av <=  0.917)
        i05  = np.where (minr | ev0s | ev1s)
        s    = 0
        s += np.sum (abs (ev [i05] - avn [i05]) ** 0.5)
        s += np.sum (av [ev1l] ** 2)
        s += np.sum (av [ev0l] ** 2)
        return s
    # end def evaluate

    def input_iter (self):
        """ Interate over all possible input values
        """
        for k in range (2 ** self.n_input):
            inp = []
            for b in range (self.n_input):
                inp.append (int (bool (k & (1 << b))))
            yield inp
    # end def input_iter

    def print_string (self, file, p, pop):
        if self.do_stop:
            print \
                ( "Evals: %s, Gen: %s"
                % (self.eval_count, self.GA_iter)
                , file = file
                )
        self.build_pheno (p, pop)
        for n, inp in enumerate (self.input_iter ()):
            inv = ' '.join (str (i) for i in reversed (inp))
            rv  = ' '.join (str (x) for x in reversed (self.expected [n]))
            inp = self.scaled_in [n]
            v   = self.predict (np.array ([inp])) [0]
            vs  = 1 / (1 + np.exp (-v))
            pv  = ' '.join ('%11.6f' % x for x in reversed (v))
            pvs = ' '.join ('%4.2f' % x for x in reversed (vs))
            print ('%s: %s [%s] [%s]' % (inv, pv, pvs, rv), file = file)
        file.flush ()
        super ().print_string (file, p, pop)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.get_evaluation (best, pga.PGA_OLDPOP)
        if eval <= 1e-10:
            self.do_stop = True
        else:
            self.do_stop = self.check_stopping_conditions ()
        return self.do_stop
    # end def stop_cond

# end class Neural_Net_Generic

class Xor (Neural_Net_Generic):
    n_input  = 2
    n_hidden = 2
    n_output = 1

    def function (self, inp):
        assert len (inp) == 2
        return [inp [0] ^ inp [1]]
    # end def function

# end class Xor

class Adder_Generic:

    def function (self, inp):
        assert len (inp) == 4
        i1 = (inp [1] << 1) + inp [0]
        i2 = (inp [3] << 1) + inp [2]
        v  = i1 + i2
        assert v <= 6
        return [((v & (1 << i)) >> i) for i in range (3)]
    # end def function

# end class Adder_Generic


class Adder_Full (Adder_Generic, Neural_Net_Generic):
    """ A 2-bit adder, fully connected layers, 4 nodes in hidden layer
    """
    n_input  = 4
    n_hidden = 4
    n_output = 3

# end class Adder_Full

class Coder_424 (Neural_Net_Generic):
    """ A 4-to-2 Encoder coupled with a 2-to-4 Decoder, so the network
        needs to learn something like the binary representation.
    """

    n_hidden = 2
    n_input  = 2 ** n_hidden
    n_output = n_input

    def function (self, inp):
        return inp
    # end def function

    def input_iter (self):
        """ Interate over all possible input values
            We only have 1 bit set for each pattern.
        """
        for k in range (self.n_input):
            v = 2 ** k
            yield [((v & (1 << i)) >> i) for i in range (self.n_input)]
    # end def input_iter

# end class Coder_424

class Adder_Sparse (Adder_Generic, Keras_Predict, Neural_Net_Generic):
    """ Sparse neural network for adder
        Originally from Rumelhart, Hinton, Williams 1986 p. 342
        (see README.rst)
    """
    n_input  = 4
    n_hidden = n_output = None
    length   = 18

    def __init__ (self, *args, **kw):
        activation = kw.get ('activation', 'tanh')
        if 'activation' in kw:
            del kw ['activation']
        dt      = tf.float64
        input   = tf.keras.layers.Input (shape = (4,), dtype = dt)
        lh1     = Select_Layer ([[0, 2]], dtype = dt) (input)
        act1    = ly.Activation (activation, dtype = dt) (lh1)
        c1      = ly.concatenate ([input, act1], axis = 1, dtype = dt)
        lh2     = Select_Layer ([[1, 3, 4]], dtype = dt) (c1)
        act2    = ly.Activation (activation, dtype = dt) (lh2)
        c2      = ly.concatenate ([input, act1, act2], axis = 1, dtype = dt)
        output  = Select_Layer ([[0, 2, 4], [1, 3, 4, 5], [5]], dtype = dt) (c2)
        self.nn = tf.keras.Model (input, output)
        super ().__init__ (*args, **kw)
    # end def __init__

    def build_pheno (self, p, pop):
        offset = 0
        c = np.array \
            ( [ self.get_float (p, pop, offset)
              , 0
              , self.get_float (p, pop, offset + 1)
              , 0
              ]
            ).reshape (4, 1)
        offset += 2
        b  = np.array ([self.get_float (p, pop, offset)])
        offset += 1
        self.nn.layers [1].set_weights ([c, b])

        c = np.array \
            ( [ 0
              , self.get_float (p, pop, offset)
              , 0
              , self.get_float (p, pop, offset + 1)
              , self.get_float (p, pop, offset + 2)
              ]
            ).reshape (5, 1)
        offset += 3
        b  = np.array ([self.get_float (p, pop, offset)])
        offset += 1
        self.nn.layers [4].set_weights ([c, b])

        c = np.array \
            ( [ self.get_float (p, pop, offset)
              , 0
              , self.get_float (p, pop, offset + 1)
              , 0
              , self.get_float (p, pop, offset + 2)
              , 0
              #
              , 0
              , self.get_float (p, pop, offset + 3)
              , 0
              , self.get_float (p, pop, offset + 4)
              , self.get_float (p, pop, offset + 5)
              , self.get_float (p, pop, offset + 6)
              #
              , 0
              , 0
              , 0
              , 0
              , 0
              , self.get_float (p, pop, offset + 7)
              ]
            ).reshape (3, 6).T
        offset += 8
        b = np.array ([self.get_float (p, pop, i + offset) for i in range (3)])
        offset += 3
        self.nn.layers [7].set_weights ([c, b])
    # end def build_pheno

# end class Adder_Sparse

def cmd_opt (argv):
    if MLPRegressor is not None:
        warnings.filterwarnings('ignore', category = ConvergenceWarning)
    de_variants = ('best', 'rand')
    # Adder_Sparse needs to be the last one
    problems    = ['Xor', 'Adder_Full', 'Coder_424', 'Adder_Sparse']
    backends    = dict (scikit_learn = Scikitlearn_Mixin, keras = Keras_Mixin)
    if tf is None:
        del backends ['keras']
        del problems [-1]
    if MLPRegressor is None:
        del backends ['scikit_learn']
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( "-b", "--binary"
        , help    = "Use binary code for gene (implies bit-gene)"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "-B", "--backend"
        , help    = "Neuronal network backend to use, default=%(default)s"
        , default = list (sorted (backends)) [-1]
        )
    cmd.add_argument \
        ( "-d", "--differential-evolution"
        , help    = "Use differential evolution (DE), implies float gene"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( "--de-crossover-prob"
        , help    = "Crossover probability for DE, default=%(default)s"
        , type    = float
        , default = 0.9
        )
    cmd.add_argument \
        ( "--de-dither"
        , help    = "Dither for DE, default=%(default)s"
        , type    = float
        , default = 0.2
        )
    cmd.add_argument \
        ( "--de-jitter"
        , help    = "Jitter for DE, default=%(default)s"
        , type    = float
        , default = 0.001
        )
    cmd.add_argument \
        ( "--de-variant"
        , help    = "Variant for DE (rand or best), default=%(default)s"
        , default = 'best'
        )
    cmd.add_argument \
        ( "-F", "--de-scale-factor"
        , help    = "Scale factor F for DE, default=%(default)s"
        , type    = float
        , default = 0.85
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
        ( "--max-weight"
        , help    = "Maximum/minimum weight, default=%(default)s"
        , type    = float
        , default = 127
        )
    cmd.add_argument \
        ( "-O", "--output-file"
        , help    = "Output file for progress information"
        )
    cmd.add_argument \
        ( "-p", "--pop-size"
        , help    = "Population size, default=%(default)s"
        , type    = int
        , default = 30
        )
    cmd.add_argument \
        ( "--print-rate"
        , help    = "Print rate in generations, default=%(default)s"
        , type    = int
        , default = 10
        )
    cmd.add_argument \
        ( "-P", "--problem"
        , help    = "Problem to optimize, default=%(default)s"
        , default = 'Xor'
        )
    cmd.add_argument \
        ( "-R", "--random-seed"
        , help    = "Seed random number generator, default=%(default)s"
        , type    = int
        , default = 42
        )
    args = cmd.parse_args (argv)
    if args.de_variant not in de_variants:
        sys.exit \
            ( "Invalid DE variant: %s, allowed are %s"
            % (args.de_variant, ', '.join (de_variants))
            )
    if args.binary and args.gray_code:
        sys.exit ("Either specify binary or gray-code, not both")
    if (args.binary or args.gray_code) and args.differential_evolution:
        sys.exit ("Binary/gray code cannot be combined with DE")
    if args.problem not in problems:
        sys.exit ("Invalid problem, use one of %s" % ', '.join (problems))
    if args.backend not in backends:
        sys.exit \
            ( "Invalid backend: %s use one of [%s]"
            % (args.backend, ', '.join (backends))
            )
    if args.problem == 'Adder_Sparse':
        Problem = Adder_Sparse
    else:
        class Problem (backends [args.backend], globals () [args.problem]):
            pass
    return args, Problem
# end def cmd_opt

def main (argv = sys.argv [1]):
    args, Problem = cmd_opt (argv)
    pg = Problem (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
