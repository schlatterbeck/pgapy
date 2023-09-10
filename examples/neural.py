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
import warnings
import numpy as np
try:
    import tensorflow as tf
except ImportError:
    tf = None

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
        if self.n_output != 1:
            return p [0]
        return p
    # end def predict

# end class Scikitlearn_Mixin

class Keras_Mixin:

    def __init__ (self, *args, **kw):
        dt = tf.float64
        input  = tf.keras.layers.Input (shape = (self.n_input,), dtype = dt)
        hidden = tf.keras.layers.Dense \
            (self.n_hidden, activation='tanh',   dtype = dt) (input)
        output = tf.keras.layers.Dense \
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

    def predict (self, *v):
        p = self.nn (*v)
        return p.numpy () [0]
    # end def predict
    
# end class Keras_Mixin

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
        s = 0
        for inp in self.input_iter ():
            v   = self.function (inp)
            inp = [x * 2 - 1 for x in inp]
            r1  = self.predict (np.array ([inp]))
            r2  = 1 / (1 + np.exp (-r1))
            for ev, av, avn in zip (v, r1, r2):
                if -0.917 <= av <= 0.917:
                    s += abs (ev - avn) ** 0.5
                elif ev:
                    if av < -0.917:
                        s += av ** 2
                    else:
                        s += abs (ev - avn) ** 0.5
                else:
                    if av > 0.917:
                        s += av ** 2
                    else:
                        s += abs (ev - avn) ** 0.5
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
            print ("Evals: %s" % self.eval_count, file = file)
        self.build_pheno (p, pop)
        for inp in self.input_iter ():
            inv = ' '.join (str (i) for i in reversed (inp))
            rv  = ' '.join (str (x) for x in self.function (inp))
            inp = [x * 2 - 1 for x in inp]
            v   = self.predict (np.array ([inp]))
            vs  = 1 / (1 + np.exp (-v))
            pv  = ' '.join ('%11.6f' % x for x in v)
            pvs = ' '.join ('%4.2f' % x for x in vs)
            print ('%s: %s [%s] [%s]' % (inv, pv, pvs, rv), file = file)
        super ().print_string (file, p, pop)
    # end def print_string

    def stop_cond (self):
        best = self.get_best_index (pga.PGA_OLDPOP)
        eval = self.evaluate (best, pga.PGA_OLDPOP)
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

class Adder_Full (Neural_Net_Generic):
    """ A 2-bit adder, fully connected layers, 4 nodes in hidden layer
    """
    n_input  = 4
    n_hidden = 4
    n_output = 3

    def function (self, inp):
        assert len (inp) == 4
        i1 = (inp [1] << 1) + inp [0]
        i2 = (inp [3] << 1) + inp [2]
        v  = i1 + i2
        assert v <= 6
        return [int (bool (v & 4)), int (bool (v & 2)), v & 1]
    # end def function

# end class Adder_Full

def main (argv):
    if MLPRegressor is not None:
        warnings.filterwarnings('ignore', category = ConvergenceWarning)
    de_variants = ('best', 'rand')
    problems    = ('Xor', 'Adder_Full')
    backends    = dict (scikit_learn = Scikitlearn_Mixin, keras = Keras_Mixin)
    if tf is None:
        del backends ['keras']
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
        , default = 'scikit_learn'
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
        sys.exit ("Invalid problem use one of %s" % ', '.join (problems))
    if args.backend not in backends:
        sys.exit \
            ( "Invalid backend: %s use one of [%s]"
            % (args.backend, ', '.join (backends))
            )
    class Problem (backends [args.backend], globals () [args.problem]):
        pass
    pg = Problem (args)
    pg.run ()
# end def main

if __name__ == '__main__':
    main (sys.argv [1:])
