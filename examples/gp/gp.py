#!/usr/bin/python3

from functools import reduce
from operator  import mul
from copy      import deepcopy
import numpy as np
import random

# Example of Genetic Programming (GP)
# [1] John Koza. Genetic Programming: On the Programming of Computers
#     by Means of Natural Selection. MIT Press, Cambridge, January 1992.

class Offspring:

    def __init__ (self, node, random):
        self.node = deepcopy (node)
        assert not isinstance (self, Terminal)
        self.use_terminal = random.random () < 0.1
        if self.use_terminal:
            self.idx       = random.randrange (self.node.n_terminals)
            self.component = self.node.get_terminal (self.idx)
        else:
            self.idx       = random.randrange (self.node.n_func)
            self.component = self.node.get_function (self.idx)
    # end def __init__

    def crossover (self, other):
        # Refuse to be replaced by a terminal
        if not self.component.parent and isinstance (other, Terminal):
            return self.node
        if self.component.parent:
            parent = self.component.parent
            pidx   = parent.child_idx (self.component)
            assert pidx is not None
            # Check maximum depth
            if self.component.d_root + other.depth > self.node.max_depth:
                return self.node
            parent.replace_child (pidx, deepcopy (other))
            node   = self.node
        else:
            node   = deepcopy (other)
            node.parent = None
        #node.invariant (down = True)
        return node
    # end def crossover

# end class Offspring

class Node:
    max_depth = 17

    def __init__ (self):
        self._parent = None
        self.d_root  = 0
    # end def __init__

    @property
    def parent (self):
        return self._parent
    # end def parent

    @parent.setter
    def parent (self, parent):
        self._parent = parent
        if parent is None:
            self.d_root = 0
        else:
            self.d_root = parent.d_root + 1
            assert self in self._parent.children
            self._parent.update_depth ()
    # end def parent

    def invariant (self, down = False):
        if isinstance (self, Terminal):
            return
        if self.n_func != 1 + sum (c.n_func for c in self.children):
            import pdb; pdb.set_trace ()
        if down:
            for c in self.children:
                c.invariant (down = True)
        else:
            if self.parent:
                self.parent.invariant ()
    # end def invariant

    def crossover (self, other, random):
        offspring = [Offspring (self, random), Offspring (other, random)]
        crossed = \
            [ offspring [i].crossover (offspring [1-i].component)
              for i in range (2)
            ]
        return crossed
    # end def crossover

    def get_terminal (self, n):
        assert n < self.n_terminals
        if n == 0:
            if isinstance (self, Terminal):
                return self
        count = 0 # this node is no terminal
        for c in self.children:
            if n < count + c.n_terminals:
                return c.get_terminal (n - count)
            count += c.n_terminals
        assert None
    # end def get_terminal

    def get_function (self, n):
        assert n < self.n_func
        if n == 0:
            assert isinstance (self, Function)
            return self
        count = 1
        for c in self.children:
            if n < count + c.n_func:
                return c.get_function (n - count)
            count += c.n_func
        assert None
    # end def get_function

# end class Node

class Function (Node):
    """ A nonterminal (in the middle of the tree)
        Needs to be subclassed with a function and a name.
        Optionally a new arity can be used.
    """
    arity = 2
    fmt   = None

    def __init__ (self, *children):
        assert len (children) <= self.arity
        if not self.fmt:
            self.fmt = '%s (%s)' % (self.name, ', '.join (('%s',) * self.arity))
        self.children = list (children)
        self.n_terminals = sum (c.n_terminals for c in self.children)
        self.n_func      = sum (c.n_func      for c in self.children) + 1
        self.depth       = 1
        if self.children:
            self.depth += max (c.depth for c in self.children)
        super ().__init__ ()
    # end def __init__

    @property
    def name (self):
        return self.__class__.__name__.split ('_', 1) [-1]
    # end def name

    def add_child (self, child):
        assert len (self.children) < self.arity
        self.children.append (child)
        self.update_n (child.n_terminals, child.n_func)
        child.parent = self
    # end def add_child

    def format (self):
        return self.fmt % tuple (c.format () for c in self.children)
    # end def format
    __str__ = __repr__ = format

    def child_idx (self, child):
        for n, c in enumerate (self.children):
            if c == child:
                return n
        return None
    # end def child_idx

    def replace_child (self, index, child):
        assert index < self.arity
        d_n_terminals = -self.children [index].n_terminals
        d_n_func      = -self.children [index].n_func
        self.children [index] = child
        d_n_terminals += self.children [index].n_terminals
        d_n_func      += self.children [index].n_func
        self.update_n (d_n_terminals, d_n_func)
        child.parent = self
    # end def replace_child

    def update_depth (self):
        depth   = max (c.depth for c in self.children) + 1
        if depth != self.depth:
            self.depth = depth
            if self.parent:
                self.parent.update_depth ()
    # end def update_depth

    def update_n (self, d_n_terminals, d_n_func):
        if d_n_terminals == 0 and d_n_func == 0:
            return
        self.n_terminals += d_n_terminals
        self.n_func      += d_n_func
        if self.parent:
            self.parent.update_n (d_n_terminals, d_n_func)
    # end def update_n

# end def Function

class Terminal (Node):
    """ A Terminal (a leaf of the tree)
        Needs to be instantiated with a name and optionally an index.
    """
    n_terminals = 1
    n_func      = 0
    depth       = 1
    children    = ()
    values      = {}

    def __init__ (self, name, index = None):
        self.name  = name
        self.index = index
        self.value = None
        super ().__init__ ()
    # end def __init__

    def __call__ (self):
        """ A terminal symbol can create a copy of itself (with the same
            parameters). This way it can be used like a class and an
            existing individual can be used to set the current value for
            all copies via the class values dict.
        """
        return self.__class__ (self.name, self.index)
    # end def __call__

    @property
    def value (self):
        return self.values [self.format ()]
    # end def value

    @value.setter
    def value (self, value):
        """ This sets value by key in the class.
            So all nonterminals with same name and index have same value
        """
        self.values [self.format ()] = value
    # end def value

    def format (self):
        if self.index is None:
            return self.name
        return '%s_%s' % (self.name, self.index)
    # end def format
    __str__ = __repr__ = format

    def eval (self):
        return self.value
    # end def eval

# end def Terminal

class F_add (Function):
    fmt = '(%s + %s)'

    def eval (self):
        return sum (c.eval () for c in self.children)
    # end def eval

# end class F_add

class F_sub (Function):
    fmt = '(%s - %s)'

    def eval (self):
        assert len (self.children) == 2
        v = [c.eval () for c in self.children]
        return v [0] - v [1]
    # end def eval

# end class F_sub

class F_mul (Function):
    fmt = '(%s * %s)'

    def eval (self):
        return reduce (mul, (c.eval () for c in  self.children), 1)
    # end def eval

# end class F_mul

class F_div (Function):
    fmt = '(%s / %s)'

    def eval (self):
        """ [1], p.82: Return 1 in case of division by 0
        """
        v = [c.eval () for c in  self.children]
        assert len (v) == 2
        if v [1] == 0:
            return 1.0
        return v [0] / v [1]
    # end def eval

# end class F_div

class F_sin (Function):
    arity = 1

    def eval (self):
        return np.sin (self.children [0].eval ())
    # end def eval
# end class F_sin

class F_cos (Function):
    arity = 1

    def eval (self):
        return np.cos (self.children [0].eval ())
    # end def eval
# end class F_sin

class F_sqrt (Function):
    arity = 1

    def eval (self):
        """ [1] p.83: use sqrt (abs (v))
        """
        assert len (self.children) == 1
        v = self.children [0].eval ()
        return np.sqrt (abs (v))
    # end def eval
# end class F_sqrt

class F_sqrt3 (Function):
    arity = 1

    def eval (self):
        """ [1] p.83: use sqrt (abs (v))
        """
        assert len (self.children) == 1
        v = self.children [0].eval ()
        return abs (v) ** (1/3)
    # end def eval
# end class F_sqrt3

class F_sqrt5 (Function):
    arity = 1

    def eval (self):
        """ [1] p.83: use sqrt (abs (v))
        """
        assert len (self.children) == 1
        v = self.children [0].eval ()
        return abs (v) ** (1/5)
    # end def eval
# end class F_sqrt5

class F_sqrt7 (Function):
    arity = 1

    def eval (self):
        """ [1] p.83: use sqrt (abs (v))
        """
        assert len (self.children) == 1
        v = self.children [0].eval ()
        return abs (v) ** (1/7)
    # end def eval
# end class F_sqrt7

class F_log (Function):
    arity = 1

    def eval (self):
        """ [1] p.83: return 0 if input == 0, log (abs (v)) otherwise
        """
        assert len (self.children) == 1
        v = self.children [0].eval ()
        return np.log (abs (v))
    # end def eval
# end class F_log

class F_and (Function):
    def eval (self):
        assert len (self.children) == 2
        v = [c.eval () for c in self.children]
        return v [0] and v [1]
    # end def eval
# end class F_and

class F_or (Function):
    def eval (self):
        assert len (self.children) == 2
        v = [c.eval () for c in self.children]
        return v [0] or v [1]
    # end def eval
# end class F_or

class F_not (Function):
    arity = 1

    def eval (self):
        assert len (self.children) == 1
        v = self.children [0].eval ()
        return not v
    # end def eval
# end class F_not

class F_nand (Function):
    def eval (self):
        assert len (self.children) == 2
        v = [c.eval () for c in self.children]
        return not (v [0] and v [1])
    # end def eval
# end class F_nand

class F_nor (Function):
    def eval (self):
        assert len (self.children) == 2
        v = [c.eval () for c in self.children]
        return not (v [0] or v [1])
    # end def eval
# end class F_nor

class Genetic_Programming:
    """ Methods used for Genetic Programming, this can be used as a Mixin.
        Note that we asume an instance of a python Random class as
        self.random.
    """

    def __init__ (self, functions, terminals):
        self.functions = functions
        self.terminals = terminals
    # end def __init__

    def grow_tree (self, tree, depth, fulldepth = False):
        flen = len (self.functions)
        tlen = len (self.terminals)
        if depth == 1:
            for i in range (tree.arity):
                idx = self.random.randrange (tlen)
                tree.add_child (self.terminals [idx] ())
        else:
            l = flen
            if not fulldepth:
                l += tlen
            for i in range (tree.arity):
                idx = self.random.randrange (l)
                if idx >= flen:
                    tree.add_child (self.terminals [idx - flen] ())
                else:
                    child = self.functions [idx] ()
                    tree.add_child (child)
                    self.grow_tree (child, depth - 1, fulldepth = fulldepth)
    # end def grow_tree

    def random_tree (self, depth, fulldepth = False):
        """ Return a random tree with maximum depth depth.
            If fulldepth is true, the tree will have maximum depth in
            all branches (we select only nonterminals until the last
            depth level). Otherwise it might not even have the maximum
            depth. The first variant is what [1], p.92f calls the 'full'
            method, the second is what [1] calls the 'grow' method.
        """
        idx = self.random.randrange (len (self.functions))
        tree = self.functions [idx] ()
        self.grow_tree (tree, depth = depth - 1, fulldepth = fulldepth)
        return tree
    # end def random_tree

    def ramped_half_and_half (self, n, depth):
        """ The ramped half-and-half strategy [1] for producing a random
            population of n individuals. This uses a ramped depth from 2
            to depth with equal number of individuals in each ramped
            depth. It creates half of the individuals with
            fulldepth = True and the other half with fulldepth = False.
        """
        pop = []
        k = n // (depth - 1)
        r = n - (k * (depth - 1))
        h = [k//2, k//2 + k - (k // 2) * 2]
        lh = [(k+r) // 2, (k+r) // 2 + (k+r) - (k+r) // 2 * 2]
        assert k
        for d in range (depth - 1):
            ll = h if d != depth - 2 else lh
            for x, nh in enumerate (ll):
                for i in range (nh):
                    pop.append (self.random_tree (d + 2, fulldepth = x))
        return pop
    # end def ramped_half_and_half

# end class Genetic_Programming

if __name__ == '__main__':
    random.seed (23)
    class GP (Genetic_Programming):
        random = random
    # end class GP

    gp = GP ((F_mul, F_div, F_add, F_sub, F_sin), (Terminal ('x'),))
    print (gp.random_tree (4))
    print (gp.random_tree (4, fulldepth = True))
    t = gp.random_tree (5)
    print (t)
    print ( "depth: %s n_term: %s n_fun: %s"
          % (t.depth, t.n_terminals, t.n_func)
          )
    pop = gp.ramped_half_and_half (10, 4)
    for k in pop:
        print (k)
    c1, c2 = pop [-2].crossover (pop [-1], random)
    print (c1)
    print (c2)
