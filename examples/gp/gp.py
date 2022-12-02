#!/usr/bin/python3

from functools import reduce
from operator  import mul
from copy      import deepcopy
import numpy as np

# Example of Genetic Programming (GP)
# [1] John Koza. Genetic Programming: On the Programming of Computers
#     by Means of Natural Selection. MIT Press, Cambridge, January 1992.

class Offspring:
    max_tries = 5

    def __init__ (self, node, random):
        self.random    = random
        self.node      = deepcopy (node)
        self.max_depth = self.node.max_depth
        self.depth     = self.node.depth
        assert not isinstance (self, Terminal)
        self.compute_component ()
    # end def __init__

    def compute_component (self, depth = None):
        random = self.random
        for i in range (self.max_tries):
            if depth is not None and depth + self.depth > self.max_depth:
                self.component = self.node.get_by_maxdepth (depth, random)
            elif self.depth <= 4:
                self.idx = random.randrange \
                    (self.node.n_terminals + self.node.n_funcs)
                if self.idx >= self.node.n_funcs:
                    self.component = self.node.get_terminal \
                        (self.idx - self.node.n_funcs)
                else:
                    self.component = self.node.get_function (self.idx)
            else:
                # Prefer functions over terminals for deeper trees
                self.use_terminal = random.random () < 0.1
                if self.use_terminal:
                    self.idx       = random.randrange (self.node.n_terminals)
                    self.component = self.node.get_terminal (self.idx)
                else:
                    self.idx       = random.randrange (self.node.n_funcs)
                    self.component = self.node.get_function (self.idx)
            if self.component.depth < self.node.max_depth:
                break
    # end def compute_component

    def crossover (self, other):
        node = self.node
        for i in range (self.max_tries):
            if i:
                self.compute_component (depth = other.depth)
            # Refuse to be replaced by a terminal
            if not self.component.parent and isinstance (other, Terminal):
                continue
            if self.component.parent:
                parent = self.component.parent
                pidx   = parent.child_idx (self.component)
                assert pidx is not None
                # Check maximum depth
                if self.component.d_root + other.depth > self.node.max_depth:
                    continue
                parent.replace_child (pidx, deepcopy (other))
                node   = self.node
            else:
                node   = deepcopy (other)
                node.parent = None
            break
        assert node.depth <= self.max_depth
        if node.debug:
            node.invariant (down = True)
        return node
    # end def crossover

# end class Offspring

class Node:
    max_depth = 17

    def __init__ (self):
        self._parent = None
        self.d_root  = 0
        self.p_eval  = None # eval of parent
        self.evalue  = None # own eval
        self.debug   = False
    # end def __init__

    def __hash__ (self):
        return hash (self.format ())
    # end def __hash__

    @property
    def parent (self):
        return self._parent
    # end def parent

    @parent.setter
    def parent (self, parent):
        self._parent = parent
        self.update_d_root ()
        if parent is not None:
            assert self in self._parent.children
            self._parent.update_depth ()
    # end def parent

    @property
    def symbol (self):
        return getattr (self, 'sym', self.name)
    # end def symbol

    def as_dot (self, first = True):
        """ Output graphviz dot format
        """
        r = []
        if first:
            r.append ('digraph "%s" {' % id (self))
        r.append ('"%s" [label = "%s"];' % (id (self), self.symbol))
        for c in self.children:
            r.append ('"%s" -> "%s";' % (id (self), id (c)))
        for c in self.children:
            r.append (c.as_dot (first = False))
        if first:
            r.append ('}')
        return '\n'.join (r)
    # end def as_dot

    def crossover (self, other, random):
        offspring = [Offspring (self, random), Offspring (other, random)]
        crossed = \
            [ offspring [i].crossover (offspring [1-i].component)
              for i in range (2)
            ]
        for n, c in enumerate (crossed):
            c.p_eval = offspring [n].node.evalue
        return crossed
    # end def crossover

    def _get_by_maxdepth (self, depth, idx):
        if idx == 0:
            assert self.d_root + depth <= self.max_depth
            return self
        idx -= 1
        if self.d_root + depth == self.max_depth:
            return idx
        for c in self.children:
            r = c._get_by_maxdepth (depth, idx)
            if isinstance (r, Node):
                return r
            idx = r
        return idx
    # end def _get_by_maxdepth

    def get_by_maxdepth (self, depth, random):
        if depth == self.max_depth:
            return self
        count = self.get_maxdepth_count (depth)
        assert count
        idx   = random.randrange (count)
        return self._get_by_maxdepth (depth, idx)
    # end def get_by_maxdepth

    def get_function (self, n):
        assert n < self.n_funcs
        if n == 0:
            assert isinstance (self, Function)
            return self
        count = 1
        for c in self.children:
            if n < count + c.n_funcs:
                return c.get_function (n - count)
            count += c.n_funcs
        assert None
    # end def get_function

    def get_maxdepth_count (self, depth):
        """ Return number of nodes with given d_root + depth <= maxdepth
        """
        count = 0
        if self.d_root + depth >= self.max_depth:
            return count
        count += 1
        for c in self.children:
            count += c.get_maxdepth_count (depth)
        return count
    # end def get_maxdepth_count

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

    def invariant (self, down = False):
        if not self.parent and self.d_root != 0:
            #import pdb; pdb.set_trace ()
            assert 0
        if self.parent and self.d_root != self.parent.d_root + 1:
            #import pdb; pdb.set_trace ()
            assert 0
        if isinstance (self, Terminal):
            return
        if self.n_funcs != 1 + sum (c.n_funcs for c in self.children):
            #import pdb; pdb.set_trace ()
            assert 0
        if down:
            for c in self.children:
                c.invariant (down = True)
        else:
            if self.parent:
                self.parent.invariant ()
    # end def invariant

    def update_d_root (self):
        if self.parent is None:
            self.d_root = 0
        else:
            self.d_root = self.parent.d_root + 1
        for c in self.children:
            c.update_d_root ()
    # end def update_d_root

# end class Node

class Function (Node):
    """ A nonterminal (in the middle of the tree)
        Needs to be subclassed with a function and a name.
        Optionally a new arity can be used.
    """
    arity = 2
    fmt   = None

    def __init__ (self, *children, debug = False):
        assert len (children) <= self.arity
        if not self.fmt:
            if getattr (self, 'sym', None) is not None:
                self.fmt = '(%%s %s %%s)' % self.symbol
            else:
                param = ', '.join (('%s',) * self.arity)
                self.fmt = '%s (%s)' % (self.name, param)
        self.children = list (children)
        self.n_terminals = sum (c.n_terminals for c in self.children)
        self.n_funcs     = sum (c.n_funcs     for c in self.children) + 1
        self.depth       = 1
        if self.children:
            self.depth += max (c.depth for c in self.children)
        super ().__init__ ()
        self.debug = debug
    # end def __init__

    @property
    def name (self):
        return self.__class__.__name__.split ('_', 1) [-1]
    # end def name

    def add_child (self, child):
        assert len (self.children) < self.arity
        self.children.append (child)
        self.update_n (child.n_terminals, child.n_funcs)
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
        d_n_funcs     = -self.children [index].n_funcs
        self.children [index] = child
        d_n_terminals += self.children [index].n_terminals
        d_n_funcs     += self.children [index].n_funcs
        self.update_n (d_n_terminals, d_n_funcs)
        child.parent = self
        child.update_d_root ()
    # end def replace_child

    def update_depth (self):
        depth = max (c.depth for c in self.children) + 1
        if depth != self.depth:
            self.depth = depth
            if self.parent:
                self.parent.update_depth ()
    # end def update_depth

    def update_n (self, d_n_terminals, d_n_funcs):
        if d_n_terminals == 0 and d_n_funcs == 0:
            return
        self.n_terminals += d_n_terminals
        self.n_funcs     += d_n_funcs
        if self.parent:
            self.parent.update_n (d_n_terminals, d_n_funcs)
    # end def update_n

# end def Function

class _Terminal (Node):
    """ A Terminal (a leaf of the tree)
        Needs to be instantiated with a name and optionally an index.
    """
    n_terminals = 1
    n_funcs     = 0
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

    def format (self):
        if self.index is None:
            return self.name
        return '%s_%s' % (self.name, self.index)
    # end def format
    __str__ = __repr__ = format

    def eval (self):
        return self.value
    # end def eval

# end def _Terminal

class Terminal (_Terminal):

    @_Terminal.value.setter
    def value (self, value):
        """ This sets value by key in the class.
            So all nonterminals with same name and index have same value
        """
        self.values [self.format ()] = value
    # end def value

# end class Terminal

class Const (Terminal):
    """ A constant """

    def __init__ (self, name, index = None):
        self.name  = name
        assert index is None
        self.index = None
        self.values [self.format ()] = float (name)
    # end def __init__

# end class Const

class F_add (Function):
    sym = '+'

    def eval (self):
        return sum (c.eval () for c in self.children)
    # end def eval

# end class F_add

class F_sub (Function):
    sym = '-'

    def eval (self):
        assert len (self.children) == 2
        v = [c.eval () for c in self.children]
        return v [0] - v [1]
    # end def eval

# end class F_sub

class F_mul (Function):
    sym = '*'

    def eval (self):
        return reduce (mul, (c.eval () for c in  self.children), 1)
    # end def eval

# end class F_mul

class F_div (Function):
    sym = '/'

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

    random_tree_depth = 6

    def __init__ (self, functions, terminals, debug = False):
        self.debug     = debug
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
                    child = self.functions [idx] (debug = self.debug)
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
            depth (except that in lower depth there may not be enough
            different individuals). It creates half of the individuals with
            fulldepth = True and the other half with fulldepth = False.
            The number if filled with trees of the maximum depth with
            fulldepth = False.
        """
        ntries = 5
        treedict = {}
        k = n // (depth - 1)
        r = n - (k * (depth - 1))
        h = [k//2, k//2 + k - (k // 2) * 2]
        lh = [(k+r) // 2, (k+r) // 2 + (k+r) - (k+r) // 2 * 2]
        assert k
        for d in range (depth - 1):
            ll = h if d != depth - 2 else lh
            for x, nh in enumerate (ll):
                for i in range (nh):
                    for tr in range (ntries):
                        t = self.random_tree (d + 2, fulldepth = x)
                        f = t.format ()
                        if f not in treedict:
                            treedict [f] = t
                            break
        miss = 0
        while len (treedict) < n:
            t = self.random_tree (depth, fulldepth = False)
            f = t.format ()
            if f not in treedict:
                treedict [f] = t
            elif miss > 500:
                t = self.random_tree (depth, fulldepth = True)
                f = t.format ()
                if f not in treedict:
                    treedict [f] = t
                else:
                    miss += 1
            else:
                miss += 1
            if miss > 1000:
                depth += 1
                miss   = 0
        return list (treedict.values ())
    # end def ramped_half_and_half

    # methods needed for the GA

    def check_duplicate (self, p1, pop1, p2, pop2):
        g1 = self.get_gene (p1, pop1)
        g2 = self.get_gene (p2, pop2)
        return g1.format () == g2.format ()
    # end def check_duplicate

    def crossover (self, p1_i, p2_i, ppop, c1_i, c2_i, cpop):
        p1 = self.get_gene (p1_i, ppop)
        p2 = self.get_gene (p2_i, ppop)
        if p1.evalue is None:
            assert self.get_evaluation_up_to_date (p1_i, ppop)
            r = self.get_evaluation (p1_i, ppop)
            if np.isscalar (r):
                p1.evalue = r
            else:
                p1.evalue = r [0]
        if p2.evalue is None:
            assert self.get_evaluation_up_to_date (p2_i, ppop)
            r = self.get_evaluation (p2_i, ppop)
            if np.isscalar (r):
                p2.evalue = r
            else:
                p2.evalue = r [0]
        c1, c2 = p1.crossover (p2, self.random)
        if c1.format () == p1.format () or c1.format () == p2.format ():
            c1 = self.mutate (c1)
        if c2.format () == p1.format () or c2.format () == p2.format ():
            c2 = self.mutate (c2)
        self.set_gene (c1_i, cpop, c1)
        self.set_gene (c2_i, cpop, c2)
    # end def crossover

    def initstring (self, p, pop):
        if not self.randpop:
            self.randpop = self.ramped_half_and_half \
                (self.popsize, self.random_tree_depth)
        self.set_gene (p, pop, self.randpop.pop ())
    # end def initstring

    def mutate (self, tree):
        """ We mutate by grafting a random tree into the to-be-mutated
            individual or vice-versa.
        """
        random = self.random
        g = deepcopy (tree)
        t = self.random_tree (self.random_tree_depth, fulldepth = False)
        if g.depth + t.depth > g.max_depth:
            p = g.get_by_maxdepth (g.max_depth - t.depth, random)
            if p.parent:
                cidx = p.parent.child_idx (p)
                p.parent.replace_child (cidx, t)
                mutated = g
            else:
                mutated = t
        else:
            # 50/50 chance of grafting t into g or vice-versa
            if random.random () < 0.5:
                idx = random.randrange (t.n_funcs + t.n_terminals)
                if idx >= t.n_funcs:
                    p = t.get_terminal (idx - t.n_funcs)
                else:
                    p = t.get_function (idx)
                if p.parent:
                    cidx = p.parent.child_idx (p)
                    p.parent.replace_child (cidx, g)
                    mutated = t
                else:
                    mutated = t
            else:
                idx = random.randrange (g.n_funcs + g.n_terminals)
                if idx >= g.n_funcs:
                    p = g.get_terminal (idx - g.n_funcs)
                else:
                    p = g.get_function (idx)
                if p.parent:
                    cidx = p.parent.child_idx (p)
                    p.parent.replace_child (cidx, t)
                    mutated = g
                else:
                    mutated = t
        return mutated
    # end def mutate

    def mutation (self, p, pop, pm):
        """ If we want no duplicates, the GA will repeatedly call
            mutation until a string is found that is not in the
            population.
        """
        if self.random.random () < pm:
            g = self.get_gene (p, pop)
            t = self.mutate (g)
            self.set_gene (p, pop, t)
            return 1
        return 0
    # end def mutation

# end class Genetic_Programming
