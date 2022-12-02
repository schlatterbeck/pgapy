Genetic Programming Examples
============================

.. |opt_xor.py| replace:: ``opt_xor.py``
.. |opt_parity3.py| replace:: ``opt_parity3.py``
.. |opt_integral.py| replace:: ``opt_integral.py``

This directory contains examples of Genetic Programming (GP) [1]_, [2]_,
[3]_, [4]_. We do symbolic regression. The example |opt_xor.py|_
optimizes a boolean expression, the goal is to express the exclusive or
operation with only negated and (NAND). Also with only NAND,
|opt_parity3.py|_ searches to express the 3-variable parity
problem. Finally |opt_integral.py|_ does a numeric integration of the
given expression (for 50 numeric values) and tries to find a symbolic
expression that most closely matches the given numbers.

All problems use the difference between the known expected outcome and
the values that are actually computed by the evolved function.
Expressions are represented by a tree (similar to what a compiler would
produce as a parse tree). For the boolean problems it is practical to
compare the output of the expression for *all* combinations in the truth
table of the searched function. For the numeric integration problems a
representative number of values is used. In all cases the sum of
absolute values of the differences between expected and actual outcome
are used as the evaluation function of the genetic algorithm.

Using parse trees for genes is an example of user defined data types in
PGAPy. The latest version of PGAPy supports these user defined data
types.

.. [1] John R. Koza. *Genetic Programming: On the Programming of
   Computers by Means of Natural Selection.* MIT Press, Cambridge,
   January 1992.
.. [2] John R. Koza. *Genetic Programming II: Automatic Discovery
   of Reusable Programs.* MIT Press, Cambridge, 1994.
.. [3] John R. Koza, F. H. Bennett, D. Andre, and M. A. Keane.
   *Genetic Programming III: Darwinian Invention and Problem Solving.*
   Morgan Kaufmann, 1999.
.. [4] John R. Koza, Martin A. Keane, Matthew J. Streeter, William
   Mydlowec, Jessen Yu, and Guido Lanza. *Genetic Programming: Routine
   Human-Competitive Machine Intelligence,* volume IV of Genetic
   Programming. Springer, 2003.
.. _opt_xor.py:
    https://github.com/schlatterbeck/pgapy/blob/master/examples/gp/opt_xor.py
.. _opt_parity3.py:
    https://github.com/schlatterbeck/pgapy/blob/master/examples/gp/opt_parity3.py
.. _opt_integral.py:
    https://github.com/schlatterbeck/pgapy/blob/master/examples/gp/opt_integral.py
