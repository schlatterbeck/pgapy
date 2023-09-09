Examples
========

This directory contains various examples accumulated over the years.
Most of these are also tested in the test directory with known good
outputs.

The ``neural.py`` example tries to optimize neural network weights using
a genetic algorithm. The default problem is 'Xor' and the file was
previously named ``xor.py``. The previous incarnation used the fann_
neural network library but since it lacks some methods in the python
wrapper (notably we cannot read back the weights) I've switched this to
the much more mainstream `scikit.learn`_. This example now also has
another problem: Optimizing a two bit adder. This network has four
inputs (for two 2-bit variables) and three outputs (the sum of two 2-bit
numbers can be at most decimal 6 represented by 3 bits).

It is well known that optimizing the xor problem with genetic algorithms
is easy, see, e.g. [1]_, [2]_. The adder problem is difficult (when
using an architecture that can easily be trained with backpropagation),
see also [1]_ and [2]_. In fact so far (with a population size of 200
and differential evolution) I've found only a single solution with::

 neural.py -R 1 --diff --problem=Adder_Full --pop-size=200

.. [1] Darrell Whitley. The GENITOR algorithm and selection pressure:
       Why rank-based allocation of reproductive trials is best. In
       Schaffer [3]_, pages 116–121.
.. [2] Darrell Whitley and Thomas Hanson. Optimizing neural networks
       using faster, more accurate genetic search. In Schaffer [3]_,
       pages 391–396.
.. [3] J. David Schaffer, editor. Proceedings of the Third International
       Conference on Genetic Algorithms (ICGA). Morgan Kaufmann,
       Fairfax, Virginia, June 1989.

.. _fann: http://leenissen.dk/fann/wp/
.. _`scikit.learn`: https://scikit-learn.org/
