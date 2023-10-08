Examples
========

This directory contains various examples accumulated over the years.
Most of these are also tested in the test directory with known good
outputs.

Neural network optimization
---------------------------

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
see also [1]_ and [2]_. In fact (with a population size of 200 and
differential evolution) there are few solutions below 1000 generations.
An example:

 neural.py -R 41 --diff --problem=Adder_Full --pop-size=200

Royal Road function by John Holland
-----------------------------------

This was posed as a challenge to ICGA '93 and later posted to the
Genetic Algorithms mailing list [1]. A better description can be found
in [2]. The original paper proposing Royal Road functions is [3]. It may
be relevant today to check published variants of genetic algorithms that
were tested against this function.

This function is explicitly designed to work better with one-point
crossover than, e.g. two-point crossover: There are gaps of bits that
are unused (so-called introns) between the relevant blocks of bits. This
reduces disruption for one-point crossover because, as long as the
crossover is in one of the intron ranges, no disruption of existing
pattern occurs. For two-point crossover the chance of disruption is much
higher because *both* crossover points need to hit an intron range to
avoid disruption.

For historical reasons the default mixing type in PGAPack is mutation
*or* crossover (so if crossover is performed no mutation occurs). The
Royal Road function needs more mutation because not all relevant pattern
are in the original population (at a population size 512 recommended by
Holland). So the traditional crossover (PGA_MIX_TRADITIONAL) works much
better with Royal Road.

With traditional mixing we can also increase the crossover rate to 1.0
(always), since we are using elitist replacement schemes it makes sense
to do as much crossover as possible.

Since the GA converges quite quickly with Royal Road, a diversity
preserving mechanism can do a lot. With population replacement RTR we
typically find better results. In 100 generations (about 44k function
evaluations) we can often reach the 4th level with only one or two
blocks not converged. When avoiding duplicates (in addition to using
RTR) we can sometimes reach the last level (all ones for the non-intron
bits).

Bibliography
------------

.. [1] Darrell Whitley. The GENITOR algorithm and selection pressure:
       Why rank-based allocation of reproductive trials is best. In
       Schaffer [3]_, pages 116–121.
.. [2] Darrell Whitley and Thomas Hanson. Optimizing neural networks
       using faster, more accurate genetic search. In Schaffer [3]_,
       pages 391–396.
.. [3] J. David Schaffer, editor. Proceedings of the Third International
       Conference on Genetic Algorithms (ICGA). Morgan Kaufmann,
       Fairfax, Virginia, June 1989.
.. [4] D. E. Rumelhart, G. E. Hinton, and R. J. Williams. Learning
       internal representations by error propagation. In Parallel
       Distributed Processing: Explorations in the Microstructure of
       Cognition [5], chapter 8, pages 318–362.
.. [5] David E. Rumelhart and James L. McClelland. Parallel Distributed
       Processing: Explorations in the Microstructure of Cognition,
       volume 1: Foundations. MIT Press, Cambridge, Massachusetts, 1986.
.. [6] John Holland. Royal road functions. In Alan C. Schultz,
       moderator. `Mail archives of the GA-List mailing list`_. Digest
       Volume 7 Issue 22, August 1993.
.. [7] Terry Jones. A description of Holland's Royal Road function.
       Evolutionary Computation, 2(4):409–415, December 1994. `Available
       from ACM`_.
.. [8] Melanie Mitchell, Stephanie Forrest, and John H. Holland.
       The royal road for genetic algorithms: Fitness landscapes
       and GA performance. In Varela and Bourgine [9], pages 245–254.
.. [9] Francisco J. Varela and Paul Bourgine, editors. Toward a
       Practice of Autonomous Systems, Proceedings of the First
       European Conference on Artificial Life. MIT Press, Paris,
       France, December 1991. Proceedings published April 1992.

.. _fann: http://leenissen.dk/fann/wp/
.. _`scikit.learn`: https://scikit-learn.org/
.. _`Mail archives of the GA-List mailing list`:
    http://www.cs.cmu.edu/afs/cs/project/ai-repository/ai/areas/genetic/ga/mail/0.html
.. _`Available from ACM`: https://dl.acm.org/doi/pdf/10.1162/evco.1994.2.4.409
