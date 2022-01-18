PGAPy: Python Wrapper for PGAPack Parallel Genetic Algorithm Library
====================================================================

.. |--| unicode:: U+2013   .. en dash

:Author: Ralf Schlatterbeck <rsc@runtux.com>

News 01-2022: This version wraps multiple evaluation with NSGA-III (note
the additional 'I').

News 12-2021: This version wraps multiple evaluation values from your
objective function: Now you can return more than one value to either use
it for constraints (that must be fulfilled before the objective is
optimized) or for multi-objective optimization with the Nondominated
Sorting Genetic Algorithm V.2 (NSGA-II). You can combine both,
multi-objective optimization and constraints.

News: This version wraps the Differential Evolution method (that's quite
an old method but is newly implemented in pgapack).

PGAPy is a wrapper for PGAPack, the parallel genetic algorithm library
(see `PGAPack Readme`_), a powerfull genetic algorithm library by
D. Levine, Mathematics and Computer Science Division Argonne National
Laboratory. The library is written in C. PGAPy wraps this library for
use with Python. The original PGAPack library is already quite old but
is one of the most complete and accurate (and fast, although this is not
my major concern when wrapping it to python) genetic algorithm
implementations out there with a lot of bells and whistles for
experimentation. It also has shown a remarkably small number of bugs
over the years. It supports parallel execution via the message
passing interface MPI_ in addition to a normal "serial" version. That's
why I wanted to use it in Python, too.

To get started you need the PGAPack library, although
it now comes bundled with PGApy, to install a *parallel* version you
currently need a pre-installed PGApack compiled for the MPI library of
choice. See `Installation`_ section for details.

There currently is not much documentation for PGAPy.
You really, absolutely need to read the documentation that comes
with PGAPack.
The PGAPack user guide is now shipped together with PGAPy. It is
installed together with some examples in share/pgapy, wherever the
Python installer decides to place this in the final installation
(usually ``/usr/local/share`` on Linux).

The original PGAPack library can still be downloaded from the PGAPack_
ftp site, it is written in ANSI C and therefore *should* run on most
platforms. Note that this version is not very actively maintained. I've
started a `PGAPack fork on github`_ where I've ported the library to the
latest version of the MPI_ standard and have fixed some minor
inconsistencies in the documentation. I've also implemented some new
features (notably enhancements in selection schemes and a new replacement
strategy called *restricted tournament replacement* and, more recently,
the differential evolution strategy.)

Note: When using NSGA_III replacement for multi (or many-) objective
optimization you need to either

- set reference points on the hyperplane intersecting all axes at
  offset 1. These reference points can be obtained with the convenience
  function ``pga.das_dennis``, it creates a regular set of reference points
  using an algorithm originally publised by I. Das and J. E. Dennis.
  These points are then passed as the paramter ``reference_points`` to
  the ``PGA`` constructor.

  See ``examples/dtlz2.py`` for a usage example and the user guide for
  the bibliographic reference. The function gets the dimensionality of
  the objective space (``num_eval`` minus ``num_constraint``) and the
  number of partition to use.
- Or set reference directions (in the objective space) with the
  ``reference_directions`` parameter, number of partitions for these
  directions with the ``refdir_partitions`` parameter (see
  ``das_dennis`` above, this uses Das/Dennis points internally), and a
  scale factor with the parameter ``refdir_scale``.

You can set both, these parameters are not mutually exclusive.

I have tested pgapy on Linux only and I'll currently not provide Windows
versions.  You also can find my `PGAPack fork on github`_ this
repository has the three upstream releases as versions in git and
contains some updates concerning support of newer MPI_ versions and
documentation updates.  I've also included patches in the git repository
of the Debian maintainer of the package, Dirk Eddelbuettel.

.. _`PGAPack Readme`:
   https://github.com/schlatterbeck/pgapack/blob/master/README.rst
.. _PGAPack:          http://ftp.mcs.anl.gov/pub/pgapack/
.. _`PGAPack fork on github`: https://github.com/schlatterbeck/pgapack
.. _MPI: http://mpi-forum.org/
.. _`my pgapack debian package builder`:
    https://github.com/schlatterbeck/debian-pgapack

To get you started, I've included some very simple examples in
``examples``, e.g., ``one-max.py`` implements the "Maxbit" example
similar to one in the PGAPack documentation. The examples were inspired
by the book "Genetic Algorithms in Python" but are written from scratch
and don't include any code from the book. The examples illustrates
several points:

- Your class implementing the genetic algorithm needs to inherit from
  pga.PGA (pga is the PGAPy wrapper module).
- You need to define an evaluation function called ``evaluate`` that
  returns a sequence of numbers indicating the fitness of the gene given.
  It gets the parameters ``p`` and ``pop`` that can be used to fetch allele
  values from the gene using the ``get_allele`` method, for more details
  refer to the PGAPack documentation. The number of evaluations returned
  by your function is defined with the constructor parameter
  ``num_eval``, the default for this parameter is 1. If your evaluation
  function does not return multiple evaluations (with the default
  setting of ``num_eval``) you can either return a one-element sequence
  or a single return value.
- When using multiple evaluations, these can either be used for
  constraints (the default) or for multi-objective optimization. In the
  latter case, the number of constraints (which by default is one less
  than the number of evaluations set with the parameter ``num_eval``)
  must be set to a number that leaves at least two evaluations for
  objectives. The number of constraints can be set with the parameter
  ``num_constraint``. When using multi-objective optimization, you need
  one of the two replacement-types ``PGA_POPREPL_NSGA_II`` or
  ``PGA_POPREPL_NSGA_III``, set this with the ``pop_replace_type`` parameter.
- You *can* define additional functions overriding built-in functions
  of the PGAPack library, illustrated by the example of
  ``print_string``.  Note that we could call the original print_string
  method of our PGA superclass.  In the same way you can implement,
  e.g., your own crossover method.
- The constructor of the class needs to define the Gene type, in the
  examples we use int and bool built-in datatypes.
- The length of the gene needs to be given in the constructor.
- We often want to maximize the numbers returned by our evaluation
  function, set the parameter ``maximize`` to False if you want to
  minimize.
- For non-binary genes we can define an array of init values, each entry
  containing a sequence with lower and upper bound. The array has to
  have the length of the gene. Note that the upper bound is *included*
  in the range of possible values (unlike the python range operator but
  compatible with the PGAPack definition).
- In the constructor of the class we can add parameters of the genetic
  algorithm. Not all parameters of PGAPack are wrapped yet, currently
  you would need to consult the sourcecode of PGAPy to find out which
  parameters are wrapped. In the example we define several print
  options.
- Finally the genetic algorithm is started with the ``run`` method.

Naming conventions in PGAPy
---------------------------

When you extend PGAPy |--| remember not all functions of PGAPack are
wrapped yet and you may need additional functions |--| you should stick to
my naming conventions when making changes.
The following naming conventions were used for the wrapper:

- Constants of PGAPack like ``PGA_REPORT_STRING`` are used as-is in
  uppercase. These constants can be directly imported from the wrapper
  module. Not all constants are wrapped so far, if you need more, add
  them to the constdef array in pgamodule.c and send_ me a patch.
- For methods of the pga.PGA class I've removed the ``PGA`` prefix used
  throughout PGAPack and converted the method to lowercase with
  underscores between uppercase words in the original function name, so
  ``PGARun`` becomes ``run``, ``PGACheckStoppingConditions`` becomes
  ``check_stopping_conditions``. An exception of the lowercase-rule is
  whenever a name contains "GA" (for "genetic algorithm"), So
  ``PGASetMaxGAIterValue`` becomes ``max_GA_iter``.
- Where possible I've made a single class method where PGAPack needs a
  separate function for each datatype, so ``PGAGetBinaryAllele``,
  ``PGAGetCharacterAllele``, ``PGAGetIntegerAllele``, ``PGAGetRealAllele`` all
  become ``get_allele``. Same holds true for ``set_allele``.
- Whenever a name in PGApack has a "Value" or "Flag" suffix, I've left
  this out, so ``PGAGetFitnessCmaxValue`` becomes ``fitness_cmax``
  and ``PGAGetMutationAndCrossoverFlag`` becomes
  ``mutation_and_crossover``, the only exception to this rule is for the
  two functions ``PGAGetMutationRealValue`` and
  ``PGAGetMutationIntegerValue`` which become ``mutation_value`` not
  just ``mutation``.
- Some fields can take multiple values (they are implemented by ORing
  integer constants, in python they are specified as a list or tuple of
  constants). These are converted to plural (if not already plural in
  PGApack), e.g., ``PGASetStoppingRuleType`` becomes ``stopping_rule_types``.
- Internal method names in the wrapper program have a leading PGA\_ |--| so
  the class method ``set_allele`` is implemented by the C-function
  ``PGA_set_allele`` in ``pgamodule.c``.

Constructor Parameters
----------------------

PGApack has a lot of ``PGASet`` and ``PGAGet`` functions for setting
parameters. These are reflected in constructor parameters on the one hand
and in read-only properties of a ``PGA`` object on the other hand. The
following table gives an overview of all the original PGApack names and
the names of the python wrapper. For the PGApack name I've only listed
the ``PGASet`` function, in many cases there is a corresponding
``PGAGet`` function. If a corresponding read-only property exists for a
constructor parameter this is indicated in the "Prop" column. In some
cases properties are missing because no corresponding ``PGAGet`` function
is implemented in PGApack, in other cases returning a numeric value that
has a symbolic constant in PGApy doesn't make much sense.
The properties have the same name as the constructor parameter.
There are Properties that don't have a corresponding constructor
parameter, namely the ``eval_count`` property (returning the count of
function evaluations) and the
``GA_iter`` property that returns the current GA generation. In the type
column I'm listing the Python type. If the type is followed by a number,
more than one item of that type is specified (a sequence in Python). Some
entries contain "sym", these are integer values with a symbolic constant,
the value "msym" indicates that several values denoted by a list of
symbolic constants can be given. A special case are the
``PGASetRealInitRange``, ``PGASetRealInitPercent``,
``PGASetIntegerInitRange`` functions. These take two values for *each
allele* of the gene. In python this is a sequence of 2-tuples.
Note that this means that you can have different ranges of allowed values
for each allele.

The ``num_eval`` property is special: Due to limitations of the C
programming language, for multiple evaluations in C the first evaluation
is returned as the function return-value of the ``evaluate`` function
and all other parameters are returned in an auxiliary array. PGApack
specifies the number of auxiliary evaluations to be returned. In Python
the evaluation function can always return a sequence of evaluation
values and the ``num_eval`` is one more than ``PGAGetNumAuxEval`` would
return. The default for ``num_eval`` is 1.

The first two (mandatory) constructor parameters are the type of the gene
(this takes a Python type, e.g., ``bool`` for a binary genome or ``int``
for an integer genome) and the length. Note that the ``string_length`` is
implicitly set with the ``length`` parameter. The ``string_length`` is
also available as the length of the ``PGA`` object using the Python
built-in ``len`` function.

==================================== ================================= ====== ====
PGApack name                         Constructor parameter             Type   Prop
==================================== ================================= ====== ====
``PGASetCrossoverBoundedFlag``       ``crossover_bounded``             int    yes
``PGASetCrossoverBounceBackFlag``    ``crossover_bounce_back``         int    yes
``PGASetCrossoverSBXEta``            ``crossover_SBX_eta``             float  yes
``PGASetCrossoverSBXOncePerString``  ``crossover_SBX_once_per_string`` int    yes
``PGASetCrossoverProb``              ``crossover_prob``                float  yes
``PGASetCrossoverType``              ``crossover_type``                sym    no
``PGAGetEvalCount``                  ``eval_count``                    int    yes
``PGASetFitnessCmaxValue``           ``fitness_cmax``                  float  yes
``PGASetFitnessType``                ``fitness_type``                  sym    no
``PGAGetGAIterValue``                ``GA_iter``                       int    yes
``PGASetIntegerInitPermute``         ``integer_init_permute``          int2   no
``PGASetIntegerInitRange``           ``init``                                 no
``PGASetMaxFitnessRank``             ``max_fitness_rank``              float  yes
``PGASetMaxGAIterValue``             ``max_GA_iter``                   int    yes
``PGASetMaxNoChangeValue``           ``max_no_change``                 int    no
``PGASetMaxSimilarityValue``         ``max_similarity``                int    no
``PGASetMutationAndCrossoverFlag``   ``mutation_and_crossover``        int    yes
``PGASetMutationBounceBackFlag``     ``mutation_bounce_back``          int    yes
``PGASetMutationBoundedFlag``        ``mutation_bounded``              int    yes
``PGASetMutationIntegerValue``       ``mutation_value``                int    yes
``PGASetMutationOrCrossoverFlag``    ``mutation_or_crossover``         int    yes
``PGASetMutationPolyEta``            ``mutation_poly_eta``             float  yes
``PGASetMutationPolyValue``          ``mutation_poly_value``           float  yes
``PGASetMutationProb``               ``mutation_prob``                 float  yes
``PGASetMutationRealValue``          ``mutation_value``                float  yes
``PGASetMutationType``               ``mutation_type``                 sym    no
``PGASetNoDuplicatesFlag``           ``no_duplicates``                 int    no
``PGASetNumAuxEval``                 ``num_eval``                      int    yes
``PGASetNumConstraint``              ``num_constraint``                int    yes
``PGASetNumReplaceValue``            ``num_replace``                   int    yes
``PGASetPopSize``                    ``pop_size``                      int    yes
``PGASetPopReplaceType``             ``pop_replace_type``              sym    no
``PGASetPrintFrequencyValue``        ``print_frequency``               int    yes
``PGASetPrintOptions``               ``print_options``                 msym   no
``PGASetPTournamentProb``            ``p_tournament_prob``             float  yes
``PGASetRandomizeSelect``            ``randomize_select``              int    yes
``PGASetRandomSeed``                 ``random_seed``                   int    yes
``PGASetRealInitRange``              ``init``                                 no
``PGASetRealInitPercent``            ``init_percent``                         no
``PGASetRestartFlag``                ``restart``                       int    yes
``PGASetRestartFrequencyValue``      ``restart_frequency``             int    yes
``PGASetRTRWindowSize``              ``rtr_window_size``               int    yes
``PGASetSelectType``                 ``select_type``                   sym    no
``PGASetStoppingRuleType``           ``stopping_rule_types``           msym   no
``PGASetStringLength``               ``string_length``                 int    yes
``PGASetSumConstraintsFlag``         ``sum_constraints``               int    yes
``PGASetTournamentSize``             ``tournament_size``               int    yes
``PGASetTournamentWithReplacement``  ``tournament_with_replacement``   int    yes
``PGASetTruncationProportion``       ``truncation_proportion``         float  yes
``PGASetUniformCrossoverProb``       ``uniform_crossover_prob``        float  yes
==================================== ================================= ====== ====

PGA Object Methods
------------------

The following are the methods that can be used during the run of the
genetic search. The ``run`` method is used to start the search. This can
be used, to, e.g., set an allele during hill-climbing in a custom
``endofgen`` method. Note that some methods only apply to certain gene
types, e.g. the ``encode_int_`` methods can only be used on binary
alleles (they encode an integer value as a binary or gray code
representation into the gene). Other methods take or return different
types depending on the type of gene, e.g. ``get_allele`` or
``set_allele``, they call different backend functions depending on the
gene type. With the ``set_random_seed`` method, the random number
generator can be re-seeded. It is usually best to seed the generator
once at (before) the beginning by specifying ``random_seed`` in the
constructor. For further details consult the user guide.
The method ``get_evaluation`` will return a double for a single
evaluation and a tuple of double for multiple evaluations (when num_eval
is >1)

============================= ================== ===========================
Method                        Parameters         Return
============================= ================== ===========================
``check_stopping_conditions``                    True if stop should occur
``encode_int_as_binary``      *p, pop,*          None
                              *frm, to, val*
``encode_int_as_gray_code``   *p, pop,*          None
                              *frm, to, val*
``encode_real_as_binary``     *p, pop, frm, to*  None
                              *l, u, val*
``encode_real_as_gray_code``  *p, pop, frm, to*  None
                              *l, u, val*
``euclidian_distance``        *p1, pop1*         float
                              *p2, pop2*
``fitness``                   *pop*              None
``get_allele``                *p, pop, index*    allele value
``get_best_index``            *pop*              index of best string
``get_best_report_index``     *pop, idx*         index of best eval with idx
``get_evaluation``            *p, pop*           evaluation of *p*
``get_evaluation_up_to_date`` *p, pop*           True if up-to-date
``get_fitness``               *p, pop*           fitness of *p* (float)
``get_int_from_binary``       *p, pop, frm, to*  int
``get_int_from_gray_code``    *p, pop, frm, to*  int
``get_iteration``                                deprecated, use ``GA_iter``
``get_real_from_binary``      *p, pop,*          float
                              *frm, to, l, u*
``get_real_from_gray_code``   *p, pop,*          float
                              *frm, to, l, u*
``random01``                                     float between 0 and 1
``random_flip``               *probability*      0 or 1
``random_gaussian``           *mean, stddev*     float
``random_interval``           *l, r*             int between l, r
``random_uniform``            *l, r*             float between l, r
``run``                                          None
``select_next_index``         *pop*              index selected individual
``set_allele``                *p, pop, i, value* None
``set_evaluation``            *p, pop, value*    None
``set_evaluation_up_to_date`` *p, pop, status*   None
``set_random_seed``           *seed*             None (use constructor!)
============================= ================== ===========================

User-Methods
------------

PGApack has the concept of user functions. These allow customization of
different areas of a genetic algorihm. In Python they are implemented as
methods that can be changed in a derived class. One of the methods that
*must* be implemented in a derived class is the ``evaluate`` function
(although technically it is not a user function in PGApack). It
interprets the gene and returns an evaluation value or a sequence of
evaluation values if you set the ``num_eval`` constructor parameter.
PGApack computes a fitness from the raw evaluation value. For some
methods an up-call into the PGA class is possible, for some methods this
is not possible (and in most cases not reasonable). Note that for the
``stop_cond`` method, the standard check for stopping conditions can be
called with::

  self.check_stopping_conditions()

The following table lists the overridable methods with their parameters
(for the function signature the first parameter *self* is omitted). Note
that in PGApack there are additional user functions that are needed for
user-defined data types which are currently not exposed in Python. In the
function signatures *p* denotes the index of the individual and *pop*
denotes the population. If more than one individual is specified (e.g.,
for crossover) these can be followed by a number. For crossover *c1* and
*c2* denote the destination individuals (children). The *propability* for
the mutation method is a floating-point value between 0 and 1. Remember
to count the number of mutations that happen, and return that value for
the mutation method!

=================== ============================== ================= =======
Method              Call Signature                 Return Value      Up-Call
=================== ============================== ================= =======
``check_duplicate`` *p1, pop1, p2, pop2*           True if dupe      no
``stop_cond``                                      True to stop      no
``crossover``       *p1, p2, p_pop, c1, c2, c_pop* None              no
``endofgen``                                       None              no
``evaluate``        *p, pop*                       sequence of float no
``gene_difference`` *p1, pop1, p2, pop2*           float             no
``initstring``      *p, pop*                       None              no
``mutation``        *p, pop, propability*          #mutations        no
``pre_eval``        *pop*                          None              no
``print_string``    *file, p, pop*                 None              yes
=================== ============================== ================= =======

Constants
---------

The following PGApack constants are available:

========================== ===========================================
Constant                   Description
========================== ===========================================
PGA_CROSSOVER_ONEPT        One-point Crossover
PGA_CROSSOVER_SBX          Simulated Binary Crossover
PGA_CROSSOVER_TWOPT        Two-point Crossover
PGA_CROSSOVER_UNIFORM      Uniform Crossover
PGA_FITNESSMIN_CMAX        Map fitness by subtracting worst
PGA_FITNESSMIN_RECIPROCAL  Map fitness via reciprocal
PGA_FITNESS_NORMAL         Linear normalization of fitness
PGA_FITNESS_RANKING        Linear fitness ranking
PGA_FITNESS_RAW            Identity fitness function
PGA_MUTATION_CONSTANT      Mutation by adding/subtracting constant
PGA_MUTATION_GAUSSIAN      Mutation by selecting from Gaussian distribution
PGA_MUTATION_PERMUTE       Mutation swaps two random genes
PGA_MUTATION_POLY          Polynomial Mutation
PGA_MUTATION_RANGE         Replace gene with uniform selection from init range
PGA_MUTATION_UNIFORM       Mutation uniform from interval
PGA_NEWPOP                 Symbolic constant for new population
PGA_OLDPOP                 Symbolic constant for old population
PGA_POPREPL_BEST           Population replacement best strings
PGA_POPREPL_NSGA_II        Use NSGA-II replacement for multi-objective opt.
PGA_POPREPL_NSGA_III       Use NSGA-III replacement for multi-objective opt.
PGA_POPREPL_PAIRWISE_BEST  Compare same index in old and new population
PGA_POPREPL_RANDOM_NOREP   Population replacement random no replacement
PGA_POPREPL_RANDOM_REP     Population replacement random with replacement
PGA_POPREPL_RTR            Restricted Tournament Replacement
PGA_REPORT_AVERAGE         Report average evaluation
PGA_REPORT_HAMMING         Report hamming distance
PGA_REPORT_OFFLINE         Report offline
PGA_REPORT_ONLINE          Report online
PGA_REPORT_STRING          Report the string
PGA_REPORT_WORST           Report the worst evaluation
PGA_SELECT_LINEAR          Return individuals in population order
PGA_SELECT_PROPORTIONAL    Fitness-proportional selection
PGA_SELECT_PTOURNAMENT     Binary probabilistic tournament selection
PGA_SELECT_SUS             Stochastic universal selection
PGA_SELECT_TOURNAMENT      Tournament selection
PGA_SELECT_TRUNCATION      Truncation selection
PGA_STOP_MAXITER           Stop on max iterations
PGA_STOP_NOCHANGE          Stop on max number of generations no change
PGA_STOP_TOOSIMILAR        Stop when individuals too similar
========================== ===========================================


Missing Features
----------------

As already mentioned, not all functions and constants of PGAPack are
wrapped yet |--| still for many applications the given set should be
enough. If you need additional functions, you may want to wrap these and
send_ me a patch.

Another feature of PGAPack is currently not implemented in the wrapper,
the usage of custom datatypes. With PGAPack you can define your own
datatypes complete with their custom implementations of the genetic
algorithm functionality like crossover, mutation, etc. I don't expect
problems implementing these, though.

Reporting Bugs
--------------

Please use the `Sourceforge Bug Tracker`_  or the `Github Bug Tracker`_ and

- give a short description of what you think is the correct behaviour
- give a description of the observed behaviour
- tell me exactly what you did.
- if you can publish your source code this makes it a lot easier to
  debug for me

.. _`Sourceforge Bug Tracker`:
    http://sourceforge.net/tracker/?group_id=152022&atid=782852
.. _`Github Bug Tracker`:
    https://github.com/schlatterbeck/pgapy/issues
.. _send: mailto:rsc@runtux.com

Resources
---------

Project information and download from `Sourceforge main page`_

.. _`Sourceforge main page`: http://sourceforge.net/projects/pgapy/

or checkout from Github_

.. _`Github`: http://github.com/schlatterbeck/pgapy

or directly install via pypi.

Installation
------------

PGApy, as the name suggests, supports parallelizing the evaluation
function of the genetic algorithm. This uses the Message Passing
Interface (MPI_) standard.

To install a *serial* version (without parallel programming using MPI_)
you can simply install from pypi using ``pip``. Alternatively when you
have unpacked or checked out from sources you can install with::

 python3 setup.py install --prefix=/usr/local

If you want a parallel version using an MPI_ (Message-Passing Interface)
library you will have to install a parallel version of PGApack first.
The easiest way to do this is to use `my pgapack debian package builder`_
from github. Clone this repository, check out the branch ``debian/sid``,
install the build dependencies, they're listed in the file
``debian/control`` and build the debian packages using::

  dpkg-buildpackage -rfakeroot

This builds pgapack debian packages for *all* supported MPI libraries in
debian, currently these are ``mpich``, ``openmpi``, and ``lam``. In addition
to the MPI libraries a serial version of the pgapack library is also
built. Proceed by installing the package pgapack and the MPI backend
library of choice. If you don't have a preference for an MPI library,
``libpgapack-openmpi`` is the package that uses the Debians default
preferences of an MPI library.

Once a parallel version of PGApack is installed, you can install PGApy
as follows: You set environment variables for the ``PGA_PARALLEL_VARIANT``
(one of ``mpich``, ``openmpi``, or ``lam``) and set the ``PGA_MODULE`` to
``module_from_parallel_install``. Finally you envoke the setup, e.g.::

 export PGA_PARALLEL_VARIANT=openmpi
 export PGA_MODULE=module_from_parallel_install
 python3 setup.py install --prefix=/usr/local

If your MPI library is installed in a different place you should study
the *Extension* configurations in ``setup.py`` to come up with an
Extension definition that fits your installation. If your installation
is interesting to more people, feel free to submit a patch that adds
your Extension-configuration to the standard ``setup.py``.

Changes
-------

Version 1.2: Many-objective optimization with NSGA-III

- Implement NSGA-III

Version 1.1.6: Polynomial mutation and simulated binary crossover (SBX)

- Simulated binary crossover (SBX)
- Polynomial mutation

Version 1.1.1-1.1.5: Small PGAPack updates, fixes for non-debian

- Fix setup.py for non-debian systems
- Update to latest PGAPack with small changes

Version 1.1: Add multi-objective optimization with NSGA-II

- Wrap latest pgapack version 1.4
- This add multi-objective optimization using the Nondominated Sorting
  Genetic Algorithm version 2 (NSGA-II) by Deb et. al. This makes use of
  the previously-introduced option to return more than one value in the
  objective function. To use the feature you need to set the
  num_constraint parameter to a value that leave some of the function
  values returned by your evaluation function as objective function
  values (and not as constraints). See example in examples/multi.py.

Version 1.0: Add constraint handling

- Wrap latest pgapack version 1.3
- This adds auxiliary evaluations. Now your evaluation function can
  return *multiple* floating-point values as a sequence if you set the
  num_eval paramter >1 in the constructor. Currently additional
  evaluation values are used for constraint handling. Constraint values
  are minimized.  Once they reach zero or a negative value they no
  longer count: The sum of all positive constraints is the overall
  constraint violation.  For details see paper by Deb, 2000, see user
  guide for citation. If you're not using constraints, nothing in your
  code needs changes.
- This release may change the path an optimization takes. So for the
  same seed of the random number generator you will get a different
  result, at least if during the search there are individuals with the
  same evaluation (and different genetic material). This is due to a
  change of the sort function in pgapack (it switched to a stable sort
  from the C standard library).

Version 0.9: Allow installation of parallel version

- Pass argv (or sys.argv) to PGACreate
- Add a stanza to setup.py to allow a parallel installation with a given
  pgapack variant compiled for an MPI library. This currently needs a
  pre-installed pgapack debian package.

Version 0.8: Bugfix in real mutation

- Fix a core-dump in the latest pgapack

Version 0.7: Major changes in wrapping

- Now Differential Evolution is implemented, see the minfloat example
  and the user guide of pgapack.

Version 0.6: Major changes in wrapping

- Now the wrapping uses the standard Python recommendations on how to
  create a custom class.
- Update documentation
- Rename ``fitness_cmax`` (from ``fitness_cmax_value``)
- Better error checking of parameters

Version 0.5: Bug-fix release

- Now the ``setup.py`` works, previous version had an encoding problem
- Wrap some minor new methods
- Bug-fix in PGAPack truncation selection

Version 0.4: Bundle PGAPack

- The PGAPack package is now included as a git submodule. By default we
  build against this library
- License fixes: The module long shipped a ``COPYING`` file that includes
  the 2-clause BSD license. But the headers of ``setup.py`` and ``pgamodule.c``
  still included another license. This has been corrected.

Version 0.3: Feature enhancements, Bug fixes

Port to Python3, Python2 is still supported, license change.

- C-Code of wrapper updated to support both, Python2 and Python3
- Update documentation
- Fix some memory leaks that could result when errors occurred during
  some callback methods
- License change: We now have the 2-clause BSD license (similar to the
  MPICH license of PGAPack), this used to be LGPL.

Version 0.2: Feature enhancements, Bug fixes

64 bit support, more PGAPack functions and attributes wrapped,
Readme-update: Sourceforge logo, Changes chapter.

- Bug-fixes for 64 bit architectures
- More functions and attributes of PGAPack wrapped
- Add a build-rule to ``setup.py`` to allow building for standard-install
  of PGAPack |--| this currently needs editing of ``setup.py`` |--| should use
  autodetect here but this would require that I set up a machine with
  standard install of PGAPack for testing.
- Add Sourceforge logo as required
- Add Changes chapter for automagic releases
- Add the ``__module__`` string to class ``PGA`` in module ``pga``. Now
  calling:: ``help (pga)`` in python works as expected, previously no
  help-text was given for the included module

Version 0.1: Initial freshmeat announcement

PGAPy is a wrapper for PGAPack, the parallel genetic algorithm library,
a powerful genetic algorithm library. PGAPy wraps this library for use
with Python. Pgapack is one of the most complete and accurate genetic
algorithm implementations out there with a lot of features for
experimentation.

- Initial Release
