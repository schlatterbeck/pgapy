PGAPy: Python wrapper for pgapack parallel genetic algorithm library
====================================================================

.. |--| unicode:: U+2014   .. em dash

:Author: Ralf Schlatterbeck <rsc@runtux.com>

PGAPy is a wrapper for pgapack, the parallel genetic algorithm library
(see `pgapack Readme`_), a powerfull genetic algorithm library by
D. Levine, Mathematics and Computer Science Division Argonne National
Laboratory. The library is written in C. PGAPy wraps this library for
use with Python. The original pgapack library is already quite old but
is one of the most complete and accurate (and fast, although this is not
my major concern when wrapping it to python) genetic algorithm
implementations out there with a lot of bells and whistles for
experimentation. It also has shown a remarkably small number of bugs
over the years. It supports parallel execution via the message
passing interface MPI_ in addition to a normal "serial" version. That's
why I wanted to use it in Python, too.

There currently is not much documentation for PGAPy.
You really, absolutely need to read the documentation that comes
with pgapack |--| and of course you need the pgapack library.

The original pgapack library can still be downloaded from the pgapack_
ftp site, it is written in ANSI C and therefore *should* run on most
platforms. Note that this version is not very actively maintained. I've
started a `pgapack fork on github`_ where I've ported the library to the
latest version of the MPI_ standard and have fixed some minor
inconsistencies in the documentation. I've also implemented some new
features (notably enhancements in selection schemes and a new replacement
strategy called *restricted tournament replacement*)

I have tested pgapy on Linux only and I'll currently not provide Windows
versions.  You also can find my `pgapack fork on github`_ this
repository has the three upstream releases as versions in git and
contains some updates concerning support of newer MPI_ versions and
documentation updates.  I've also included patches in the git repository
of the Debian maintainer of the package, Dirk Eddelbuettel.

.. _`pgapack Readme`: http://ftp.mcs.anl.gov/pub/pgapack/README
.. _pgapack:          http://ftp.mcs.anl.gov/pub/pgapack/
.. _`pgapack fork on github`: https://github.com/schlatterbeck/pgapack
.. _MPI: http://mpi-forum.org/

To get you started, I've included some very simple examples in
``examples``, e.g., ``one-max.py`` implements the "Maxbit" example
similar to one in the pgapack documentation. The examples were inspired
by the book "Genetic Algorithms in Python" but are written from scratch
and don't include any code from the book. The examples illustrates
several points:

 - Your class implementing the genetic algorithm needs to inherit from
   pga.PGA (pga is the PGAPy wrapper module).
 - You need to define an evaluation function called ``evaluate`` that
   returns a number indicating the fitness of the gene given with the
   parameters ``p`` and ``pop`` that can be used to fetch allele values from
   the gene using the ``get_allele`` method, for more details refer to the
   pgapack documentation.
 - You *can* define additional functions overriding built-in functions
   of the pgapack library, illustrated by the example of
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
   compatible with the pgapack definition).
 - In the constructor of the class we can add parameters of the genetic
   algorithm. Not all parameters of pgapack are wrapped yet, currently
   you would need to consult the sourcecode of PGAPy to find out which
   parameters are wrapped. In the example we define several print
   options.
 - Finally the genetic algorithm is started with the ``run`` method.

Naming conventions in PGAPy
---------------------------

When you extend PGAPy |--| remember not all functions of pgapack are
wrapped yet and you may need additional functions |--| you should stick to
my naming conventions when making changes.
The following naming conventions were used for the wrapper:

 - Constants of pgapack like ``PGA_REPORT_STRING`` are used as-is in
   uppercase. These constants can be directly imported from the wrapper
   module. Not all constants are wrapped so far, if you need more, add
   them to the constdef array in pgamodule.c and send_ me a patch.
 - For methods of the pga.PGA class I've removed the ``PGA`` prefix used
   throughout pgapack and converted the method to lowercase with
   underscores between uppercase words in the original function name, so
   ``PGARun`` becomes ``run``, ``PGACheckStoppingConditions`` becomes
   ``check_stopping_conditions``.
 - Where possible I've made a single class method where pgapack needs a
   separate function for each datatype, so ``PGAGetBinaryAllele``,
   ``PGAGetCharacterAllele``, ``PGAGetIntegerAllele``, ``PGAGetRealAllele`` all
   become ``get_allele``. Same holds true for ``set_allele``.
 - Internal method names in the wrapper program have a leading PGA\_ |--| so
   the class method ``set_allele`` is implemented by the C-function
   ``PGA_set_allele`` in ``pgamodule.c``.

Missing Features
----------------
As already mentioned, not all functions and constants of pgapack are
wrapped yet |--| still for many applications the given set should be
enough. If you need additional functions, you may want to wrap these and
send_ me a patch.

Another feature of pgapack is currently not implemented in the wrapper,
the usage of custom datatypes. With pgapack you can define your own
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

Changes
-------

Version 0.4: Bundle pgapack

  - The pgapack package is now included as a git submodule. By default we
    build against this library
  - License fixes: The module long shipped a COPYING file that includes
    the 2-clause BSD license. But the headers of setup.py and pgamodule.c
    still included another license. This has been corrected.

Version 0.3: Feature enhancements, Bug fixes

Port to Python3, Python2 is still supported, license change.

  - C-Code of wrapper updated to support both, Python2 and Python3
  - Update documentation
  - Fix some memory leaks that could result when errors occurred during
    some callback methods
  - License change: We now have the 2-clause BSD license (similar to the
    MPICH license of pgapack), this used to be LGPL.

Version 0.2: Feature enhancements, Bug fixes

64 bit support, more pgapack functions and attributes wrapped,
Readme-update: Sourceforge logo, Changes chapter.

 - Bug-fixes for 64 bit architectures
 - More functions and attributes of pgapack wrapped
 - Add a build-rule to setup.py to allow building for standard-install
   of pgapack |--| this currently needs editing of setup.py |--| should use
   autodetect here but this would require that I set up a machine with
   standard install of pgapack for testing.
 - Add Sourceforge logo as required
 - Add Changes chapter for automagic releases

Version 0.1: Initial freshmeat announcement

PGAPy is a wrapper for pgapack, the parallel genetic algorithm library,
a powerful genetic algorithm library. PGAPy wraps this library for use
with Python. Pgapack is one of the most complete and accurate genetic
algorithm implementations out there with a lot of features for
experimentation.

 - Initial Release
