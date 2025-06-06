# To use this Makefile, get a copy of my Release Tools
# git clone git@github.com:schlatterbeck/releasetool.git
# or from sourceforge:
# git clone git://git.code.sf.net/p/sfreleasetools/code sfreleasetools
# And point the environment variable RELEASETOOLS to the checkout
ifeq (,${RELEASETOOLS})
    RELEASETOOLS=../releasetools
endif

README=README.rst
EXAMPLES=$(patsubst %.py,examples/%.py,cards.py cards_mutate.py    \
    constraint.py dtlz2.py fourbar.py gears.py hello_world_char.py \
    hello_world_int.py himmelblau.py magic_prio.py magic_square.py \
    minfloat.py multi.py namefull.py one_max.py sort_numbers.py    \
    twobar.py vibr.py neural.py)
SRC=Makefile MANIFEST.in setup.py $(README) README.html \
    pgamodule.c $(EXAMPLES)
PGAPACK_DOC=pgapack/docs
USERGUIDE=$(PGAPACK_DOC)/user_guide.pdf

LASTRELEASE:=$(shell $(RELEASETOOLS)/lastrelease -n)
VERSIONH=Version.h
VERSIONPY=pga/Version.py
VERSION=$(VERSIONH) $(VERSIONPY)

USERNAME=schlatterbeck
PROJECT=pgapy
PACKAGE=pgapy
CHANGES=changes
NOTES=notes

all: $(VERSION) $(USERGUIDE)

$(VERSION): $(SRC)

test: all
	python3 setup.py build_ext --inplace
	python3 -m pytest test

clean:
	rm -f MANIFEST $(VERSION) Version.pyc default.css README.html
	rm -rf ${CLEAN} PGAPy.egg-info pga/pga.cpython* pga/__pycache__
	make -C pgapack clobber

$(USERGUIDE):
	make -C $(PGAPACK_DOC)

dist: $(USERGUIDE) $(VERSION)

include $(RELEASETOOLS)/Makefile-pyrelease

# vim: set noet :
