# To use this Makefile, get a copy of my SF Release Tools
# git clone git://git.code.sf.net/p/sfreleasetools/code sfreleasetools
# And point the environment variable RELEASETOOLS to the checkout
ifeq (,${RELEASETOOLS})
    RELEASETOOLS=../releasetools
endif

README=README.rst
EXAMPLES=$(patsubst %.py,examples/%.py,cards.py hello-world.py \
    one-max.py sort-numbers.py)
SRC=Makefile MANIFEST.in setup.py $(README) README.html \
    pgamodule.c $(EXAMPLES)
PGAPACK_DOC=pgapack/docs
USERGUIDE=$(PGAPACK_DOC)/user_guide.pdf

LASTRELEASE:=$(shell $(RELEASETOOLS)/lastrelease -n)
VERSIONH=Version.h
VERSIONPY=Version.py
VERSION=$(VERSIONH) $(VERSIONPY)

USERNAME=schlatterbeck
PROJECT=pgapy
PACKAGE=pgapy
CHANGES=changes
NOTES=notes

all: $(VERSION) $(USERGUIDE)

$(VERSION): $(SRC)

clean:
	rm -f MANIFEST Version.h Version.py Version.pyc default.css README.html
	rm -rf ${CLEAN}
	make -C pgapack clean

$(USERGUIDE):
	make -C $(PGAPACK_DOC)

include $(RELEASETOOLS)/Makefile-pyrelease
