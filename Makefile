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

LASTRELEASE:=$(shell $(RELEASETOOLS)/lastrelease -n)
VERSIONH=Version.h
VERSIONPY=Version.py
VERSION=$(VERSIONH) $(VERSIONPY)

USERNAME=schlatterbeck
PROJECT=pgapy
PACKAGE=pgapy
CHANGES=changes
NOTES=notes

all: $(VERSION)

$(VERSION): $(SRC)

dist: all
	python setup.py sdist --formats=gztar,zip
	python setup.py bdist

clean:
	rm -f MANIFEST Version.h Version.py Version.pyc default.css README.html
	rm -rf dist build

include $(RELEASETOOLS)/Makefile-sf
