PKG=pgapy
SRC=Makefile MANIFEST.in setup.py README README.html \
    pgamodule.c test.py

LASTRELEASE:=$(shell ../svntools/lastrelease -n)
VERSIONH=Version.h
VERSIONPY=Version.py
VERSION=$(VERSIONH) $(VERSIONPY)

USERNAME=schlatterbeck
HOSTNAME=shell.sourceforge.net
PROJECTDIR=/home/groups/p/pg/pgapy/htdocs

all: $(VERSION)

$(VERSION): $(SRC)

dist: all
	python setup.py sdist --formats=gztar,zip
	python setup.py bdist

clean:
	rm -f MANIFEST Version.h Version.py Version.pyc default.css README.html
	rm -rf dist build

include ../make/Makefile-sf
