PKG=pgapy
SRC=Makefile setup.py pgamodule.c test.py*

VERSION=Version.h
LASTRELASE:=$(shell if x=`lastrelease -d` ;then echo $$x ;else echo 'NO_TAG' ;fi)

all: $(VERSION)

$(VERSION): $(SRC)

dist: all
	python setup.py sdist

%.py: %.v
	sed -e 's/RELEASE/$(LASTRELASE)/' $< > $@

%.h: %.py
	python $< > $@

clean:
	rm -f MANIFEST Version.h Version.py Version.pyc
	rm -rf dist build
