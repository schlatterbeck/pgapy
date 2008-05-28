PKG=pgapy
SRC=Makefile MANIFEST.in setup.py README README.html default.css \
    pgamodule.c test.py

VERSION=Version.h Version.py
LASTRELASE:=$(shell ../svntools/lastrelease -n)

USERNAME=schlatterbeck
HOSTNAME=shell.sourceforge.net
PROJECTDIR=/home/groups/p/pg/pgapy/htdocs

all: $(VERSION)

$(VERSION): $(SRC)

dist: all
	python setup.py sdist --formats=gztar,zip
	python setup.py bdist

README.html: README
	rst2html $< > $@

default.css: ../../content/html/stylesheets/default.css
	cp ../../content/html/stylesheets/default.css .

%.py: %.v
	sed -e 's/RELEASE/$(LASTRELASE)/' $< > $@

%.h: %.py
	python $< > $@

upload_homepage: all
	scp README.html $(USERNAME)@$(HOSTNAME):$(PROJECTDIR)/index.html
	scp default.css $(USERNAME)@$(HOSTNAME):$(PROJECTDIR)

clean:
	rm -f MANIFEST Version.h Version.py Version.pyc default.css README.html
	rm -rf dist build

include ../make/Makefile-sf
