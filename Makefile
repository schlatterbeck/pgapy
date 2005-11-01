PKG=pgapy
SRC=Makefile setup.py pgamodule.c test.py* README README.html default.css

VERSION=Version.h Version.py
LASTRELASE:=$(shell if x=`lastrelease -d` ;then echo $$x ;else echo 'NO_TAG' ;fi)

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

default.css: ../../html/stylesheets/default.css
	cp ../../html/stylesheets/default.css .

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
