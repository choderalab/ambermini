include config.h

PROGS=antechamber 

all: prep $(PROGS)
	$(MAKE) antechamber
	$(MAKE) leap
	$(MAKE) sqm

install: all
	cd bin && /bin/mv $(PROGS) $(PREFIX)/bin
	cp -r share/amber $(PREFIX)/share

prep:
	mkdir -p bin lib include

antechamber::
	cd antechamber && $(MAKE) install

leap::
	cd leap && $(MAKE) install

sqm::
	cd sqm && $(MAKE) install
