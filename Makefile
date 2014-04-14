include config.h

PROGS=antechamber 

all: prep libs
	$(MAKE) antechamber
	$(MAKE) leap
	$(MAKE) sqm

install: all
	cd bin && /bin/mv $(PROGS) $(PREFIX)/bin
	cp -r share/amber $(PREFIX)/share

prep:
	mkdir -p bin

antechamber::
	cd antechamber && $(MAKE) install

leap::
	cd leap && $(MAKE) install

sqm::
	cd sqm && $(MAKE) install

libs::
	cd cifparse && $(MAKE) install

clean:
	-cd antechamber && $(MAKE) clean
	-cd leap && $(MAKE) clean
	-cd sqm && $(MAKE) clean
	-cd lib && $(MAKE) clean
	-cd cifparse && $(MAKE) clean
