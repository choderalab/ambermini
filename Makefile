include config.h

PROGS=antechamber acdoctor am1bcc antechamber_pvt atomtype bondtype charmmgen \
		database espgen parmcal parmchk parmchk2 prepgen residuegen sqm teLeap \
		tleap translate parmchk_pvt parmchk2_pvt

all: prep libs
	$(MAKE) antechamber
	$(MAKE) leap
	$(MAKE) sqm

install: all
	mkdir -p $(PREFIX)/bin $(PREFIX)/share
	cd bin && /bin/mv $(PROGS) $(PREFIX)/bin
	cp -r share/amber/dat $(PREFIX)

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
	cd blas && $(MAKE) install
	cd lapack && $(MAKE) install
	cd arpack && $(MAKE) install

clean:
	-cd antechamber && $(MAKE) clean
	-cd leap && $(MAKE) clean
	-cd sqm && $(MAKE) clean
	-cd lib && $(MAKE) clean
	-cd cifparse && $(MAKE) clean

uninstall:
	cd $(PREFIX)/bin && rm -f $(PROGS)
	rm -fr $(PREFIX)/dat
