include config.h

ifeq ($(OS),Windows_NT)
WRAPPER_SFX=.bat
else
WRAPPER_SFX=
endif

PROGS=antechamber$(WRAPPER_SFX) acdoctor$(SFX) am1bcc$(SFX) antechamber_pvt$(SFX) \
		atomtype$(SFX) bondtype$(SFX) charmmgen$(SFX) database$(SFX) espgen$(SFX) \
		match$(SFX) match_atomname$(SFX) parmcal$(SFX) parmchk$(WRAPPER_SFX) \
		paramfit$(WRAPPER_SFX) parmchk2$(WRAPPER_SFX) prepgen$(SFX) residuegen$(SFX) respgen$(SFX) sqm$(SFX) \
		teLeap$(SFX) tleap$(WRAPPER_SFX) translate$(SFX)  \
		paramfit_pvt$(SFX) \
		parmchk_pvt$(SFX) parmchk2_pvt$(SFX)

all: prep libs
	$(MAKE) antechamber
	$(MAKE) leap
	$(MAKE) paramfit
	$(MAKE) sqm

install: all
	$(MKDIR) -p $(PREFIX)/bin $(PREFIX)/share
	cd bin && mv $(PROGS) $(PREFIX)/bin
	cp -r share/amber/dat $(PREFIX)
	cp -r test $(PREFIX)
	cp config.h $(PREFIX)/test

prep:
	$(MKDIR) -p bin

antechamber::
	cd antechamber && $(MAKE) install

leap::
	cd leap && $(MAKE) install

paramfit::
	cd paramfit && $(MAKE) install

sqm::
	cd sqm && $(MAKE) install

libs::
	cd cifparse && $(MAKE) install
	cd blas && $(MAKE) install
	cd lapack && $(MAKE) install
	cd arpack && $(MAKE) install

clean:
	cd antechamber && $(MAKE) clean
	cd leap && $(MAKE) clean
	cd sqm && $(MAKE) clean
	cd lib && $(MAKE) clean
	cd cifparse && $(MAKE) clean

uninstall:
	cd $(PREFIX)/bin && rm -f $(PROGS)
	$(RM) -fr $(PREFIX)/dat
