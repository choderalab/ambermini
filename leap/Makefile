include ../config.h

ifeq ($(OS),Windows_NT)
TLEAP_WRAPPER=tleap.bat
else
TLEAP_WRAPPER=tleap
endif


install:
	@echo "Compiling basic programs"
	cd src/leap && $(MAKE) installteLeap
	@echo "Copying shell script" 
	cp $(TLEAP_WRAPPER) $(BINDIR)/$(TLEAP_WRAPPER)
	chmod +x $(BINDIR)/$(TLEAP_WRAPPER)
	$(MAKE) $(MAKE_XLEAP)

uninstall:
	cd src/leap && $(MAKE) uninstallteLeap
	-rm -f $(BINDIR)/tleap
	-cd src/leap && $(MAKE) uninstall
	-rm -f $(BINDIR)/xleap

install_xleap:
	cd src/Wc && $(MAKE)
	cd src/Xmu && $(MAKE)
	cd src/Xpm && $(MAKE)
	cd src/Xraw && $(MAKE)
	cd src/leap && $(MAKE) install
	@echo "Copying shell script" 
	cp xleap $(BINDIR)/xleap
	chmod +x $(BINDIR)/xleap
	@echo "Copying pixmaps"
	-test -d $(DATDIR)/pixmap || mkdir $(DATDIR)/pixmap
	cd pixmap; /bin/cp * $(DATDIR)/pixmap
	@echo "Finished with the install of LEaP"

skip_xleap:
	@echo "NOTE: skipping the building of xleap as X11 libraries are not installed or"
	@echo "NOTE: -nox11 was specified as an option to configure."

clean:
	cd src/Wc && $(MAKE) clean
	cd src/Xmu && $(MAKE) clean
	cd src/Xpm && $(MAKE) clean
	cd src/Xraw && $(MAKE) clean
	cd src/leap && $(MAKE) clean

