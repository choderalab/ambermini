# makefile for antechamber

include ../config.h

PROGS = \
	$(BINDIR)/am1bcc$(SFX) \
	$(BINDIR)/acdoctor$(SFX) \
	$(BINDIR)/antechamber_pvt$(SFX) \
	$(BINDIR)/antechamber$(SFX) \
	$(BINDIR)/atomtype$(SFX) \
	$(BINDIR)/bondtype$(SFX) \
	$(BINDIR)/charmmgen$(SFX) \
	$(BINDIR)/database$(SFX) \
	$(BINDIR)/espgen$(SFX) \
	$(BINDIR)/parmcal$(SFX) \
	$(BINDIR)/parmchk$(SFX) \
	$(BINDIR)/parmchk_pvt$(SFX) \
	$(BINDIR)/parmchk2$(SFX) \
	$(BINDIR)/parmchk2_pvt$(SFX) \
	$(BINDIR)/prepgen$(SFX) \
	$(BINDIR)/residuegen$(SFX) \
	$(BINDIR)/respgen$(SFX) \
	$(BINDIR)/translate$(SFX)

TESTS = test-penalties test-int-part test-combination

#############################################################################

all: $(PROGS)

tests: $(TESTS)

install: all

clean:
	-rm -f *.o $(TESTS)

uninstall:
	-rm -f $(PROGS)

#############################################################################

#_____________________________________________________________________________
# Dependencies, from ./makedepend:

acdoctor.o: common.h
acdoctor.o: define.h
acdoctor.o: atom.h
acdoctor.o: utility.c
acdoctor.o: common.c
acdoctor.o: ring.c
acdoctor.o: rotate.c
acdoctor.o: ac.c
acdoctor.o: charmm.c
acdoctor.o: mol2.c
acdoctor.o: mopcrt.c
acdoctor.o: divcrt.c
acdoctor.o: mopint.c
acdoctor.o: mopout.c
acdoctor.o: divout.c
acdoctor.o: gcrt.c
acdoctor.o: gzmat.c
acdoctor.o: gout.c
acdoctor.o: pdb.c
acdoctor.o: csd.c
acdoctor.o: mdl.c
acdoctor.o: alc.c
acdoctor.o: hin.c
acdoctor.o: prep.c
acdoctor.o: rst.c
acdoctor.o: jzmat.c
acdoctor.o: jcrt.c
acdoctor.o: jout.c
am1bcc.o: common.h
am1bcc.o: define.h
am1bcc.o: atom.h
am1bcc.o: utility.c
am1bcc.o: common.c
am1bcc.o: rotate.c
am1bcc.o: ac.c
am1bcc.o: pdb.c
antechamber.o: common.h
antechamber.o: define.h
antechamber.o: atom.h
antechamber.o: utility.c
antechamber.o: common.c
antechamber.o: ring.c
antechamber.o: rotate.c
antechamber.o: equatom.c
antechamber.o: ac.c
antechamber.o: charmm.c
antechamber.o: mol2.c
antechamber.o: mmcif.c
antechamber.o: mopcrt.c
antechamber.o: divcrt.c
antechamber.o: sqmcrt.c
antechamber.o: mopint.c
antechamber.o: mopout.c
antechamber.o: divout.c
antechamber.o: gcrt.c
antechamber.o: gzmat.c
antechamber.o: gout.c
antechamber.o: pdb.c
antechamber.o: csd.c
antechamber.o: mdl.c
antechamber.o: alc.c
antechamber.o: hin.c
antechamber.o: prep.c
antechamber.o: rst.c
antechamber.o: jzmat.c
antechamber.o: jcrt.c
antechamber.o: jout.c
antechamber.o: gesp.c
antechamber.o: charge.c
atomtype.o: common.h
atomtype.o: define.h
atomtype.o: atom.h
atomtype.o: utility.c
atomtype.o: common.c
atomtype.o: ring.c
atomtype.o: rotate.c
atomtype.o: ac.c
atomtype.o: mol2.c
bondtype.o: ncsu-penalties.h
bondtype.o: common.h
bondtype.o: define.h
bondtype.o: atom.h
bondtype.o: utility.c
bondtype.o: common.c
bondtype.o: ring.c
bondtype.o: ac.c
bondtype.o: mol2.c
charmmgen.o: common.h
charmmgen.o: define.h
charmmgen.o: atom.h
charmmgen.o: utility.c
charmmgen.o: common.c
charmmgen.o: rotate.c
charmmgen.o: ac.c
charmmgen.o: mol2.c
charmmgen.o: prep.c
crdgrow.o: common.h
crdgrow.o: define.h
crdgrow.o: atom.h
crdgrow.o: utility.c
crdgrow.o: rotate.c
crdgrow.o: common.c
crdgrow.o: pdb.c
crdgrow.o: ac.c
crdgrow.o: prep.c
database.o: utility.c
delphigen.o: common.h
delphigen.o: define.h
delphigen.o: atom.h
delphigen.o: common.c
delphigen.o: ac.c
delphigen.o: pdb.c
espgen.o: common.h
espgen.o: define.h
espgen.o: atom.h
ncsu-penalties.o: 
ncsu-penalties.o: ncsu-penalties.h
parmcal.o: common.h
parmcal.o: define.h
parmcal.o: utility.c
parmchk.o: common.h
parmchk.o: define.h
parmchk.o: atom.h
parmchk.o: utility.c
parmchk.o: common.c
parmchk.o: rotate.c
parmchk.o: ac.c
parmchk.o: mol2.c
parmchk.o: prep.c
parmchk2.o: common.h
parmchk2.o: define.h
parmchk2.o: atom.h
parmchk2.o: utility.c
parmchk2.o: common.c
parmchk2.o: rotate.c
parmchk2.o: ac.c
parmchk2.o: mol2.c
parmchk2.o: prep.c
parmjoin.o: common.h
parmjoin.o: define.h
parmjoin.o: utility.c
prepgen.o: common.h
prepgen.o: define.h
prepgen.o: atom.h
prepgen.o: utility.c
prepgen.o: common.c
prepgen.o: rotate.c
prepgen.o: ring.c
prepgen.o: ac.c
residuegen.o: common.h
residuegen.o: define.h
residuegen.o: atom.h
residuegen.o: utility.c
residuegen.o: common.c
residuegen.o: ac.c
respgen.o: common.h
respgen.o: define.h
respgen.o: atom.h
respgen.o: utility.c
respgen.o: common.c
respgen.o: ac.c
respgen.o: equatom.c
rotate.o: common.h
smiles.o: common.h
smiles.o: define.h
smiles.o: atom.h
smiles.o: globals.h
test-penalties.o: ncsu-penalties.h
translate.o: common.h
translate.o: define.h
translate.o: atom.h
translate.o: utility.c
translate.o: common.c
translate.o: rotate.c
translate.o: ac.c
translate.o: pdb.c
translate.o: mol2.c
translate.o: prep.c
translate.o: lsfit.c

.c.o:
	$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

#_____________________________________________________________________________
# Programs

$(BINDIR)/am1bcc$(SFX): am1bcc.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) am1bcc.o $(LM) -o $(BINDIR)/am1bcc$(SFX)

$(BINDIR)/acdoctor$(SFX): acdoctor.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) acdoctor.o $(LM) -o $(BINDIR)/acdoctor$(SFX)

$(BINDIR)/antechamber_pvt$(SFX): antechamber.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) antechamber.o \
    -L$(LIBDIR) -lcifparse $(LM) -o $(BINDIR)/antechamber_pvt$(SFX)

ifeq ($(OS),Windows_NT)
$(BINDIR)/antechamber$(SFX):
	cp -L antechamber.bat $(BINDIR)
else
$(BINDIR)/antechamber$(SFX):
	cp -L antechamber $(BINDIR)
endif

$(BINDIR)/atomtype$(SFX): atomtype.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) atomtype.o $(LM) -o $(BINDIR)/atomtype$(SFX)

$(BINDIR)/bondtype$(SFX): bondtype.o ncsu-penalties.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) ncsu-penalties.o bondtype.o $(LM) -o $(BINDIR)/bondtype$(SFX)

$(BINDIR)/bondtype.orig.o: bondtype.cpp
	$(CXX) -c $(CXXNOOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) -o bondtype.orig.o bondtype.cpp

$(BINDIR)/bondtype.orig$(SFX): bondtype.orig.o
	$(CXX) $(CXXFLAGS) $(AMBERCXXFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) bondtype.orig.o $(LM) -o $(BINDIR)/bondtype.orig$(SFX)

$(BINDIR)/charmmgen$(SFX): charmmgen.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) charmmgen.o $(LM) -o $(BINDIR)/charmmgen$(SFX)

$(BINDIR)/crdgrow$(SFX): crdgrow.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) crdgrow.o $(LM) -o $(BINDIR)/crdgrow$(SFX)

$(BINDIR)/database$(SFX): database.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) database.o $(LM) -o $(BINDIR)/database$(SFX)

$(BINDIR)/delphigen$(SFX): delphigen.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) delphigen.o -o $(BINDIR)/delphigen$(SFX)

$(BINDIR)/espgen$(SFX): espgen.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) espgen.o $(LM) -o $(BINDIR)/espgen$(SFX)

$(BINDIR)/parmcal$(SFX): parmcal.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) parmcal.o $(LM) -o $(BINDIR)/parmcal$(SFX)

ifeq ($(OS),Windows_NT)
$(BINDIR)/parmchk$(SFX):
	cp -L parmchk.bat $(BINDIR)
else
$(BINDIR)/parmchk$(SFX):
	cp -L parmchk $(BINDIR)
endif

$(BINDIR)/parmchk_pvt$(SFX): parmchk.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) parmchk.o $(LM) -o $(BINDIR)/parmchk_pvt$(SFX)

ifeq ($(OS),Windows_NT)
$(BINDIR)/parmchk2$(SFX):
	cp -L parmchk2.bat $(BINDIR)
else
$(BINDIR)/parmchk2$(SFX):
	cp -L parmchk2 $(BINDIR)
endif

$(BINDIR)/parmchk2_pvt$(SFX): parmchk2.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) parmchk2.o $(LM) -o $(BINDIR)/parmchk2_pvt$(SFX)

$(BINDIR)/parmjoin$(SFX): parmjoin.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) parmjoin.o -o $(BINDIR)/parmjoin$(SFX)

$(BINDIR)/prepgen$(SFX): prepgen.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) prepgen.o $(LM) -o $(BINDIR)/prepgen$(SFX)

$(BINDIR)/residuegen$(SFX): residuegen.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) residuegen.o $(LM) -o $(BINDIR)/residuegen$(SFX)

$(BINDIR)/respgen$(SFX): respgen.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) respgen.o $(LM) -o $(BINDIR)/respgen$(SFX)

$(BINDIR)/translate$(SFX): translate.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) translate.o $(LM) -o $(BINDIR)/translate$(SFX)

$(BINDIR)/top2ff$(SFX): top2ff.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) top2ff.o $(LM) -o $(BINDIR)/top2ff$(SFX)

$(BINDIR)/top2mol2$(SFX): top2mol2.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) top2mol2.o $(LM) -o $(BINDIR)/top2mol2$(SFX)

test-penalties: ncsu-penalties.o test-penalties.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) -o test-penalties $(LDFLAGS) $(AMBERLDFLAGS) ncsu-penalties.o test-penalties.o

test-int-part: ncsu-penalties.h ncsu-penalties.c
	$(CC)  $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o test-int-part ncsu-penalties.c $(LDFLAGS) $(AMBERLDFLAGS) \
   -DENABLE_INT_PART_MAIN=1

test-combination: ncsu-penalties.h ncsu-penalties.c
	$(CC) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o test-combination ncsu-penalties.c $(LDFLAGS) $(AMBERLDFLAGS) \
   -DENABLE_COMBINATION_MAIN=1

# DO NOT DELETE
