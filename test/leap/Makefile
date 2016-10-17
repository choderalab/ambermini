# Makefile for tleap tests

test.ol3:
	@-cd OL3 && ./Run.OL3.sh $(OPT)

test.yil:
	@-cd Yildirim && ./Run.yildirim.sh $(OPT)
	@-cd Yildirim && ./Run.GACU.sh $(OPT)

test.ff14sb.ol3.bsc0:
	@-cd 4PY5 && ./Run.4PY5.sh $(OPT)

test.bsc0:
	@-cd BSC0 && ./Run.BSC0.sh $(OPT)

test.bsc1:
	@-cd BSC1 && ./Run.BSC1.sh $(OPT)

test.ol15:
	@-cd OL15 && ./Run.OL15.sh $(OPT)

test.pbradii:
	@-cd pbradii && ./Run.pbradii.sh $(OPT)

test.leap:
	@-./Run.tleap
	@-./Run.parm7
	@-./Run.cmds
	@-./Run.mbondi3
	@-./Run.addions
	@-cd glycam/06j && ./Run.glycam
	@-cd glycam/06EPb && ./Run.glycam
	@-cd glycam/06j_10 && ./Run.glycam
	@-cd glycam/06j_12SB && ./Run.glycam


COMPLETE_TESTS=test.ol3 \
               test.yil \
               test.ff14sb.ol3.bsc0 \
               test.bsc0 \
               test.bsc1 \
               test.ol15 \
               test.pbradii

test: test.all

test.all:
	$(MAKE) test.complete test.leap

test.complete: $(COMPLETE_TESTS)

clean:
	@echo "Cleaning all test directories:"
	(find . -name 'leap.log' -o -name 'leap.out' | \
	 while read dif ;\
	 do \
	     rm -f $$dif ;\
	 done ;\
	)
	$(MAKE) test.complete OPT=clean

