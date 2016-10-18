include config.h

test.antechamber:
	-cd antechamber/bondtype && ./Run.bondtype
	-cd antechamber/c60 && ./Run.c60
	-cd antechamber/charmm && ./Run.charmm
	-cd antechamber/chemokine && ./Run.1b2t
	-cd antechamber/tp && ./Run.tp
	-cd antechamber/tp && ./Run.tp2
	-cd antechamber/ash && ./Run.ash
	# Disabled as grms_tol flag is not working with -ek flag
	# -cd antechamber/sustiva && ./Run.sustiva
	-cd antechamber/fluorescein && ./Run.fluorescein
	-cd antechamber/guanine_amber && ./Run.guanine_amber
	# Disabled as resp is not built
	# -cd antechamber/fpph && ./Run.fpph_resp
	-cd antechamber/residuegen && ./Run.ala
	-cd antechamber/parmchk2 && ./Run.parmchk2

test.leap:
	cd leap && $(MAKE) -k test

test.sqm:
	-cd ./sqm/c2h2 && ./Run
	-cd ./sqm/ch3f && ./Run
	-cd ./sqm/nma && ./Run
	-cd ./sqm/znh+ && ./Run
	-cd ./sqm/h2s && ./Run
	-cd ./sqm/dhplus && ./Run
	-cd ./sqm/fecp2 && ./Run
	-cd ./sqm/external_charges && ./Run

finished:
	@echo ""
	@echo "Finished test suite for AmberTools at `date`."
	@echo ""

is_amberhome_defined:
	@(if [ -z "$(AMBERHOME)" ] ; then \
	    echo "Error: AMBERHOME should be defined or else some tests will fail !" ;\
	    exit 2 ;\
	fi ;\
	)

clean:
	-./test_at_clean.sh
	-cd mmpbsa_py && $(MAKE) -k clean
	-cd ../src/mm_pbsa/Examples && ./Run.mmpbsa.clean
	-cd ../src/FEW/examples/test && ./Run.few.clean
	-cd cpptraj && $(MAKE) -k clean
	-cd leap && $(MAKE) -k clean
	-cd leap/glycam/06j && ./Run.glycam clean
	-cd leap/glycam/06EPb && ./Run.glycam clean
	-cd leap/glycam/06j_10 && ./Run.glycam clean
	-cd leap/glycam/06j_12SB && ./Run.glycam clean

check:
	(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	)

