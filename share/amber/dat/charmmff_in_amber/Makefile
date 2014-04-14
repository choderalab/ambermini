#
# Script originally written by George Seibel for Amber 3.0 Rev A
# copyright (c) 1989 Regents of the University of California
#
#   Generate the Amber 4.0 standard database
#
PREP=/thr/slater/amber5/exe.R10000/prep

db4.ch22.dat : all22.in all22nt.in all22ct.in dna22.5p.in dna22.3p.in dna22.in
	/bin/rm -f output params
	@echo "Making CHARM-22 compatible database:"
	@echo ""
	@echo "    All-atom residues:"
	$(PREP) -i all22.in
	@/bin/rm -f output params
	@echo "    All-atom charged N-terminal AA residues:"
	$(PREP) -i all22nt.in	
	@/bin/rm -f output params
	@echo "    All-atom charged C-terminal AA residues:"
	$(PREP) -i all22ct.in
	@/bin/rm -f output params
	@echo "    dna residues"
	$(PREP) -i dna22.in
	@/bin/rm -f output params
	@echo "    dna residues at 3p end"
	$(PREP) -i dna22.3p.in
	@/bin/rm -f output params
	@echo "    dna residues at 5p end"
	$(PREP) -i dna22.5p.in
	@/bin/rm -f output params
	@echo done

