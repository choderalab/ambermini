logfile ff94_prm.log

set default OldPrmtopFormat on
source oldff/leaprc.ff94

logfile ff94_prm.log

x = loadpdb ff94/all_aminoan.p
saveamberparm x ./all_aminoan94.top ./all_aminoan94.crd

ncres = { NALA CALA NPRO CPRO }
x = loadpdbusingseq ff94/all_aminonc.p ncres
saveamberparm x ./all_aminonc94.top ./all_aminonc94.crd

x = loadpdb ff94/all_dna94.p
saveamberparm x ./all_dna94.top ./all_dna94.crd
strand = { RC5 RG RU RA3 }
x = loadpdbusingseq ff94/all_rna94.p strand
saveamberparm x ./all_rna94.top ./all_rna94.crd

quit
