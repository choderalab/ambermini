logfile ff94_iowatprm.log

set default OldPrmtopFormat on
source ../../dat/leap/cmd/leaprc

logfile ff94_iowatprm.log

#
#  add ions in vacuum  --test alias too (ai==addions)
#
x = loadpdb ff94/all_dna94.p
ai x IB 0
saveamberparm x ./all_dnaio94.top ./all_dnaio94.crd

#
#  just solvate
#
x = loadpdb ff94/all_dna94.p
alignaxes x
solvatebox x WATBOX216 10
saveamberparm x ./all_dnawat94.top ./all_dnawat94.crd

#
#  add ions/solvate
#
x = loadpdb ff94/all_dna94.p
alignaxes x
addions x Na+ 0
solvatebox x WATBOX216 10
saveamberparm x ./all_dnaiowat94.top ./all_dnaiowat94.crd

#
#  solvate/add ions
#
x = loadpdb ff94/all_dna94.p
alignaxes x
solvatebox x WATBOX216 10
addions x Na+ 0
saveamberparm x ./all_dnawatio94.top ./all_dnawatio94.crd

quit
