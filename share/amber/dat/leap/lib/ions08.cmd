clearVariables
logFile ions08.log
#
#	Monovalent, monoatomic ions using more rational atom types
#

i = createAtom   Li+  Li+  1.0
set i    element Li
set i    position { 0 0 0 }
r = createResidue Li+
add r i
Li+ = createUnit Li+
add Li+ r
saveOff Li+ ./ions08.lib

i = createAtom   Na+  Na+  1.0
set i    element Na
set i    position { 0 0 0 }
r = createResidue Na+
add r i
Na+ = createUnit Na+
add Na+ r
saveOff Na+ ./ions08.lib

i = createAtom   K+   K+   1.0
set i    element K
set i    position { 0 0 0 }
r = createResidue K+
add r i
K+ = createUnit K+
add K+ r
saveOff K+ ./ions08.lib

i = createAtom   Rb+  Rb+  1.0
set i    element Rb
set i    position { 0 0 0 }
r = createResidue Rb+
add r i
Rb+ = createUnit Rb+
add Rb+ r
saveOff Rb+ ./ions08.lib

i = createAtom   Cs+  Cs+  1.0
set i    element Cs
set i    position { 0 0 0 }
r = createResidue Cs+
add r i
Cs+ = createUnit Cs+
add Cs+ r
saveOff Cs+ ./ions08.lib

i = createAtom   F-  F-  -1.0
set i    element F
set i    position { 0 0 0 }
r = createResidue F-
add r i
F- = createUnit F-
add F- r
saveOff F- ./ions08.lib

i = createAtom   Cl-  Cl-  -1.0
set i    element Cl
set i    position { 0 0 0 }
r = createResidue Cl-
add r i
Cl- = createUnit Cl-
add Cl- r
saveOff Cl- ./ions08.lib

i = createAtom   Br-  Br-  -1.0
set i    element Br
set i    position { 0 0 0 }
r = createResidue Br-
add r i
Br- = createUnit Br-
add Br- r
saveOff Br- ./ions08.lib

i = createAtom   I-  I-  -1.0
set i    element I
set i    position { 0 0 0 }
r = createResidue I-
add r i
I- = createUnit I-
add I- r
saveOff I- ./ions08.lib

i = createAtom   Mg+  MG  2.0
set i    element Mg
set i    position { 0 0 0 }
r = createResidue Mg+
add r i
Mg+ = createUnit Mg+
add Mg+ r
saveOff Mg+ ./ions08.lib

quit
