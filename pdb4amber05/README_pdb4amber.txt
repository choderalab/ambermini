COPYRIGHT:

#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials
#       provided with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical
#       Research Inc.nor the names of its contributors may be used to
#       endorse or promote 
#       products derived from this software without specific prior
#       written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

INSTALLATION:

(1) Copy the Python routine pdb4amber to a directory in you PATH (or
remember where it is and call it later with PATH/pdb4amber).
(2) Make it executable (chmod ugo+x PATH/pdb4amber)


PURPOSE:

pdb4amber analyses PDB files and cleans them for further usage,
especially with the LeaP programs of Amber. It does NOT use any
information in the original PDB file other than that contained in ATOM
and HETATM records. The final output files are stripped of everything
not directly related to ATOM or HETATM records. 

RUNNING:

Type pdb4amber or pdb4amber -h to get further instructions:

--------------------------------------------------------
 pdb4amber version 0.5 (December 2013)
--------------------------------------------------------
Usage: pdb4amber [options]

Options:
  -h, --help           show this help message and exit
  -i FILE, --in=FILE   PDB input file                      (no default)
  -o FILE, --out=FILE  PDB output file                     (no default)
  -y, --nohyd          remove all hydrogen atoms           (default: no)
  -d, --dry            remove all water molecules          (default: no)
  -p, --prot           keep only Amber-compatible residues (default: no)


-i or --in specifies the PDB input file (required)
-o or --out specifies the name of the new PDB output file (required)

NOTE: Input and output file names cannot be the same. Be aware that on
some systems with non case-sensitive formats, the automatic check of
this requirement might fail if input and output files just differ in
upper- and lower-case. 

-y or --nohyd requires that all hydrogens are removed (NOT default)
   Usually PDB standard PDB files from x-ray have no hydrogens, but
   already processed PDB files might have hydrogens attached. Since
   LeaP is strict on hydrogen atom names, it is better to remove all
   hydrogens prior to running a PDB file through LeaP.

-d or --dry removes all water molecules (NOT default)
   Water molecules in PDB files are useful in many occasions. However,
   you might want to remove them prior to explicit solvent simulations
   or when running 3D-RISM. Removed waters are stored in a separate
   file by pdb4amber (see under OUTPUT below)

-p or --prot keeps only the actual protein + water + recognized CAPS
   like ACE, NME, and NHE (NOT default)
   It is advised to separate pure protein and ligands prior to
   preparing files for Amber simulations. With this option, anything
   not pure protein is removed by pdb4amber and stored in a separate
   file for further processing later (see under OUTPUT below). Without
   this option, everything is kept also in the newly generated PDB
   file. 

In general, you should use the '-ypd' combination of options to
generate a protein-only file.

OUTPUT:

The new output file (specified with -o or --out) is a standard PDB
file with all residues sequentially re-numbered from 1 to N. In
addition, several other files are created automatically:

(a) A text file with the output PDB file name and _renum.txt added. This
    is a table to help convert the renumbered residues into the original
    ones.
(b) A PDB file with the output PDB file name and _nonprot.pdb appended.
    This is a PDB file that contains only non-protein residues (apart
    from water), i.e., mainly ligands and other stuff.
(c) When using -d (--dry), a PDB file with the output file name plus
    _water.pdb added. This file contains exclusively the water that
    has been stripped from the original PDB file.
(d) A text file with the output PDB file name and _sslink attached, if
    disulfide bonds have been detected by pdb4amber. This file might
    be used by the 'pytleap' script (part of the 'amberlite' package
    in AmberTools 13) to generate the correct disulfide bonds between
    cysteines. 

The following information is written to the screen (but can also be
captured into a text file by ending the command line with '>', e.g.:
pdb4amber -i pdbin.pdb -o pdbout.pdb [-options] > some_file_name.log

Chains:
-------
All chain indicators in the PDB file are listed. This is useful
especially in cases where the x-ray unit cell contains more than one
image of a protein (or complex). In many cases, one is only interested
in one main peptide chain. A long list of different chains may
indicate that the PDB file should cleaned manually prior to 
using pdb4amber.

Insertions:
-----------
Insertions are mostly 'artificial' residue numbers to keep specific
key residue numbers in large protein families constant. 
pdb4amber discards insertion codes and re-numbers all residues from 1
to N. But the insertions are listed to the screen and also included in
the _renum.txt file. 

Histidines:
-----------
By default, Amber routines assume that all histidines are epsilon
tautomers (HIE). pdb4amber lists all HIS residues (with the final
renumbered residue numbers) to allow users to easily locate HIS and
change some to delta (HID), protonated (HIP), or histidinium bound
to zinc (HIN) if required by the local environment in the protein. The
residue numbers refer to the renumbered scheme!

Non-Standard Residues:
----------------------
Non-standard residues (i.e., residues not automatically recognized by
Amber) are listed. Mostly they are ligands (sometimes co-factors,
detergent, buffer components, etc.). The user must take care of these
separately. These residues are also found in the _nonprot.pdb file
mentioned above. They are removed from the final output PDB file if
the -p (--prot) option was chosen. Otherwise they are left also in the
output PDB file.

Cysteines in Disulfide Bonds:
-----------------------------
pdb4amber locates possible (most probable) disulfide bonds by
checking the distance between SG (gamma sulfur) atoms in cysteines. If
a distance SG-SG less than 2.5 Angstroem is found between the SG atoms
of two CYS, a disulfide bond is assumed. The respective CYS residues
are renamed to CYX (required for Amber) in the final PDB output
file. In that case, another small text file (original file name +
_sslink) is also created. This file can be used with the *pytleap*
utility later to create the correct Amber topology file. The residue
numbers of the CYX residues refer to the renumbered scheme!

Gaps:
-----
pdb4amber tries hard (and mostly succeeds) in locating 'gaps', i.e.,
missing residues in the PDB file. This is done by checking distances
of consecutive C-alpha atoms. If such a distance is larger than 5
Angstroem, pdb4amber considers that there is a gap between the two
residues and reports the gap to the screen. The listed residue numbers
refer to the renumbered scheme! 
It is up to user to decide how to handle the gaps. 
DOING NOTHING AT ALL WILL MOST PROBABLY LEAD TO TROUBLE LATER!
By simply introducing a TER record at the gap, Amber (LeaP) will later
introduce the charged N (NH3+) or C (COO-) terminals at the gap
borders. If far from the binding site, this might be OK (except in
long and unconstrained MD, where such unnatural charges will
inevitably lead to unrealistic behavior).
The better solution is to introduce ACE or NME caps at the correct
positions (in addition to a TER record separating the gap
residues). This can be done in various ways (e.g. with PyMol). The
correct names of the newly introduced residues (ACE or NME) and atoms
(CH3 for the methyl carbon, C, N, O for the others) must be observed!

Missing Atoms:
--------------
pdb4amber tries to determine missing heavy atoms atoms in standard
amino acids and reports these. Residue numbers refer to the renumbered
sequence. Note that this has no implcations on further usage of the
file with LeaP since missing atoms are added automatically anyway. In
some cases, this addition may lead to clashes however and it might be
useful to know which residues are actually affected by LeaP.


FINAL REMARK:

The PDB file format and the non-respect of the (in principle
well-defined rules) is a constant source of trouble. There is no
guarantee that pdb4amber works for all possible variants and format
violations in PDB files. If you encounter problems with this routine,
please report them to me so that I can find a solution. You will never
loose (overwrite) your original PDB file. The worst that can happen is
that the resulting output is flawed and not usable.

Romain M. Wolf
Basel, December 2013
romain.wolf (at) gmail.com
