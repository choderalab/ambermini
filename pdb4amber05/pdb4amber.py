#!/usr/bin/env python

# Romain M. Wolf, NIBR Basel, December 2013
# with revisions by Pawel Janowski & Jason Swails, Rutgers U., Feb. 2014
#    & Jan. 2015

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
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
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

__version__ = "1.1"
__date__ = "January 2015"

# PDB analyzer to prepare protein(-ligand) PDB files for Amber simulations.

import os, sys
from cStringIO import StringIO
from optparse import OptionParser
from math import sqrt
import subprocess

import signal
def sigint_handler(*args, **kwargs):
  print >> sys.stderr, "Interrupt signal caught. Exiting"
  sys.exit(1)

signal.signal(signal.SIGINT, sigint_handler)

# Global constants
RESPROT = ('ALA', 'ARG', 'ASN', 'ASP',
           'CYS', 'GLN', 'GLU', 'GLY',
           'HIS', 'ILE', 'LEU', 'LYS',
           'MET', 'PHE', 'PRO', 'SER',
           'THR', 'TRP', 'TYR', 'VAL',
           'HID', 'HIE', 'HIN', 'HIP',
           'CYX', 'ASH', 'GLH', 'LYH',
           'ACE', 'NME', 'GL4', 'AS4',
           'WAT','HOH')

#=============================================
def pdb_read(pdbin, noter, model):
#=============================================
# only records starting with the following strings are kept...
  if noter:
    ACCEPTED =  ('ATOM','HETATM')
  else:
    ACCEPTED =  ('ATOM','HETATM', 'TER', 'MODEL', 'ENDMDL')
# records starting with the following strings are considered 'dividers'
# and they are cleaned for the rest of the line...
  DIVIDERS =  ('TER', 'MODEL', 'ENDMDL')

  if pdbin == 'stdin':
    records = sys.stdin
  elif hasattr(pdbin, 'readline'):
    records = pdbin
  else:
    records = open(pdbin, 'r')
  print >> sys.stderr, "\n=================================================="
  print >> sys.stderr, "Summary of pdb4amber for file %s"%pdbin
  print >> sys.stderr, "=================================================="
  '''
  This PDB reader first splits PDB lines (records) into individual (named) fields,
  then later re-assembles the fields in records. This might be done more elegantly,
  but we keep this scheme for now for clarity (given the messy PDB file format)!

    record_type     0    atom_number     1    blank1           2
    atom_name       3    alt_loc_ind     4    residue_name     5
    blank2          6    chain_id        7    residue_number   8
    insertion_code  9    blank3         10    x               11
    y               12   z              13    occupancy       14
    bfactor         15   blank4         16    element         17
    charge          18

  '''

  record_type = []; atom_number = []; blank1 = []
  atom_name = []; alt_loc_ind = []; residue_name = []; blank2 = []
  chain_id = []; residue_number = []; insertion_code = []; blank3 = []
  x = []; y = []; z = []; occupancy = []
  bfactor = []; blank4 = []; element = []; charge = []
# keep everything in 'ACCEPTED'
  lines = records.readlines()

  if model != 0:
    import re
    start, end = -1,-1
    for i, line in enumerate(lines):
      if re.search(r'^MODEL\s+%d\s' %model, line):
        start = i
        break
    for i in range (start, len(lines)):
      if lines[i][0:6]=='ENDMDL':
        end = i
        break
    assert start !=-1, "Requested MODEL %d not found." %model
    assert end !=-1, "ENDMDL line after MODEL %d not found." %model
    lines = lines[start:end]
    # import code; code.interact(local=dict(globals(), **locals()))

  for line in lines:
    # make all lines 80 characters long (fill up with blanks if necessary)
    # so that the line parser will not fail on shorter lines...
    line = (line.rstrip() + (80-len(line)) * ' ')
    if line[0:6].rstrip() not in ACCEPTED:
      #~ print >> sys.stderr, '%-6sk' %line[0:6].rstrip()
      continue
# make clean divider lines without additional stuff that might hurt later
    elif line[0:6].rstrip() in DIVIDERS:
      if line[0:5] == 'MODEL':
        line = line.rstrip()
      else:
        line = line[0:6].rstrip()
      line = (line + (80-len(line)) * ' ')
    else:
      pass
# split the line into records
    record_type.append('%-6s' % line[0:6])
    atom_number.append(line[6:11])
    blank1.append(line[11:12])
    atom_name.append(line[12:16])
    alt_loc_ind.append(line[16:17]); residue_name.append(line[17:20])
    blank2.append(line[20:21])
    chain_id.append(line[21:22])
    residue_number.append(line[22:26])
    insertion_code.append(line[26:27])
    blank3.append(line[27:30])
    x.append(line[30:38]); y.append(line[38:46]); z.append(line[46:54])
    occupancy.append(line[54:60]); bfactor.append(line[60:66])
    blank4.append(line[66:76]); element.append(line[76:78])
    charge.append(line[78:80])
  insert_resnums = []; insert_resnames = []; chains = []
  recordlist = []
  for i, record in enumerate(record_type):
# determine insertion code
    if insertion_code[i] != ' ' and residue_number[i]+insertion_code[i] not in insert_resnums:
      insert_resnums.append(residue_number[i]+insertion_code[i])
      insert_resnames.append(residue_name[i])
# determine chain id and record it, if not yet found before
    if chain_id[i] != ' ' and chain_id[i] not in chains:
      chains.append(chain_id[i])
    record = [record_type[i], atom_number[i], blank1[i], atom_name[i], alt_loc_ind[i],
                  residue_name[i], blank2[i], chain_id[i], residue_number[i],
                  insertion_code[i], blank3[i], x[i], y[i], z[i],
                  occupancy[i], bfactor[i], blank4[i], element[i], charge[i]]
# append the accepted record to the overlall record list
    recordlist.append(record)
# report findings so far to the screen
  if chains:
    print >> sys.stderr, "\n----------Chains"
    print >> sys.stderr, "The following (original) chains have been found:"
    for chain in chains:
      print >> sys.stderr, chain
  if insert_resnums:
    print >> sys.stderr, "\n----------Insertions"
    print >> sys.stderr, "The following (original-number) residues were considered as insertions:"
    print >> sys.stderr, "They will be renumbered 'normally' in the final 1-N sequence." 
    for i in range(0,len(insert_resnums)):
      print >> sys.stderr, "%s%s"%(insert_resnames[i],(insert_resnums[i]))

  return recordlist

#==================================================
def pdb_write(recordlist, filename, cnct=''):
#==================================================
# uses a record list as created in pdb_read and writes it out to the filename
# using the format below
  if filename == 'stdout':
    pdbout = sys.stdout
  else:  
    pdbout  = open(filename, 'w')
  format = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
  for i, record in enumerate(recordlist):
    pdbout.write(format % tuple(record))
  pdbout.write(cnct)
  pdbout.write('END'+(77 * ' ')+'\n')  
  if pdbout != sys.stdout:
    pdbout.close()

#==================================================
def prot_only(recordlist):
#==================================================
# this strips any residues not recognized by Amber libraries...
# in a personalized Amber installation with additional libraries, 
# you might consider extending this list
  global RESPROT
  protlist=[]
  for record in recordlist:
    if record[5] not in RESPROT:
      continue
    else:
      protlist.append(record)
  return protlist

#==================================================
def remove_hydrogens(recordlist):
#==================================================
  nohlist = []

  for record in recordlist:
    if record[3][0] == 'H' or record[3][1] == 'H':
      continue
    else:
      nohlist.append(record)

# return the record list with all hydrogens removed
  return nohlist

#========================================
def remove_water(recordlist, filename):
#========================================
# removes all water molecules of option -d was specified
  drylist = []; waterlist = []; nwaters = 0
  is_prev_water = 0
# format = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
# watpdb = open(filename+'_water.pdb', 'w')
  for record in recordlist:
# if previous record was water, remove TER card
    if is_prev_water and record[0]=='TER   ':
      is_prev_water=0      
# if oxygen, then count water
    elif (record[5] == 'HOH' or record[5] == 'WAT') and 'O' in record[3]:
      nwaters +=1
      waterlist.append(record)
      is_prev_water=1
#     watpdb.write(format % tuple(record))
      continue
# if not oxygen, just remove, but do not count, since this is probably hydrogen
    elif  record[5] == 'HOH' or record[5] == 'WAT':
      is_prev_water=1
      continue
    else:
      drylist.append(record)
      is_prev_water=0
 
# report the water removal to the screen
  print >> sys.stderr, "\n---------- Water"
  print >> sys.stderr, "%d water molecules have been removed"%nwaters
  print >> sys.stderr, "and stored in the file %s_water.pdb"%filename
# return the dry record list with all water removed
  pdb_write(waterlist, filename+'_water.pdb')
  return drylist

#========================================
def remove_mostpop_altloc(recordlist,filename):
#========================================
  noaltlist = []
  minor_altloc = open(filename+'_minor_altloc.pdb', 'w')
  pdbformat = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
  report_list=[]
          
#keep most populous conformation
  import collections
  #n_altlocs is a dict where key is a unique atom identifier and value is
  # a two element list. First element contains list of occupancy values 
  # of all occurrences of that atom. Second element counts the appearance
  # of the atom during the second iteration.
  n_altlocs=collections.OrderedDict()
#first iteration to set up dict with occupancies for each altloc atom  
  for record in recordlist:
    id = "%s_%s" %(record[8],record[3])
    if record[4] != ' ':
      if id not in n_altlocs.keys():
        n_altlocs[id] = [ [float(record[14])], 0 ]
      else:
        n_altlocs[id][0].append( float(record[14]) )
        
#now iterate again.      
  for record in recordlist:
    id = "%s_%s" %(record[8],record[3])
    if id in n_altlocs.keys():
      # find index of highest occupancy occurence
      alt_max = n_altlocs[id][0].index( max(n_altlocs[id][0]) )
      # if current atom has highest occupancy altloc, add to list
      if n_altlocs[id][1] == alt_max:
        record[4] = ' '
        noaltlist.append(record)
        report_list.append(record)
      # otherwise write to minor_altloc file and discard  
      else:
        minor_altloc.write(pdbformat % tuple(record))
      n_altlocs[id][1] += 1  
    #if current atom has no altlocs, add directly to list    
    else:
      record[4] = ' '
      noaltlist.append(record)
      
  if n_altlocs:
    print >> sys.stderr, "\n---------- Alternate Locations (Original Residues!)"
    print >> sys.stderr, "The following atoms had alternate locations:"
    for record in report_list:
      print >> sys.stderr, "%s_%-4s %s"%(record[5], record[8].strip(), record[3])
    print >> sys.stderr, "The alternate coordinates have been discarded."
    print >> sys.stderr, "Only the highest occupancy of each atom was kept."         
    print >> sys.stderr, "Alternate conformations were printed to %s_minor_altloc.pdb" %filename
    
  minor_altloc.close()    
  return noaltlist

#========================================
def remove_altloc(recordlist):
#========================================
  noaltlist = []
  altloc_resnum = []; altloc_resname = []
  for record in recordlist:
# we accept only altlocs 'A' and '1'
    if record[4] != ' ' and record[4] != 'A' and record[4] != '1':
      if record[8] not in altloc_resnum:
        altloc_resnum.append(record[8])
        altloc_resname.append(record[5])
      continue

    else:
      record[4] = ' '
      noaltlist.append(record)
  if altloc_resnum:
    print >> sys.stderr, "\n---------- Alternate Locations (Original Residues!)"
    print >> sys.stderr, "The following residues had alternate locations:"

    for i, rname in enumerate(altloc_resname):
      print >> sys.stderr, "%s_%d"%(rname, int(altloc_resnum[i]))

    print >> sys.stderr, "The alternate coordinates have been discarded."
    print >> sys.stderr, "Only the first occurrence for each atom was kept."

  return noaltlist

#==================================================
def atom_wrap(recordlist):
#==================================================
# !!! this function should always be called !!!

# wraps 4-letter hydrogens
  wraplist = []

  for record in recordlist:
    if record[0] != 'ATOM  ' and record[0] != 'HETATM':
      wraplist.append(record)
      continue

# shifts 3-letter atoms if needed
    elif record[3][0] != ' ' and record[3][3] == ' ':
      atomname = record[3][3] + record[3][0:3]
      record[3] = atomname
      wraplist.append(record)
      continue

    else:
      wraplist.append(record)
      continue

  return wraplist

#========================================
def renumber(recordlist, filename):
#========================================
  table = open('%s_renum.txt'%filename, 'w')
  renumbered = []; current = -100; iatom = 1
  original = []; oriresname = []; final = []; finresname = []

  for record in recordlist:
    if not 'ATOM' in record[0] and not 'HETATM' in record[0]:
      renumbered.append(record)

    elif current == -100:
      actualnum = record[8]+record[9]
      actualname = record[5]
      if record[3] == ' CA ' or record[3] == 'CH3':
        original.append(record[8]+record[9])
        oriresname.append(record[5])

      record[8] = 1
      record[9] = ' '
      current = 1
      record[1] = iatom
      iatom += 1

      if record[3] == ' CA ' or record[3] == ' CH3':
        final.append(record[8])
        finresname.append(record[5])

      renumbered.append(record)

    elif record[8]+record[9] == actualnum and record[5] == actualname:

      if record[3] == ' CA ' or record[3] == ' CH3':
        original.append(record[8]+record[9])
        oriresname.append(record[5])

      record[8] = current
      record[9] = ' '
      record[1] = iatom
      iatom += 1

      if record[3] == ' CA ' or record[3] == ' CH3':
        final.append(record[8])
        finresname.append(record[5])
      renumbered.append(record)

    elif record[8]+record[9] != actualnum or record[5] != actualname :
      actualnum = record[8]+record[9]
      actualname = record[5]
      if record[3] == ' CA ' or record[3] == ' CH3':
        original.append(record[8]+record[9])
        oriresname.append(record[5])

      current += 1
      record[8] = current
      record[9] = ' '
      record[1] = iatom
      iatom += 1

      if record[3] == ' CA ' or record[3] == ' CH3':
        final.append(record[8])
        finresname.append(record[5])
      renumbered.append(record)


  for i in range(0, len(original)):
    table.write("%3s %5s    %3s %5s\n" %(oriresname[i], (original[i]),
                                    finresname[i], final[i]) )

  return renumbered

#========================================
def non_standard(recordlist, filename):
#========================================
# define the common AA and less common AA names that make up proteins
# and that are recognized by Amber routines in ATOM (or HETATM) records
  RES = ('A', 'A3', 'A5', 'ACE',
         'ALA', 'AN', 'ARG', 'ASH',
         'ASN', 'ASP', 'Br-', 'C',
         'C3', 'C5', 'CN', 'CYM',
         'CYS', 'CYX', 'Cl-', 'Cs+',
         'DA', 'DA3', 'DA5', 'DAN',
         'DC', 'DC3', 'DC4', 'DC5',
         'DCN', 'DG', 'DG3', 'DG5',
         'DGN', 'DT', 'DT3', 'DT5',
         'DTN', 'F-', 'G', 'G3',
         'G5', 'GLH', 'GLN', 'GLU',
         'GLY', 'GN', 'HID', 'HIE',
         'HIP', 'HIS', 'HOH', 'HYP',
         'I-', 'ILE', 'K+', 'LEU',
         'LYN', 'LYS', 'Li+', 'MET',
         'Mg+', 'NHE', 'NME', 'Na+',
         'OHE', 'PHE', 'PL3', 'PRO',
         'Rb+', 'SER', 'SPC', 'SPF',
         'SPG', 'T4E', 'THR', 'TP3',
         'TP4', 'TP5', 'TPF', 'TRP',
         'TYR', 'U', 'U3', 'U5',
         'UN', 'VAL', 'WAT', 'U5',
         'UN', 'VAL', 'WAT')

  hetero = open(filename+'_nonprot.pdb', 'w')
  pdbformat = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
  ns_resname = []

  for record in recordlist:
    if record[5].strip() not in RES and record[5] != '   ':
      hetero.write(pdbformat % tuple(record))
      if record[5] not in ns_resname:
          ns_resname.append(record[5])

  if ns_resname:
    print >> sys.stderr, "\n---------- Non-Standard Residues"
    print >> sys.stderr, "The following non-standard residue names in the original PDB file"
    print >> sys.stderr, "are not recognized by Amber and have been written to the separate"
    print >> sys.stderr, "file %s_nonprot.pdb"%filename
    print >> sys.stderr, "\n".join(ns_resname)

  return ns_resname

#========================================
def non_standard_elbow(recordlist):
#========================================
# define the common AA and less common AA names that make up proteins
# and that are recognized by Amber routines in ATOM (or HETATM) records
  RES = ('A', 'A3', 'A5', 'ACE',
         'ALA', 'AN', 'ARG', 'ASH',
        'ASN', 'ASP', 'BA', 'BR',
        'C', 'C3', 'C5', 'CA',
        'CD', 'CL', 'CN', 'CO',
        'CS', 'CU', 'CYM', 'CYS',
        'CYX', 'DA', 'DA3', 'DA5',
        'DAN', 'DC', 'DC3', 'DC4',
        'DC5', 'DCN', 'DG', 'DG3',
        'DG5', 'DGN', 'DT', 'DT3',
        'DT5', 'DTN', 'EU', 'F',
        'FE2', 'G', 'G3', 'G5',
        'GLH', 'GLN', 'GLU', 'GLY',
        'GN', 'HG', 'HID', 'HIE',
        'HIP', 'HIS', 'HOH', 'HYP',
        'ILE', 'IOD', 'K', 'LEU',
        'LI', 'LYN', 'LYS', 'MET',
        'MG', 'MN', 'NA', 'NHE',
        'NI', 'NME', 'OHE', 'PB',
        'PD', 'PHE', 'PL3', 'PRO',
        'PT', 'RB', 'SER', 'SPC',
        'SPF', 'SPG', 'SR', 'T4E',
        'THR', 'TP3', 'TP4', 'TP5',
        'TPF', 'TRP', 'TYR', 'U',
        'U3', 'U5', 'UN', 'V2+',
        'VAL', 'WAT', 'YB2', 'ZN')

  
  pdbformat = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
  ns_resname = []

  for record in recordlist:
    if record[5].strip() not in RES and record[5] != '   ':
      # if 1st instance of this residue, add atom and get chain/resid
      if record[5].strip() not in ns_resname:
        ns_resname.append(record[5].strip())
        try: f.close()
        except: pass  
        f=open('4antechamber_%s.pdb' %(record[5].strip()), 'w')
        resid=record[8]
        chain=record[7]
        f.write(pdbformat % tuple(record))
      else: 
        #if next atom in the 1st instance of the residue, add atom
        if record[8]==resid and record[7]==chain:
          f.write(pdbformat % tuple(record))
        else:
          #if next instance of the residue, close f and continue
          if not f.closed: f.close()  

  return ns_resname

#========================================
def find_his(recordlist):
#========================================

  amber_hist = {}
  standard_hist = {}
  res_to_change = {}
  for record in recordlist:
    if record[5] in ('HID', 'HIP', 'HIE') and record[3] == ' CA ':
      amber_hist[record[8]] = record[5]
    if record[5] == 'HIS':
        if record[8] not in standard_hist:
          standard_hist[record[8]] = [record[3].strip()]
        else:
          standard_hist[record[8]].append(record[3].strip())

  print >> sys.stderr, "\n---------- Histidines (Renumbered Residues!)"
  if len(amber_hist) == 0 and len(standard_hist) == 0:
    print >> sys.stderr, "No histidine residues found."
    return recordlist

  if len(amber_hist) > 0:
    print >> sys.stderr, "The following histidine residues are already named according to Amber convention."
    for rnum, rname in amber_hist.iteritems():
      print >> sys.stderr, '%s_%d' % (rname, rnum)
  
  if len(standard_hist) > 0:
    for rnum,atoms in standard_hist.items():
      if 'HD1' in atoms and 'HE2' in atoms:
        res_to_change[rnum] = 'HIP'
        del standard_hist[rnum]
      elif 'HD1' in atoms and 'HE2' not in atoms:
        res_to_change[rnum] = 'HID'
        del standard_hist[rnum]
      elif 'HD1' not in atoms and 'HE2' in atoms:
        res_to_change[rnum] = 'HIE'        
        del standard_hist[rnum]

  if len(res_to_change) > 0:
    print >> sys.stderr, "The following HIS residues will be changed to Amber convention names:"
    for rnum, rname in res_to_change.items():
       print >> sys.stderr, "HIS %d --> %s." %(rnum, rname)
    for record in recordlist:
      if record[8] in res_to_change.keys():
        record[5] = res_to_change[record[8]]    

  if len(standard_hist) > 0:
    print >> sys.stderr, "It was not possible to determine the protonation state of the following HIS"
    print >> sys.stderr, "residues based on presence of hydrogens. Amber will consider them as HIE"
    print >> sys.stderr, "(epsilon-HIS) by default. If other protonation state desired change to HID"
    print >> sys.stderr, "(delta-HIS) or HIP (protonated HIS) by hand."
    for rnum in standard_hist.keys():
       print >> sys.stderr, "HIS_%d" %rnum
       
  return recordlist



#========================================
def constph(recordlist):
#========================================
  print >> sys.stderr, "\n---------- Constant pH Simulation"
  as4, gl4, hip = [],[],[]
  for record in recordlist:
    if record[5] == 'ASP':
      record[5] = 'AS4'
      as4.append(record[8])
    elif record[5] == 'GLU':
      record[5] = 'GL4'
      gl4.append(record[8])
    elif record[5] == 'HIS':
      record[5] = 'HIP'
      hip.append(record[8])
    else:
      continue

  print >> sys.stderr, "ASP --> AS4: %d" %len(set(as4))
  print >> sys.stderr, "GLU --> GL4: %d" %len(set(gl4))
  print >> sys.stderr, "HIS --> HIP: %d" %len(set(hip))

  return recordlist

#========================================
def find_disulfide(recordlist, filename):
#========================================
  cys_residues = [];  cys_sgx = []; cys_sgy = []; cys_sgz = []
  cyx_residues = [];  cys_sqn = []; ncys = 0; ncyx = 0

  print >> sys.stderr, "\n---------- Cysteines in Disulfide Bonds (Renumbered Residues!)"
  for record in recordlist:

    if 'SG' in record[3] and ('CYS' in record[5] or 'CYX' in record[5]):
      cys_residues.append(record[8])
      cys_sgx.append(record[11])
      cys_sgy.append(record[12])
      cys_sgz.append(record[13])
      cys_sqn.append(record[1])
      ncys += 1
  
  cnct=''
  if ncys > 0:
      
    sslink = open('%s_sslink'%filename, 'w')
    dist = [[0 for i in range(ncys)] for j in range(ncys)]
    for i in range(0, ncys-1):
      for j in range(i+1, ncys):
        dx = float(cys_sgx[i]) - float(cys_sgx[j])
        dx2 = dx*dx
        dy = float(cys_sgy[i]) - float(cys_sgy[j])
        dy2 = dy*dy
        dz = float(cys_sgz[i]) - float(cys_sgz[j])
        dz2 = dz*dz
        dist[i][j] = sqrt(dx2 +dy2 +dz2)
        if 2.5 > dist[i][j] > 0.1:
          cyx_residues.append(cys_residues[i])
          cyx_residues.append(cys_residues[j])
          print >> sys.stderr,("CYS_%s - CYS_%s: S-S distance = %f Ang."%(cys_residues[i],
                 cys_residues[j], dist[i][j]))
          sslink.write('%s %s\n'%(cys_residues[i], cys_residues[j]))
          ncyx += 1
          cnct += 'CONECT%5d%5d\n' %(cys_sqn[i],cys_sqn[j])

# rename the CYS to CYX for disulfide-involved cysteines
    for record in recordlist:
      if record[8] in cyx_residues:
        record[5] = 'CYX'
      else:
        continue
  if ncyx:
    print >> sys.stderr, "The above CYS have been renamed to CYX in the new PDB file."
    print >> sys.stderr, "Disulfide bond CONECT cards have been added to the new PDB file."

  else:
    print >> sys.stderr, "No disulfide bonds have been detected."
  return recordlist, cnct

#========================================
def find_gaps(recordlist):
#========================================
  global RESPROT
  ca_atoms = []
  gaplist = []
 
  def is_ter(index):
    resnum = recordlist[index][8]
    next = 1
    while True:
      if recordlist[index+next][0] in ['TER   ', 'MODEL ', 'ENDMDL']:
        return True
      elif recordlist[index+next][8]==resnum:
        next+=1
      else:
        return False
          
  for i,record in enumerate(recordlist):
    if ('CA' in record[3] or 'CH3' in record[3]) and record[5] in RESPROT:
      ca_atoms.append(i)

  nca = len(ca_atoms)
  ngaps = 0

  for i in range(nca-1):
    if is_ter(ca_atoms[i]):
      continue
    ca1 = recordlist[ca_atoms[i]]
    ca2 = recordlist[ca_atoms[i+1]]
    dx = float(ca1[11]) - float(ca2[11])
    dy = float(ca1[12]) - float(ca2[12])
    dz = float(ca2[13]) - float(ca2[13])
    gap = sqrt(dx*dx +dy*dy +dz*dz)
    
    if gap > 5.0:
      gaprecord = (gap, ca1[5], int(ca1[8]), ca2[5], int(ca2[8]))
      gaplist.append(gaprecord)
      ngaps += 1

  if ngaps > 0:
    print >> sys.stderr, "\n---------- Gaps (Renumbered Residues!)"
    cformat = "gap of %lf A between %s_%d and %s_%d"

    for i, gaprecord in enumerate(gaplist):
      print >> sys.stderr, (cformat % tuple(gaprecord))

    print >> sys.stderr, "You MUST (!!!) insert a TER record between the residues listed above and"
    print >> sys.stderr, "consider to introduce caps (ACE and NME) at the dangling N- and C-terminals."

  return()

#========================================
def find_incomplete(recordlist):
#========================================
# finds residues with missing heavy atoms in the following list of residues;
# dictionary with number of heavy atoms:
  #PAJ TODO:complete this list with nucleic acids
  HEAVY = {'ALA':5,  'ARG':11, 'ASN':8,  'ASP':8,
           'CYS':6,  'GLN':9,  'GLU':9,  'GLY':4,
           'HIS':10, 'ILE':8,  'LEU':8,  'LYS':9,
           'MET':8,  'PHE':11, 'PRO':7,  'SER':6,
           'THR':7,  'TRP':14, 'TYR':12, 'VAL':7,
           'HID':10, 'HIE':10, 'HIN':10, 'HIP':10,
           'CYX':6,  'ASH':8,  'GLH':9,  'LYH':9}
  print >> sys.stderr, '\n---------- Missing Heavy Atoms (Renumbered Residues!)' 
  resnum = []
  resname = []
  resheavy = []
  flag = 0
  nheavy = 1
  length = len(recordlist)
  for i, record in enumerate(recordlist):
    if i == length-1:
      for j, res in enumerate(resnum):
        missing = HEAVY[resname[j]] - resheavy[j]
        if missing > 0:
          flag = 1
          print >> sys.stderr, "%s_%s misses %d heavy atom(s)"%(resname[j], res, missing) 
      if flag == 0:
        print >> sys.stderr, "None"
      return()
    if not HEAVY.has_key(record[5]):
      continue
    if recordlist[i+1][8] == record[8]:
      nheavy += 1
    else:
      resnum.append(record[8])
      resname.append(record[5])
      resheavy.append(nheavy)
      nheavy = 1
      continue
  return()

def run(arg_pdbout, arg_pdbin, 
        arg_nohyd = True, 
        arg_dry   = False, 
        arg_prot  = False, 
        arg_noter = False,
        arg_constph = False,
        arg_mostpop = False,
        arg_reduce = False,
        arg_model = 0,
        arg_elbow = False
        ):
  filename, extension = os.path.splitext(arg_pdbout)
  pdbin = arg_pdbin
  
  # optionally run reduce on input file
  if arg_reduce:
    if arg_pdbin == 'stdin':
      pdbfile = sys.stdin
    else:
      pdbfile = open(arg_pdbin, 'r')
    try:
      process = subprocess.Popen(['reduce', '-BUILD', '-NUC', '-'], stdin=pdbfile,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = process.communicate()
      if process.wait():
        print >> sys.stderr, ("REDUCE returned non-zero exit status: "
                              "See reduce_info.log for more details")
        open('reduce_info.log', 'w').write(err)
      # Uncomment the following to print out the reduce log even if it worked
#     else:
#       open('reduce_info.log', 'w').write(err)
      pdbh = StringIO(out)
      recordlist = pdb_read(pdbh, arg_noter, arg_model)
    finally:
      if pdbfile is not sys.stdin: pdbfile.close()
  else:  
    recordlist = pdb_read(pdbin, arg_noter, arg_model)
  
  # wrap all atom names to pure standard (always):======================
  recordlist = atom_wrap(recordlist)
  # remove alternate locations and keep only the first one:=============
  if arg_mostpop:
    recordlist = remove_mostpop_altloc(recordlist, filename)
  else:
    recordlist = remove_altloc(recordlist)
  # remove hydrogens if option -y is used:==============================
  if arg_nohyd:
    recordlist = remove_hydrogens(recordlist)
  # find non-standard Amber residues:===================================
  non_standard(recordlist, filename)
  ns_names = None
  if arg_elbow:
    ns_names=non_standard_elbow(recordlist)
  # keep only protein:==================================================
  if arg_prot:
    recordlist = prot_only(recordlist)
  # remove water if -d option used:=====================================
  if arg_dry:
    recordlist = remove_water(recordlist, filename)
  # renumber atoms and residues:========================================
  recordlist = renumber(recordlist, filename)
  # after this call, residue numbers refer to the ***new*** PDB file
  #=====================================================================
  # find histidines that might have to be changed
  if arg_constph:
    recordlist = constph(recordlist)
  else:
    recordlist = find_his(recordlist)
  # find possible S-S in the final protein:=============================
  recordlist, cnct = find_disulfide(recordlist, filename)
  # find possible gaps:==================================================
  find_gaps(recordlist)
  # count heavy atoms
  find_incomplete(recordlist)
  # =====================================================================
  # make final output to new PDB file
  pdb_write(recordlist, arg_pdbout, cnct)
  print >> sys.stderr, ""
  return ns_names
    
#========================================main===========================
if __name__ ==  "__main__":
  parser = OptionParser(version=__version__)
  parser.add_option("-i","--in", metavar = "FILE", dest = "pdbin",
                    help = "PDB input file                      (default: stdin)",
                    default='stdin')
  parser.add_option("-o","--out", metavar = "FILE", dest = "pdbout",
                    help = "PDB output file                     (default: stdout)",
                    default='stdout')
  parser.add_option("-y","--nohyd", action = "store_true", dest = "nohyd",
                    help = "remove all hydrogen atoms           (default: no)")
  parser.add_option("-d","--dry", action = "store_true", dest = "dry",
                    help = "remove all water molecules          (default: no)")
  parser.add_option("-p", "--prot", action = "store_true", dest = "prot",
                    help = "keep only Amber-compatible residues (default: no)")
  parser.add_option("--noter", action =  "store_true", dest = "noter",
                    help = "remove TER, MODEL, ENDMDL cards     (default: no)")
  parser.add_option("--constantph", action = "store_true", dest = "constantph",
                    help = "rename GLU,ASP,HIS for constant pH simulation")
  parser.add_option("--most-populous", action = "store_true", dest = "mostpop",
                    help = "keep most populous alt. conf. (default is to keep 'A')")
  parser.add_option("--reduce", action = "store_true", dest = "reduce",
                    help = "Run Reduce first to add hydrogens.  (default: no)")     
  parser.add_option("--model", type = "int", dest = "model", default = 0,
                    help = "Model to use from a multi-model pdb file (integer).  (default: use all models)")                                        
  (opt, args) = parser.parse_args()

  if opt.pdbin == opt.pdbout:
    print >> sys.stderr, "The input and output file names cannot be the same!\n"
    sys.exit(1)

  # Make sure that if we are reading from stdin it's being directed from a pipe
  # or a file. We don't want to wait for user input that will never come.

  if opt.pdbin == 'stdin':
    if os.isatty(sys.stdin.fileno()):
      sys.exit(parser.print_help() or 1)

  run(opt.pdbout, opt.pdbin, opt.nohyd, opt.dry, opt.prot, opt.noter, 
      opt.constantph, opt.mostpop, opt.reduce, opt.model)
