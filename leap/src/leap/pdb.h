/*
 *      Copyright (c) 1989 The Regents of the University of California.
 *      All rights reserved.
 *
 *      Redistribution and use in source and binary forms are permitted
 *      provided that the above copyright notice and this paragraph are
 *      duplicated in all such forms and that any documentation,
 *      advertising materials, and other materials related to such
 *      distribution and use acknowledge that the software was developed
 *      by the University of California, San Francisco.  The name of the
 *      University may not be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *      THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 *      IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 *      WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *      $Id: pdb.h,v 10.0 2008/04/15 23:22:21 case Exp $
 *
 *      Based on Brookhaven National Laboratory Protein Data Bank, March 1989
 *
 *      C structure declarations
 */

/*      Modifications induced by the implementation of the savemol2 command
 *      Christine Cezard (2007) 
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 *      added # define        MOL2_ATOM       36
 *      the value of PDB_NUM_R  is now  37
 *      added lines 86 and 120 
*/
 
# ifndef PDB_H
# define        PDB_H

# include       <stdio.h>

# define        PDB_RECLEN      80              /* PDB record length */
# define        PDB_BUFSIZ      PDB_RECLEN + 2  /* + '\n' + '\0' */             

# define        PDB_UNKNOWN     0

/* records originally in alphabetical order */

# define        PDB_ANISOU      1
# define        PDB_ATOM        2
# define        PDB_AUTHOR      3
# define        PDB_COMPND      4
# define        PDB_CONECT      5
# define        PDB_CRYST1      6
# define        PDB_END         7
# define        PDB_FORMUL      8
# define        PDB_FTNOTE      9
# define        PDB_HEADER      10
# define        PDB_HELIX       11
# define        PDB_HET         12
# define        PDB_HETATM      13
# define        PDB_JRNL        14
# define        PDB_MASTER      15
# define        PDB_MTRIX       16
# define        PDB_OBSLTE      17
# define        PDB_ORIGX       18
# define        PDB_REMARK      19
# define        PDB_REVDAT      20
# define        PDB_SCALE       21
# define        PDB_SEQRES      22
# define        PDB_SHEET       23
# define        PDB_SIGATM      24
# define        PDB_SIGUIJ      25
# define        PDB_SITE        26
# define        PDB_SOURCE      27
# define        PDB_SPRSDE      28
# define        PDB_SSBOND      29
# define        PDB_TER         30
# define        PDB_TURN        31
# define        PDB_TVECT       32
# define        PDB_USER        33
# define        PDB_MODEL       34
# define        PDB_ENDMDL      35
# define        MOL2_ATOM       36

# define        PDB_NUM_R       37


typedef char    mol2_atom[3];
typedef char    pdb_date[10];
typedef char    pdb_aname[5];           /* atom name - NO2* */
typedef char    pdb_rname[4];           /* residue name - ALA */
typedef char    pdb_pname[5];           /* pdb name - 9lyz */
typedef char    pdb_id[4];              /* generic short id field */
typedef double  pdb_float;              /* size of floating point */

typedef struct  {                       /* residue info */
        pdb_rname       name;
        char            chain_id;
        int             seq_num;
        char            insert_code;
}       pdb_residue;

/*
 *      structures declarations for each record type
 */

struct  pdb_unknown     {
        char    junk[81];
};
struct  pdb_anisou      {
        int             serial_num;
        pdb_aname       name;
        char            alt_loc;
        pdb_residue     residue;
        int             u[6];
};
struct  pdb_atom                {
        int             serial_num;
        pdb_aname       name;
        char            alt_loc;
        pdb_residue     residue;
	mol2_atom       type_at;
        pdb_float       x, y, z;
        pdb_float       occupancy, temp_factor;
        int             ftnote_num;
};
struct  pdb_author      {
        char    data[61];
        char    continuation;
};
# define        pdb_compnd      pdb_author
struct  pdb_conect      {
        int     serial_num;
        int     covalent[4];
        struct  {
                int     hydrogen[2];
                int     salt;
        }       bonds[2];
};
struct  pdb_cryst1      {
        pdb_float       a, b, c;
        pdb_float       alpha, beta, gamma;
        char            space_grp[12];
        int             z;
};
/* no structure for PDB_END */
struct  pdb_formul      {
        int             component;
        pdb_rname       het_id;
        int             continuation;
        char            exclude;        /* * to exclude */
        char            formula[52];
};
struct  pdb_ftnote      {
        int     num;
        char    text[60];
};
struct  pdb_header      {
        char            class[41];
        pdb_date        date;
        pdb_pname       id;
};
struct  pdb_helix               {
        int             serial_num;
        pdb_id          id;
        pdb_residue     residues[2];
        int             class;
        char            comment[31];
};
struct  pdb_het         {
        pdb_residue     het_grp;
        int             num_atoms;
        char            text[41];
};
# define        pdb_hetatm      pdb_atom
# define        pdb_jrnl        pdb_author
struct  pdb_master      {
        int     num_remark;
        int     num_ftnote;
        int     num_het;
        int     num_helix;
        int     num_sheet;
        int     num_turn;
        int     num_site;
        int     num_transform;
        int     num_coordinate;
        int     num_ter;
        int     num_conect;
        int     num_seqres;
};
struct  pdb_mtrix               {
        int             row_num;
        int             serial_num;
        pdb_float       m1, m2, m3, v;
        int             given;
};
struct  pdb_obslte      {
        int             continuation;
        pdb_date        date;
        pdb_pname       old_id;
        pdb_pname       id_map[8];
};
struct  pdb_origx               {
        int             row_num;
        pdb_float       o1, o2, o3, t;
};
# define        pdb_remark      pdb_ftnote
struct  pdb_revdat      {
        int             modification;
        int             continuation;
        pdb_date        date;
        char            id[8];
        char            mod_type;
        char            corrections[31];
};
struct  pdb_scale               {
        int             row_num;
        pdb_float       s1, s2, s3, u;
};
struct  pdb_seqres      {
        int             serial_num;
        char            chain_id;
        int             count;
        pdb_rname       names[13];
};
struct  pdb_sheet               {
        int             strand_num;
        pdb_id          id;
        int             count;
        pdb_residue     residues[2];
        int             sense;
        struct  {
                pdb_aname       name;
                pdb_residue     residue;
        }               atoms[2];
};
# define        pdb_sigatm      pdb_atom
# define        pdb_siguij      pdb_anisou
struct  pdb_site                {
        int             seq_num;
        pdb_id          id;
        int             count;
        pdb_residue     residues[4];
};
# define        pdb_source      pdb_author
struct  pdb_sprsde      {
        int             continuation;
        pdb_date        date;
        pdb_pname       id;
        pdb_pname       supersede[8];
};
struct  pdb_ssbond      {
        int             seq_num;
        pdb_residue     residues[2];
        char            comment[31];
};
struct  pdb_ter         {
        int             serial_num;
        pdb_residue     residue;
};
struct  pdb_turn                {
        int             seq_num;
        pdb_id          id;
        pdb_residue     residues[2];
        char            comment[31];
};
struct  pdb_tvect               {
        int             serial_num;
        pdb_float       t1, t2, t3;
        char            comment[31];
};
struct  pdb_user        {
        char    subtype[3];
        char    text[67];
};
struct  pdb_model       {
        int     num;
};
/* no structure for PDB_ENDMDL */

 



typedef struct  pdb_record      {
        int     record_type;
        union   {
                struct  pdb_unknown     unknown;
                struct  pdb_anisou      anisou;
                struct  pdb_atom        atom;
                struct  pdb_author      author;
                struct  pdb_compnd      compnd;
                struct  pdb_conect      conect;
                struct  pdb_cryst1      cryst1;
                /* no pdb_end structure */
                struct  pdb_formul      formul;
                struct  pdb_ftnote      ftnote;
                struct  pdb_header      header;
                struct  pdb_helix       helix;
                struct  pdb_het         het;
                struct  pdb_hetatm      hetatm;
                struct  pdb_jrnl        jrnl;
                struct  pdb_master      master;
                struct  pdb_mtrix       mtrix;
                struct  pdb_obslte      obslte;
                struct  pdb_origx       origx;
                struct  pdb_remark      remark;
                struct  pdb_revdat      revdat;
                struct  pdb_scale       scale;
                struct  pdb_seqres      seqres;
                struct  pdb_sheet       sheet;
                struct  pdb_sigatm      sigatm;
                struct  pdb_siguij      siguij;
                struct  pdb_site        site;
                struct  pdb_source      source;
                struct  pdb_sprsde      sprsde;
                struct  pdb_ssbond      ssbond;
                struct  pdb_ter         ter;
                struct  pdb_turn        turn;
                struct  pdb_tvect       tvect;
                struct  pdb_user        user;
                struct  pdb_model       model;
				/* no pdb_endmdl structure */
        }       pdb;
}       pdb_record;


#ifndef BASICS_H
# include       "basics.h"
#endif

extern pdb_record      pdb_read_record(FILE *);
extern void            pdb_write_record(FILE *, pdb_record *, char *, int);
extern void            mol2_write_record(FILE *, pdb_record *, char *, int);

# endif /* PDB_H */
