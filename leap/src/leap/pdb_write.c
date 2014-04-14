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
 *      THIS SOFTWARE IS PROVIDED `AS IS' AND WITHOUT ANY EXPRESS OR
 *      IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 *      WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *      $Id: pdb_write.c,v 10.1 2009/12/28 16:58:42 case Exp $
 *
 *      subroutine for writing PDB format files
 *
 */

/*      Modifications induced by the implementation of the savemol2 command
 *      Christine Cezard (2007) 
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 *      declaration of a case MOL2_ATOM line 372
 */    
  
 
/* LINTLIBRARY */

# include       <stdio.h>
# include       <ctype.h>
# include       "pdb_int.h"

/*
 *      for each pdb record type there is a format reading in the
 *      record values and for printing them out.
 *
 *      The actual format of a line written, is the print format
 *      followed by blank padding to 72 characters, followed by
 *      8 characters of file and line information.
 */
/*
        char            *name;          if NULL, don't print 
        int             line_num;       if name == NULL, don't print 
*/

void
pdb_write_record(FILE *f, pdb_record *r, char *name, int line_num)
{
        char            buffer[PDB_BUFSIZ];
        char            *fmt;
        register struct pdb_sheet       *sh;
        register pdb_residue    *shr0, *shr1, *sha0, *sha1;

        /* convert C structure to pdb record */

        fmt = pdb_record_format[r->record_type].print_format;
        switch (r->record_type) {

        
		
		case PDB_UNKNOWN:
                pdb_sprintf(buffer, fmt, r->pdb.unknown.junk);
                break;

        case PDB_ANISOU:
        case PDB_SIGUIJ:
                pdb_sprintf(buffer, fmt, r->pdb.anisou.serial_num,
                        r->pdb.anisou.name, r->pdb.anisou.alt_loc,
                        r->pdb.anisou.residue.name,
                        r->pdb.anisou.residue.chain_id,
                        r->pdb.anisou.residue.seq_num,
                        r->pdb.anisou.residue.insert_code,
                        r->pdb.anisou.u[0], r->pdb.anisou.u[1],
                        r->pdb.anisou.u[2], r->pdb.anisou.u[3],
                        r->pdb.anisou.u[4], r->pdb.anisou.u[5]);
                break;

        case PDB_ATOM:
        case PDB_HETATM:
        case PDB_SIGATM:
/*
"ATOM  %5d %-4s%C%-3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D"
ATOM      1  N   GLY   204      16.233  16.706  11.826  1.00  0.00
ATOM      6  HA2 GLY   204      15.123  14.988  11.463  1.00  0.00

ATOM  150015 Na+  Na+  1443     -17.714  -6.411   7.833  1.00  0.00

"ATOM  %6d%-4s%C%-3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D";
ATOM  150009Na+  Na+  1439     -23.630   8.766 -18.151  1.00  0.00

ATOM  41150  H1  WAT  8081      21.760   0.196  43.791  1.00  0.00
ATOM  51150  H1  WAT  10581      21.703 -29.969  47.175  1.00  0.00
*/
        {
                int             atom_overflow = 0, residue_overflow = 0;

                if (r->pdb.atom.serial_num > 99999)
                        atom_overflow = 1;
                if (r->pdb.atom.residue.seq_num > 9999)
                        residue_overflow = 1;
                if (atom_overflow  &&  !residue_overflow) {
                        fmt =
                  "ATOM  %6d%-4s%C%-3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D";
                } else if (!atom_overflow  &&  residue_overflow) {
                        fmt =
                  "ATOM  %5d %-4s%C%-3s %C%5d%C  %8.3f%8.3f%8.3f%6.2f%6.2f %3D";
                } else if (atom_overflow  &&  residue_overflow) {
                        fmt =
                  "ATOM  %6d%-4s%C%-3s %C%5d%C  %8.3f%8.3f%8.3f%6.2f%6.2f %3D";
                }
                pdb_sprintf(buffer, fmt, r->pdb.atom.serial_num,
                        r->pdb.atom.name, r->pdb.atom.alt_loc,
                        r->pdb.atom.residue.name,
                        r->pdb.atom.residue.chain_id,
                        r->pdb.atom.residue.seq_num,
                        r->pdb.atom.residue.insert_code,
                        r->pdb.atom.x, r->pdb.atom.y, r->pdb.atom.z,
                        r->pdb.atom.occupancy, r->pdb.atom.temp_factor,
                        r->pdb.atom.ftnote_num);
                break;
        }
        case PDB_AUTHOR:
        case PDB_COMPND:
        case PDB_JRNL:
        case PDB_SOURCE:
                pdb_sprintf(buffer, fmt, r->pdb.author.continuation,
                        r->pdb.author.data);
                break;

        case PDB_CONECT:
                pdb_sprintf(buffer, fmt, r->pdb.conect.serial_num,
                        r->pdb.conect.covalent[0], r->pdb.conect.covalent[1], 
						r->pdb.conect.covalent[2]);
                break;

        case PDB_CRYST1:
                pdb_sprintf(buffer, fmt, r->pdb.cryst1.a, r->pdb.cryst1.b,
                        r->pdb.cryst1.c, r->pdb.cryst1.alpha,
                        r->pdb.cryst1.beta, r->pdb.cryst1.gamma,
                        r->pdb.cryst1.space_grp, r->pdb.cryst1.z);
                break;

        case PDB_END:
                pdb_sprintf(buffer, fmt);
                break;

        case PDB_FORMUL:
                pdb_sprintf(buffer, fmt, r->pdb.formul.component,
                        r->pdb.formul.het_id, r->pdb.formul.continuation,
                        r->pdb.formul.exclude, r->pdb.formul.formula);
                break;

        case PDB_FTNOTE:
                pdb_sprintf(buffer, fmt, r->pdb.ftnote.num, r->pdb.ftnote.text);
                break;

        case PDB_HEADER:
                pdb_sprintf(buffer, fmt, r->pdb.header.class,
                        r->pdb.header.date, r->pdb.header.id);
                break;

        case PDB_HELIX:
                pdb_sprintf(buffer, fmt, r->pdb.helix.serial_num,
                        r->pdb.helix.id,
                        r->pdb.helix.residues[0].name,
                        r->pdb.helix.residues[0].chain_id,
                        r->pdb.helix.residues[0].seq_num,
                        r->pdb.helix.residues[0].insert_code,
                        r->pdb.helix.residues[1].name,
                        r->pdb.helix.residues[1].chain_id,
                        r->pdb.helix.residues[1].seq_num,
                        r->pdb.helix.residues[1].insert_code,
                        r->pdb.helix.class, r->pdb.helix.comment);
                break;

        case PDB_HET:
                pdb_sprintf(buffer, fmt, r->pdb.het.het_grp.name,
                        r->pdb.het.het_grp.chain_id, r->pdb.het.het_grp.seq_num,
                        r->pdb.het.het_grp.insert_code, r->pdb.het.num_atoms,
                        r->pdb.het.text);
                break;

        case PDB_MASTER:
                pdb_sprintf(buffer, fmt, r->pdb.master.num_remark,
                        r->pdb.master.num_ftnote, r->pdb.master.num_het,
                        r->pdb.master.num_helix, r->pdb.master.num_sheet,
                        r->pdb.master.num_turn, r->pdb.master.num_site,
                        r->pdb.master.num_transform,
                        r->pdb.master.num_coordinate, r->pdb.master.num_ter,
                        r->pdb.master.num_conect, r->pdb.master.num_seqres);
                break;

        case PDB_MTRIX:
                pdb_sprintf(buffer, fmt, r->pdb.mtrix.row_num,
                        r->pdb.mtrix.serial_num, r->pdb.mtrix.m1,
                        r->pdb.mtrix.m2, r->pdb.mtrix.m3, r->pdb.mtrix.v,
                        r->pdb.mtrix.given);
                break;

        case PDB_OBSLTE:
                pdb_sprintf(buffer, fmt, r->pdb.obslte.continuation,
                        r->pdb.obslte.date, r->pdb.obslte.old_id,
                        r->pdb.obslte.id_map[0], r->pdb.obslte.id_map[1],
                        r->pdb.obslte.id_map[2], r->pdb.obslte.id_map[3],
                        r->pdb.obslte.id_map[4], r->pdb.obslte.id_map[2],
                        r->pdb.obslte.id_map[6], r->pdb.obslte.id_map[7]);
                break;

        case PDB_ORIGX:
                pdb_sprintf(buffer, fmt, r->pdb.origx.row_num, r->pdb.origx.o1,
                        r->pdb.origx.o2, r->pdb.origx.o3, r->pdb.origx.t);
                break;

        case PDB_REMARK:
                pdb_sprintf(buffer, fmt, r->pdb.remark.num, r->pdb.remark.text);
                break;

        case PDB_REVDAT:
                pdb_sprintf(buffer, fmt, r->pdb.revdat.modification,
                        r->pdb.revdat.continuation, r->pdb.revdat.date,
                        r->pdb.revdat.id, r->pdb.revdat.mod_type,
                        r->pdb.revdat.corrections);
                break;

        case PDB_SCALE:
                pdb_sprintf(buffer, fmt, r->pdb.scale.row_num, r->pdb.scale.s1,
                        r->pdb.scale.s2, r->pdb.scale.s3, r->pdb.scale.u);
                break;

        case PDB_SEQRES:
                pdb_sprintf(buffer, fmt, r->pdb.seqres.serial_num,
                        r->pdb.seqres.chain_id, r->pdb.seqres.count,
                        r->pdb.seqres.names[0], r->pdb.seqres.names[1],
                        r->pdb.seqres.names[2], r->pdb.seqres.names[3],
                        r->pdb.seqres.names[4], r->pdb.seqres.names[5],
                        r->pdb.seqres.names[6], r->pdb.seqres.names[7],
                        r->pdb.seqres.names[8], r->pdb.seqres.names[9],
                        r->pdb.seqres.names[10], r->pdb.seqres.names[11],
                        r->pdb.seqres.names[12]);
                break;

        case PDB_SHEET:
                sh = &r->pdb.sheet;
                shr0 = &sh->residues[0];
                shr1 = &sh->residues[1];
                sha0 = &sh->atoms[0].residue;
                sha1 = &sh->atoms[1].residue;
                pdb_sprintf(buffer, fmt, sh->strand_num,
                        sh->id, sh->count,
                        shr0->name, shr0->chain_id, shr0->seq_num,
                        shr0->insert_code,
                        shr1->name, shr1->chain_id, shr1->seq_num,
                        shr1->insert_code,
                        sh->sense,
                        sh->atoms[0].name,
                        sha0->name, sha0->chain_id, sha0->seq_num,
                        sha0->insert_code,
                        sh->atoms[1].name,
                        sha1->name, sha1->chain_id, sha1->seq_num,
                        sha1->insert_code);
                break;

        case PDB_SITE:
                shr0 = &r->pdb.site.residues[0];
                shr1 = &r->pdb.site.residues[1];
                sha0 = &r->pdb.site.residues[2];
                sha1 = &r->pdb.site.residues[3];
                pdb_sprintf(buffer, fmt, r->pdb.site.seq_num,
                        r->pdb.site.id, r->pdb.site.count,
                        shr0->name, shr0->chain_id, shr0->seq_num,
                        shr0->insert_code,
                        shr1->name, shr1->chain_id, shr1->seq_num,
                        shr1->insert_code,
                        sha0->name, sha0->chain_id, sha0->seq_num,
                        sha0->insert_code,
                        sha1->name, sha1->chain_id, sha1->seq_num,
                        sha1->insert_code);
                break;

        case PDB_SPRSDE:
                pdb_sprintf(buffer, fmt, r->pdb.sprsde.continuation,
                        r->pdb.sprsde.date, r->pdb.sprsde.id,
                        r->pdb.sprsde.supersede[0], r->pdb.sprsde.supersede[1],
                        r->pdb.sprsde.supersede[2], r->pdb.sprsde.supersede[3],
                        r->pdb.sprsde.supersede[4], r->pdb.sprsde.supersede[5],
                        r->pdb.sprsde.supersede[6], r->pdb.sprsde.supersede[7]);
                break;

        case PDB_SSBOND:
                pdb_sprintf(buffer, fmt, r->pdb.ssbond.seq_num,
                        r->pdb.ssbond.residues[0].name,
                        r->pdb.ssbond.residues[0].chain_id,
                        r->pdb.ssbond.residues[0].seq_num,
                        r->pdb.ssbond.residues[0].insert_code,
                        r->pdb.ssbond.residues[1].name,
                        r->pdb.ssbond.residues[1].chain_id,
                        r->pdb.ssbond.residues[1].seq_num,
                        r->pdb.ssbond.residues[1].insert_code,
                        r->pdb.ssbond.comment);
                break;

        case PDB_TER: {
/*
"TER   %5d      %-3s %C%4d%C"
TER    3269      GLN   203
TER   150016      WAT  35296
*/

#ifdef ELABORATE_PDB_TER

                int             atom_overflow = 0, residue_overflow = 0;

                if (r->pdb.atom.serial_num > 99999)
                        atom_overflow = 1;
                if (r->pdb.atom.residue.seq_num > 9999)
                        residue_overflow = 1;
                if (atom_overflow  &&  !residue_overflow) {
                        fmt = "TER   %6d     %-3s %C%4d%C";
                } else if (!atom_overflow  &&  residue_overflow) {
                        fmt = "TER   %5d      %-3a %C%5d%C";
                } else if (atom_overflow  &&  residue_overflow) {
                        fmt = "TER   %6d     %-3s %C%5d%C";
                }
                pdb_sprintf(buffer, fmt, r->pdb.ter.serial_num,
                        r->pdb.ter.residue.name, r->pdb.ter.residue.chain_id,
                        r->pdb.ter.residue.seq_num,
                        r->pdb.ter.residue.insert_code);
#else
                fmt = "TER   ";
                pdb_sprintf(buffer, fmt );
#endif
                break;
        }

        case PDB_TURN:
                pdb_sprintf(buffer, fmt, r->pdb.turn.seq_num,
                        r->pdb.turn.id,
                        r->pdb.turn.residues[0].name,
                        r->pdb.turn.residues[0].chain_id,
                        r->pdb.turn.residues[0].seq_num,
                        r->pdb.turn.residues[0].insert_code,
                        r->pdb.turn.residues[1].name,
                        r->pdb.turn.residues[1].chain_id,
                        r->pdb.turn.residues[1].seq_num,
                        r->pdb.turn.residues[1].insert_code,
                        r->pdb.turn.comment);
                break;

        case PDB_TVECT:
                pdb_sprintf(buffer, fmt, r->pdb.tvect.serial_num,
                        r->pdb.tvect.t1, r->pdb.tvect.t2, r->pdb.tvect.t3,
                        r->pdb.tvect.comment);
                break;

        case PDB_USER:
                pdb_sprintf(buffer, fmt, r->pdb.user.subtype, r->pdb.user.text);
                break;

        case PDB_ENDMDL:
                pdb_sprintf(buffer, fmt);
                break;

        case PDB_MODEL:
                pdb_sprintf(buffer, fmt, r->pdb.model.num);
                break;
				
        case MOL2_ATOM:
                pdb_sprintf(buffer, fmt, r->pdb.atom.serial_num,
                        r->pdb.atom.name, 
                        r->pdb.atom.x, r->pdb.atom.y, r->pdb.atom.z,
                        r->pdb.atom.type_at, r->pdb.atom.residue.seq_num,
						r->pdb.atom.residue.name, r->pdb.atom.temp_factor);
                break;
       			
				
				
        default:
                (void) sprintf(buffer, "unknown pdb record #%d",
                                                                r->record_type);
                break;
        }
        if (name != NULL) {
                if (line_num >= 10000)
                        (void) fprintf(f, "%-72.72s%-4.4s%04d\n", buffer, name,
                                                        line_num % 10000);
                else
                        (void) fprintf(f, "%-72.72s%-4.4s%4d\n", buffer, name,
                                                                line_num);
        } else {
                register        char    *s, *t;
                /* Do not shorten TER/END cards */
                if (r->record_type == PDB_TER || r->record_type == PDB_END)
                        (void) fprintf(f, "%s\n", buffer);
                else {
                        /* find last non-blank in buffer, and shorten it */
                        t = NULL;
                        for (s = buffer; *s != '\0'; s++)
                                if (!isspace(*s))
                                        t = s + 1;
                        if (t == NULL)  /* this should never happen, but ... */
                                t = buffer;
                        *t = '\0';
                        (void) fprintf(f, "%s\n", buffer);
                }
        }
}


# ifdef vms
pdb_write_dummy()
{
        pdb_fmt_dummy();
}
# endif


