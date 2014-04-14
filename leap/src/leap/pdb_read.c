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
 *      subroutine for reading PDB format files
 *
 */

/* LINTLIBRARY */

#include        "basics.h"
#include        "pdb_int.h"

/*
 *      for each pdb record type there is a format reading in the
 *      record values and for printing them out.
 *
 *      The actual format of a line written, is the print format
 *      followed by blank padding to 72 characters, followed by
 *      8 characters of file and line information.
 */

pdb_record
pdb_read_record(FILE *f)
{

        char            buffer[PDB_BUFSIZ];
        pdb_record      r;
        char            *fmt;
        int             c;
        char            *cp;
        register struct pdb_sheet       *sh;
        register pdb_residue    *sha0, *sha1;

                memset( &r, 0, sizeof(r) );

        if (fgets(buffer, PDB_BUFSIZ, f) == NULL) {
                /* at eof or error - default to eof */
                r.record_type = PDB_END;
                return r;
        }

        /* map record name to lower case */
        cp = (char*)strchr(buffer, '\n');
        if (cp != NULL)
                *cp = '\0';
        else
                /* discard extra characters since line too long */
                while ((c = getc(f)) != '\n' && c != EOF)
                                ;

        for (cp = buffer; cp < &buffer[4] && *cp != '\0'; cp++)
                if (isupper(*cp))
                        *cp = cLower(*cp);

        /* convert pdb record to C structure */

        for (c = 0; c < 4; c += 1)
                if (buffer[c] == '\0')
                        break;
        if (c < 4) {
                while (c < 4)
                        buffer[c++] = ' ';
                buffer[4] = '\0';
        }

        r.record_type = PDB_UNKNOWN;
        switch (buffer[0]) {

        case 'a':
                if (STREQN(buffer + 1, "tom", 3))
                        r.record_type = PDB_ATOM;
                else if (STREQN(buffer + 1, "uth", 3))
                        r.record_type = PDB_AUTHOR;
                else if (STREQN(buffer + 1, "nis", 3))
                        r.record_type = PDB_ANISOU;
                break;

        case 'c':
                if (STREQN(buffer + 1, "omp", 3))
                        r.record_type = PDB_COMPND;
                else if (STREQN(buffer + 1, "rys", 3))
                        r.record_type = PDB_CRYST1;
                else if (STREQN(buffer + 1, "one", 3))
                        r.record_type = PDB_CONECT;
                break;

        case 'e':
                if (STREQN(buffer + 1, "nd ", 3))
                        r.record_type = PDB_END;
                else if (STREQN(buffer + 1, "ndm", 3))
                        r.record_type = PDB_ENDMDL;
                break;

        case 'f':
                if (STREQN(buffer + 1, "tno", 3))
                        r.record_type = PDB_FTNOTE;
                else if (STREQN(buffer + 1, "orm", 3))
                        r.record_type = PDB_FORMUL;
                break;

        case 'h':
                if (STREQN(buffer + 1, "eta", 3))
                        r.record_type = PDB_HETATM;
                else if (STREQN(buffer + 1, "ead", 3))
                        r.record_type = PDB_HEADER;
                else if (STREQN(buffer + 1, "et ", 3))
                        r.record_type = PDB_HET;
                else if (STREQN(buffer + 1, "eli", 3))
                        r.record_type = PDB_HELIX;
                break;

        case 'j':
                if (STREQN(buffer + 1, "rnl", 3))
                        r.record_type = PDB_JRNL;
                break;

        case 'm':
                if (STREQN(buffer + 1, "tri", 3))
                        r.record_type = PDB_MTRIX;
                else if (STREQN(buffer + 1, "ast", 3))
                        r.record_type = PDB_MASTER;
                else if (STREQN(buffer + 1, "ode", 3))
                        r.record_type = PDB_MODEL;
                break;

        case 'o':
                if (STREQN(buffer + 1, "bsl", 3))
                        r.record_type = PDB_OBSLTE;
                else if (STREQN(buffer + 1, "rig", 3))
                        r.record_type = PDB_ORIGX;
                break;

        case 'r':
                if (STREQN(buffer + 1, "ema", 3))
                        r.record_type = PDB_REMARK;
                else if (STREQN(buffer + 1, "evd", 3))
                        r.record_type = PDB_REVDAT;
                break;

        case 's':
                switch (buffer[1]) {

                case 'c':
                        if (STREQN(buffer + 2, "al", 2))
                                r.record_type = PDB_SCALE;
                        break;

                case 'e':
                        if (STREQN(buffer + 2, "qr", 2))
                                r.record_type = PDB_SEQRES;
                        break;

                case 'h':
                        if (STREQN(buffer + 2, "ee", 2))
                                r.record_type = PDB_SHEET;
                        break;

                case 'i':
                        if (STREQN(buffer + 2, "te", 2))
                                r.record_type = PDB_SITE;
                        else if (STREQN(buffer + 2, "ga", 2))
                                r.record_type = PDB_SIGATM;
                        else if (STREQN(buffer + 2, "gu", 2))
                                r.record_type = PDB_SIGUIJ;
                        break;

                case 'o':
                        if (STREQN(buffer + 2, "ur", 2))
                                r.record_type = PDB_SOURCE;
                        break;

                case 'p':
                        if (STREQN(buffer + 2, "rs", 2))
                                r.record_type = PDB_SPRSDE;
                        break;

                case 's':
                        if (STREQN(buffer + 2, "bo", 2))
                                r.record_type = PDB_SSBOND;
                        break;
                }
                break;

        case 't':
                if (STREQN(buffer + 1, "urn", 3))
                        r.record_type = PDB_TURN;
                else if (STREQN(buffer + 1, "vec", 3))
                        r.record_type = PDB_TVECT;
                else if (STREQN(buffer + 1, "er", 2))
                        r.record_type = PDB_TER;
                break;

        case 'u':
                if (STREQN(buffer + 1, "ser", 3))
                        r.record_type = PDB_USER;
                break;
        }

        fmt = pdb_record_format[r.record_type].scan_format;
        switch (r.record_type) {

        case PDB_UNKNOWN:
                pdb_sscanf(buffer, fmt, r.pdb.unknown.junk);
                break;

        case PDB_ANISOU:
        case PDB_SIGUIJ:
                pdb_sscanf(buffer, fmt, &r.pdb.anisou.serial_num,
                        r.pdb.anisou.name, &r.pdb.anisou.alt_loc,
                        r.pdb.anisou.residue.name,
                        &r.pdb.anisou.residue.chain_id,
                        &r.pdb.anisou.residue.seq_num,
                        &r.pdb.anisou.residue.insert_code,
                        &r.pdb.anisou.u[0], &r.pdb.anisou.u[1],
                        &r.pdb.anisou.u[2], &r.pdb.anisou.u[3],
                        &r.pdb.anisou.u[4], &r.pdb.anisou.u[5]);
                break;

        case PDB_ATOM:
        case PDB_HETATM:
        case PDB_SIGATM:
        {
/*
ATOM  41150  H1  WAT  8081      21.760   0.196  43.791  1.00  0.00
ATOM  51150  H1  WAT  10581 
*/
                int atom_overflow = 0, residue_overflow = 0;

                if (isdigit(buffer[11]))
                        atom_overflow = 1;
                if (isdigit(buffer[26]))
                        residue_overflow = 1;
                if (atom_overflow  &&  !residue_overflow) {
                        fmt =
          "%6 %6d%4s%c%3s %c%4d%c   %8f%8f%8f%6f%6f %3d";
                } else if (!atom_overflow  &&  residue_overflow) {
                        fmt =
                  "%6 %5d %4s%c%3s %c%5d%c  %8f%8f%8f%6f%6f %3d";
                } else if (atom_overflow  &&  residue_overflow) {
                        fmt =
                  "%6 %6d%4s%c%3s %c%5d%c  %8f%8f%8f%6f%6f %3d";
                }
                
                pdb_sscanf(buffer, fmt, &r.pdb.atom.serial_num,
                        r.pdb.atom.name, &r.pdb.atom.alt_loc,
                        r.pdb.atom.residue.name,
                        &r.pdb.atom.residue.chain_id,
                        &r.pdb.atom.residue.seq_num,
                        &r.pdb.atom.residue.insert_code,
                        &r.pdb.atom.x, &r.pdb.atom.y, &r.pdb.atom.z,
                        &r.pdb.atom.occupancy, &r.pdb.atom.temp_factor,
                        &r.pdb.atom.ftnote_num);
                break;
        }
        case PDB_AUTHOR:
        case PDB_COMPND:
        case PDB_JRNL:
        case PDB_SOURCE:
                pdb_sscanf(buffer, fmt, &r.pdb.author.continuation,
                        r.pdb.author.data);
                break;

        case PDB_CONECT:
                pdb_sscanf(buffer, fmt, &r.pdb.conect.serial_num,
                        &r.pdb.conect.covalent[0], &r.pdb.conect.covalent[1],
                        &r.pdb.conect.covalent[2], &r.pdb.conect.covalent[3],
                        &r.pdb.conect.bonds[0].hydrogen[0],
                        &r.pdb.conect.bonds[0].hydrogen[1],
                        &r.pdb.conect.bonds[0].salt,
                        &r.pdb.conect.bonds[1].hydrogen[0],
                        &r.pdb.conect.bonds[1].hydrogen[1],
                        &r.pdb.conect.bonds[1].salt);
                break;

        case PDB_CRYST1:
                pdb_sscanf(buffer, fmt, &r.pdb.cryst1.a, &r.pdb.cryst1.b,
                        &r.pdb.cryst1.c, &r.pdb.cryst1.alpha,
                        &r.pdb.cryst1.beta, &r.pdb.cryst1.gamma,
                        r.pdb.cryst1.space_grp, &r.pdb.cryst1.z);
                break;

        case PDB_END:
                break;

        case PDB_FORMUL:
                pdb_sscanf(buffer, fmt, &r.pdb.formul.component,
                        r.pdb.formul.het_id, &r.pdb.formul.continuation,
                        &r.pdb.formul.exclude, r.pdb.formul.formula);
                break;

        case PDB_FTNOTE:
        case PDB_REMARK:
                pdb_sscanf(buffer, fmt, &r.pdb.ftnote.num, r.pdb.ftnote.text);
                break;

        case PDB_HEADER:
                pdb_sscanf(buffer, fmt, r.pdb.header.class, r.pdb.header.date,
                        r.pdb.header.id);
                break;

        case PDB_HELIX:
                pdb_sscanf(buffer, fmt, &r.pdb.helix.serial_num,
                        r.pdb.helix.id,
                        r.pdb.helix.residues[0].name,
                        &r.pdb.helix.residues[0].chain_id,
                        &r.pdb.helix.residues[0].seq_num,
                        &r.pdb.helix.residues[0].insert_code,
                        r.pdb.helix.residues[1].name,
                        &r.pdb.helix.residues[1].chain_id,
                        &r.pdb.helix.residues[1].seq_num,
                        &r.pdb.helix.residues[1].insert_code,
                        &r.pdb.helix.class, r.pdb.helix.comment);
                break;

        case PDB_HET:
                pdb_sscanf(buffer, fmt, r.pdb.het.het_grp.name,
                        &r.pdb.het.het_grp.chain_id, &r.pdb.het.het_grp.seq_num,
                        &r.pdb.het.het_grp.insert_code, &r.pdb.het.num_atoms,
                        r.pdb.het.text);
                break;

        case PDB_MASTER:
                pdb_sscanf(buffer, fmt, &r.pdb.master.num_remark,
                        &r.pdb.master.num_ftnote, &r.pdb.master.num_het,
                        &r.pdb.master.num_helix, &r.pdb.master.num_sheet,
                        &r.pdb.master.num_turn, &r.pdb.master.num_site,
                        &r.pdb.master.num_transform,
                        &r.pdb.master.num_coordinate,
                        &r.pdb.master.num_ter, &r.pdb.master.num_conect,
                        &r.pdb.master.num_seqres);
                break;

        case PDB_MTRIX:
                pdb_sscanf(buffer, fmt, &r.pdb.mtrix.row_num,
                        &r.pdb.mtrix.serial_num, &r.pdb.mtrix.m1,
                        &r.pdb.mtrix.m2, &r.pdb.mtrix.m3, &r.pdb.mtrix.v,
                        &r.pdb.mtrix.given);
                break;

        case PDB_OBSLTE:
                pdb_sscanf(buffer, fmt, &r.pdb.obslte.continuation,
                        r.pdb.obslte.date, r.pdb.obslte.old_id,
                        r.pdb.obslte.id_map[0], r.pdb.obslte.id_map[1],
                        r.pdb.obslte.id_map[2], r.pdb.obslte.id_map[3],
                        r.pdb.obslte.id_map[4], r.pdb.obslte.id_map[2],
                        r.pdb.obslte.id_map[6], r.pdb.obslte.id_map[7]);
                break;

        case PDB_ORIGX:
                pdb_sscanf(buffer, fmt, &r.pdb.origx.row_num, &r.pdb.origx.o1,
                        &r.pdb.origx.o2, &r.pdb.origx.o3, &r.pdb.origx.t);
                break;

        case PDB_REVDAT:
                pdb_sscanf(buffer, fmt, &r.pdb.revdat.modification,
                        &r.pdb.revdat.continuation, r.pdb.revdat.date,
                        r.pdb.revdat.id, &r.pdb.revdat.mod_type,
                        r.pdb.revdat.corrections);
                break;

        case PDB_SCALE:
                pdb_sscanf(buffer, fmt, &r.pdb.scale.row_num, &r.pdb.scale.s1,
                        &r.pdb.scale.s2, &r.pdb.scale.s3, &r.pdb.scale.u);
                break;

        case PDB_SEQRES:
                pdb_sscanf(buffer, fmt, &r.pdb.seqres.serial_num,
                        &r.pdb.seqres.chain_id, &r.pdb.seqres.count,
                        r.pdb.seqres.names[0], r.pdb.seqres.names[1],
                        r.pdb.seqres.names[2], r.pdb.seqres.names[3],
                        r.pdb.seqres.names[4], r.pdb.seqres.names[5],
                        r.pdb.seqres.names[6], r.pdb.seqres.names[7],
                        r.pdb.seqres.names[8], r.pdb.seqres.names[9],
                        r.pdb.seqres.names[10], r.pdb.seqres.names[11],
                        r.pdb.seqres.names[12]);
                break;

        case PDB_SHEET:
                sh = &r.pdb.sheet;
                sha0 = &sh->atoms[0].residue;
                sha1 = &sh->atoms[1].residue;
                pdb_sscanf(buffer, fmt, &r.pdb.sheet.strand_num,
                        sh->id, &r.pdb.sheet.count,
                        sh->residues[0].name, &sh->residues[0].chain_id,
                        &sh->residues[0].seq_num, &sh->residues[0].insert_code,
                        sh->residues[1].name, &sh->residues[1].chain_id,
                        &sh->residues[1].seq_num, &sh->residues[1].insert_code,
                        &sh->sense,
                        sh->atoms[0].name,
                        sha0->name, &sha0->chain_id, &sha0->seq_num,
                        &sha0->insert_code,
                        sh->atoms[1].name,
                        sha1->name, &sha1->chain_id, &sha1->seq_num,
                        &sha1->insert_code);
                break;

        case PDB_SITE:
                pdb_sscanf(buffer, fmt, &r.pdb.site.seq_num,
                        r.pdb.site.id, &r.pdb.site.count,
                        r.pdb.site.residues[0].name,
                        &r.pdb.site.residues[0].chain_id,
                        &r.pdb.site.residues[0].seq_num,
                        &r.pdb.site.residues[0].insert_code,
                        r.pdb.site.residues[1].name,
                        &r.pdb.site.residues[1].chain_id,
                        &r.pdb.site.residues[1].seq_num,
                        &r.pdb.site.residues[1].insert_code,
                        r.pdb.site.residues[2].name,
                        &r.pdb.site.residues[2].chain_id,
                        &r.pdb.site.residues[2].seq_num,
                        &r.pdb.site.residues[2].insert_code,
                        r.pdb.site.residues[3].name,
                        &r.pdb.site.residues[3].chain_id,
                        &r.pdb.site.residues[3].seq_num,
                        &r.pdb.site.residues[3].insert_code);
                break;

        case PDB_SPRSDE:
                pdb_sscanf(buffer, fmt, &r.pdb.sprsde.continuation,
                        r.pdb.sprsde.date, r.pdb.sprsde.id,
                        r.pdb.sprsde.supersede[0], r.pdb.sprsde.supersede[1],
                        r.pdb.sprsde.supersede[2], r.pdb.sprsde.supersede[3],
                        r.pdb.sprsde.supersede[4], r.pdb.sprsde.supersede[5],
                        r.pdb.sprsde.supersede[6], r.pdb.sprsde.supersede[7]);
                break;

        case PDB_SSBOND:
                pdb_sscanf(buffer, fmt, &r.pdb.ssbond.seq_num,
                        r.pdb.ssbond.residues[0].name,
                        &r.pdb.ssbond.residues[0].chain_id,
                        &r.pdb.ssbond.residues[0].seq_num,
                        &r.pdb.ssbond.residues[0].insert_code,
                        r.pdb.ssbond.residues[1].name,
                        &r.pdb.ssbond.residues[1].chain_id,
                        &r.pdb.ssbond.residues[1].seq_num,
                        &r.pdb.ssbond.residues[1].insert_code,
                        r.pdb.ssbond.comment);
                break;

        case PDB_TER:
                pdb_sscanf(buffer, fmt, &r.pdb.ter.serial_num,
                        r.pdb.ter.residue.name, &r.pdb.ter.residue.chain_id,
                        &r.pdb.ter.residue.seq_num,
                        &r.pdb.ter.residue.insert_code);
                break;

        case PDB_TURN:
                pdb_sscanf(buffer, fmt, &r.pdb.turn.seq_num, r.pdb.turn.id,
                        r.pdb.turn.residues[0].name,
                        &r.pdb.turn.residues[0].chain_id,
                        &r.pdb.turn.residues[0].seq_num,
                        &r.pdb.turn.residues[0].insert_code,
                        r.pdb.turn.residues[1].name,
                        &r.pdb.turn.residues[1].chain_id,
                        &r.pdb.turn.residues[1].seq_num,
                        &r.pdb.turn.residues[1].insert_code,
                        r.pdb.turn.comment);
                break;

        case PDB_TVECT:
                pdb_sscanf(buffer, fmt, &r.pdb.tvect.serial_num, &r.pdb.tvect.t1,
                        &r.pdb.tvect.t2, &r.pdb.tvect.t3, r.pdb.tvect.comment);
                break;

        case PDB_USER:
                pdb_sscanf(buffer, fmt, r.pdb.user.subtype, r.pdb.user.text);
                break;

        case PDB_ENDMDL:
                break;

        case PDB_MODEL:
                pdb_sscanf(buffer, fmt, &r.pdb.model.num);
                break;
        }
        
        return r;
}

# ifdef vms
pdb_read_dummy()
{
        pdb_fmt_dummy();
}
# endif
