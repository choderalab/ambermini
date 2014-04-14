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
 *      WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *      $Id: pdb_format.c,v 10.1 2008/10/16 11:20:23 case Exp $
 */

/*      Modifications induced by the implementation of the savemol2 command
*        Christine Cezard (2007) 
*       Universite de Picardie - Jules Verne, Amiens
*        http://q4md-forcefieldtools.org
*
*       declaration of a MOL2_ATOM label line 
 */    
 
  
/* LINTLIBRARY */


# include       "pdb_int.h"
/*
 *      for each pdb record type there is a format reading in the
 *      record values and for printing them out.
 *
 *      The actual format of a line written, is the print format
 *      followed by blank padding to 72 characters, followed by
 *      8 characters of file and line information.
 */

struct  pdb_format      pdb_record_format[PDB_NUM_R]    = {
        {                                       /* 0 PDB_UNKNOWN */
                "%80s", "UNKNOWN:  ??%-6.6s??" },
        {                                       /* 1 PDB_ANISOU, SIGUIJ */
                "%6 %5d %4s%c%3s %c%4d%c %7d%7d%7d%7d%7d%7d",
                "ANISOU%5d %-4s%C%-3s %C%4d%C %7d%7d%7d%7d%7d%7d"
        },
        {                                       /* 2 PDB_ATOM, HETATM, SIGATM */
                "%6 %5d %4s%c%3s %c%4d%c   %8f%8f%8f%6f%6f %3d",
                "ATOM  %5d %-4s%C%-3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D"
        },
{                                       /* 3 PDB_AUTHOR, COMPND, JRNL, SOURCE */
                "%9 %c%60s", "AUTHOR   %C%-60s" },
        {                                       /* 4 PDB_COMPND */
                "%9 %c%60s", "COMPND   %C%-60s" },
        {                                       /* 5 PDB_CONECT */
                "%6 %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",
                "%10d%10d%10d%10d" },
        {                                       /* 6 PDB_CRYST1 */
                "%6 %9f%9f%9f%7f%7f%7f %11s%4d",
                "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d" },
        {                                       /* 7 PDB_END */
                "", "END   " },
        {                                       /* 8 PDB_FORMUL */
                "%8 %2d  %3s %2d%c%51s", "FORMUL  %2D  %-3s %2D%C%-51s" },
        {                                       /* 9 PDB_FTNOTE, REMARK */
                "%7 %3d %59s", "FTNOTE %3D %-59s" },
        {                                       /* 10 PDB_HEADER */
                "%10 %40s%9s%3 %4s", "HEADER    %-40s%-12s%-4s" },
        {                                       /* 11 PDB_HELIX */
                "%7 %3d %3s %3s %c %4d%c %3s %c %4d%c%2d%30s",
                "HELIX  %3D %3s %-3s %C %4d%C %-3s %C %4d%C%2d%-30s" },
        {                                       /* 12 PDB_HET */
                "%7 %3s  %c%4d%c  %5d%5 %40s",
                "HET    %-3s  %C%4d%C  %5d     %-40s" },
        {                                       /* 13 PDB_HETATM */
                "%6 %5d %4s%c%3s %c%4d%c   %8f%8f%8f%6f%6f %3d",
                "HETATM%5d %-4s%C%-3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D"
        },
        {                                       /* 14 PDB_JRNL */
                "%9 %c%60s", "JRNL     %C%-60s" },
        {                                       /* 15 PDB_MASTER */
                "%10 %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",
                "MASTER    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d" },
        {                                       /* 16 PDB_MTRIX */
                "%5 %d %3d%10f%10f%10f%5 %10f   %2d",
                "MTRIX%d %3d%10.5f%10.5f%10.5f     %10.5f   %2D" },
        {                                       /* 17 PDB_OBSLTE */
                "%8 %2d %9s %4s%6 %4s %4s %4s %4s %4s %4s %4s %4s",
                "OBSLTE  %2D %-9s %-10s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-4s" },
        {                                       /* 18 PDB_ORIGX, SCALE */
                "%5 %d%4 %10f%10f%10f%5 %10f",
                "ORIGX%d    %10.5f%10.5f%10.5f     %10.5f" },
        {                                       /* 19 PDB_REMARK */
                "%7 %3d %59s", "REMARK %3D %-59s" },
        {                                       /* 20 PDB_REVDAT */
                "%7 %3d%2d %9s %7s %c%7 %31s",
                "REVDAT %3D%2D %-9s %-7s %c       %-31s" },
        {                                       /* 21 PDB_SCALE */
                "%5 %d%4 %10f%10f%10f%5 %10f",
                "SCALE%d    %10.5f%10.5f%10.5f     %10.5f" },
        {                                       /* 22 PDB_SEQRES */
        "%6 %4d %c %4d  %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s",
        "SEQRES%4d %C %4d  %-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-3s"
        ,},
        {                                       /* 23 PDB_SHEET */
"%6 %4d %3s%2d %3s %c%4d%c %3s %c%4d%c%2d %4s%3s %c%4d%c %4s%3s %c%4d%c",
"SHEET %4D %3s%2d %-3s %C%4d%C %-3s %C%4d%C%2d %-4s%-3s %C%4D%C %-4s%-3s %C%4D%C"
        },
        {                                       /* 24 PDB_SIGATM */
                "%6 %5d %4s%c%3s %c%4d%c   %8f%8f%8f%6f%6f %3d",
                "SIGATM%5d %-4s%C%-3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D"
        },
        {                                       /* 25 PDB_SIGUIJ */
                "%6 %5d %4s%c%3s %c%4d%c %7d%7d%7d%7d%7d%7d",
                "SIGUIJ%5d %-4s%C%-3s %C%4d%C %7D%7D%7D%7D%7D%7D"
        },
        {                                       /* 26 PDB_SITE */
        "%7 %3d %3s %2d %3s %c%4d%c %3s %c%4d%c %3s %c%4d%c %3s %c%4d%c",
        "SITE   %3d %3s %2d %-3s %C%4D%C %-3s %C%4D%C %-3s %C%4D%C %-3s %C%4D%C"
        },
        {                                       /* 27 PDB_SOURCE */
                "%9 %c%60s", "SOURCE   %C%-60s" },
        {                                       /* 28 PDB_SPRSDE */
                "%8 %2d %9s %4s%6 %4s %4s %4s %4s %4s %4s %4s %4s",
                "SPRSDE  %2D %-9s %-10s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-4s" },
        {                                       /* 29 PDB_SSBOND */
                "%7 %3d %3s %c %4d%c   %3s %c %4d%c%4 %30s",
                "SSBOND %3D %3s %C %4d%C   %3s %C %4D%C    %-30s" },
        {                                       /* 30 PDB_TER */
                "%6 %5d%6 %3s %c%4d%c",
                "TER   %5d      %-3s %C%4d%C" },
        {                                       /* 31 PDB_TURN */
                "%7 %3d %3s %3s %c%4d%c %3s %c%4d%c%4 %30s",
                "TURN   %3D %3s %-3s %C%4d%C %-3s %C%4d%C    %-30s" },
        {                                       /* 32 PDB_TVECT */
                "%7 %3d%10f%10f%10f%30s",
                "TVECT  %3D%10.5f%10.5f%10.5f%-30s" },
        {                                       /* 33 PDB_USER */
                "%4 %2s%66s", "USER%-2s%-66s" },
        {                                       /* 34 PDB_MODEL */
                "%9 %5d", "MODEL    %5d" },
        {                                       /* 35 PDB_ENDMDL */
                "", "ENDMDL" },
        {                                       /* MOL2_ATOM */
                "%3d %6s   %12.6f%12.6f%12.6f %3s %3d %3s %12.4f", 
                "%3d %-4s %11.6f %11.6f %11.6f %-4s %2d %-5s %8.4f ****" },
};

