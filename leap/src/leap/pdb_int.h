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
 *      WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *      $Id: pdb_int.h,v 10.0 2008/04/15 23:22:21 case Exp $
 */

# include       "pdb.h"

# define        STREQN(s1, s2, n)       (strncmp(s1, s2, n) == 0)

struct  pdb_format      {
        char    *scan_format;
        char    *print_format;
};

extern  struct  pdb_format      pdb_record_format[];

extern  void    pdb_sscanf(char *, char *, ...);
extern  void    pdb_sprintf(char *, char *, ...);
