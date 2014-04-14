/* 
 *      File:   help.h
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              This object manages the database of help information.
 *              Help is stored as a DICTIONARY of keywords which
 *              are associated with very large strings which contain
 *              the help information for that keyword.
 *              There is only one HELP object per system.
 */


#ifndef	HELP_H
#define	HELP_H	

typedef struct  {
	char	*sUpSubject;
	char	*sSubject;
	char	*sText;
} HELPt;

typedef HELPt*  HELP;



/*
 *      Functions
 */

extern void	HelpInitialize();
extern void	HelpShutdown();
extern void     HTInit();
extern void	HelpLoop();
extern HELP	hHelpNext();

                /* sHelpText returns pointers to VERY VERY LONG STRINGS */

extern HELP	hHelp(char *sSubject);

#define sHelpSubject(h)         ((HELP)(h))->sSubject
#define sHelpUpSubject(h)       ((HELP)(h))->sUpSubject
#define sHelpText(h)            ((HELP)(h))->sText

                /* HelpAdd should not be called by the user */
                /* It should only be called from within helptext.c */
extern void	HelpAdd(char *sUpSubject, char *sText);


#endif
