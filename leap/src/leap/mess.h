/*
 *	File:	message.h
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
 *	Description:
 *		Store messages in a long list.
 *		Keep track of what function they come from
 *		using 'function.h'.
 */


#ifndef	BASICS_H
#include	"basics.h"
#endif

typedef	struct MESS_STRUCT	{
	struct MESS_STRUCT	*mNext;
	char			cType;
	int			iFunction;
	char			*cPText;
	BOOL			bPrint;
} MESSt;

typedef	MESSt	*MESS;


extern MESS	GmMesss;
extern MESS	GmLastMess;

extern MESS	mMessAdd(char *sText);
extern MESS	mMessLoop();
extern MESS	mMessNext(MESS *mPMess);

#define	sMessText(m)		(m->cPText)
#define	bMessPrint(m)		(m->bPrint)
#define	iMessFunction(m)	(m->iFunction)
#define	MessSetPrint(m,b)	(m->bPrint = b)

