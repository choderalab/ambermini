/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmCrDataFP.c:                                                              *
*                                                                             *
*  XPM library                                                                *
*  Scan a pixmap and possibly its mask and create an XPM array                *
*                                                                             *
*  Developed by Arnaud Le Hors                                                *
\*****************************************************************************/

#include "xpmP.h"
#ifdef VMS
#include "sys$library:string.h"
#else
#if defined(SYSV) || defined(SVR4)
#include <string.h>
#else
#include <strings.h>
#endif
#endif

int
XpmCreateDataFromPixmap(display, data_return, pixmap, shapemask, attributes)
    Display *display;
    char ***data_return;
    Pixmap pixmap;
    Pixmap shapemask;
    XpmAttributes *attributes;
{
    XImage *ximage = NULL;
    XImage *shapeimage = NULL;
    unsigned int width = 0;
    unsigned int height = 0;
    int ErrorStatus;

    /* get geometry */
    if (attributes && attributes->valuemask & XpmSize) {
	width = attributes->width;
	height = attributes->height;
    }

    /* get the ximages */
    if (pixmap)
	xpmCreateImageFromPixmap(display, pixmap, &ximage, &width, &height);
    if (shapemask)
	xpmCreateImageFromPixmap(display, shapemask, &shapeimage,
				 &width, &height);

    /* create the data */
    ErrorStatus = XpmCreateDataFromImage(display, data_return, ximage,
					 shapeimage, attributes);

    /* destroy the ximages */
    if (ximage)
	XDestroyImage(ximage);
    if (shapeimage)
	XDestroyImage(shapeimage);

    return (ErrorStatus);
}
