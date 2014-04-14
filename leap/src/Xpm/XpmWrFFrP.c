/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmWrFFrP.c:                                                                *
*                                                                             *
*  XPM library                                                                *
*  Write a pixmap and possibly its mask to an XPM file                        *
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
XpmWriteFileFromPixmap(display, filename, pixmap, shapemask, attributes)
    Display *display;
    char *filename;
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

    /* write to the file */
    ErrorStatus = XpmWriteFileFromImage(display, filename, ximage, shapeimage,
					attributes);

    /* destroy the ximages */
    if (ximage)
	XDestroyImage(ximage);
    if (shapeimage)
	XDestroyImage(shapeimage);

     return (ErrorStatus);
}
