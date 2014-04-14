/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmRdFToP.c:                                                                *
*                                                                             *
*  XPM library                                                                *
*  Parse an XPM file and create the pixmap and possibly its mask              *
*                                                                             *
*  Developed by Arnaud Le Hors                                                *
\*****************************************************************************/

#include "xpmP.h"

int
XpmReadFileToPixmap(display, d, filename, pixmap_return,
		    shapemask_return, attributes)
    Display *display;
    Drawable d;
    char *filename;
    Pixmap *pixmap_return;
    Pixmap *shapemask_return;
    XpmAttributes *attributes;
{
    XImage *ximage, *shapeimage;
    int ErrorStatus;

    /* initialize return values */
    if (pixmap_return)
	*pixmap_return = 0;
    if (shapemask_return)
	*shapemask_return = 0;

    /* create the images */
    ErrorStatus = XpmReadFileToImage(display, filename, &ximage, &shapeimage,
				     attributes);

    if (ErrorStatus < 0)	/* fatal error */
        return (ErrorStatus);

    /* create the pixmaps and destroy images */
     if (ximage) {
	xpmCreatePixmapFromImage(display, d, ximage, pixmap_return);
	XDestroyImage(ximage);
    }
    if (shapeimage) {
	xpmCreatePixmapFromImage(display, d, shapeimage, shapemask_return);
	XDestroyImage(shapeimage);
    }

    return (ErrorStatus);
}
