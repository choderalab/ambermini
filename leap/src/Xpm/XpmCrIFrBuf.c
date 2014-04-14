/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmCrIFrBuf.c:                                                              *
*                                                                             *
*  XPM library                                                                *
*  Parse an Xpm buffer (file in memory) and create the image and possibly its *
*  mask                                                                       *
*  Developed by Arnaud Le Hors                                                *
\*****************************************************************************/

#include "xpmP.h"

int
XpmCreateImageFromBuffer(display, buffer, image_return,
			 shapeimage_return, attributes)
    Display *display;
    char *buffer;
    XImage **image_return;
    XImage **shapeimage_return;
    XpmAttributes *attributes;
{
    XpmImage image;
    XpmInfo info;
    int ErrorStatus;

    /* create an XpmImage from the buffer */
    if (attributes) {
	xpmInitAttributes(attributes);
	xpmSetInfoMask(&info, attributes);
	ErrorStatus = XpmCreateXpmImageFromBuffer(buffer, &image, &info);
    } else
	ErrorStatus = XpmCreateXpmImageFromBuffer(buffer, &image, NULL);

    if (ErrorStatus != XpmSuccess)
        return (ErrorStatus);

    /* create the related ximages */
    ErrorStatus = XpmCreateImageFromXpmImage(display, &image,
					     image_return, shapeimage_return,
					     attributes);
    if (attributes) {
	if (ErrorStatus >= 0)	/* no fatal error */
	    xpmSetAttributes(attributes, &image, &info);
	XpmFreeXpmInfo(&info);
    }

    /* free the XpmImage */
    XpmFreeXpmImage(&image);

    return (ErrorStatus);
}

int
XpmCreateXpmImageFromBuffer(buffer, image, info)
    char *buffer;
    XpmImage *image;
    XpmInfo *info;
{
    xpmData mdata;
    int ErrorStatus;

    /* init returned values */
    xpmInitXpmImage(image);
    xpmInitXpmInfo(info);

    /* open buffer to read */
    xpmOpenBuffer(buffer, &mdata);

    /* create the XpmImage from the XpmData */
    ErrorStatus = xpmParseData(&mdata, image, info);

    xpmDataClose(&mdata);

    return (ErrorStatus);
}
