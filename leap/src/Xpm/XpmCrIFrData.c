/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmCrIFrData.c:                                                             *
*                                                                             *
*  XPM library                                                                *
*  Parse an Xpm array and create the image and possibly its mask              *
*                                                                             *
*  Developed by Arnaud Le Hors                                                *
\*****************************************************************************/

#include "xpmP.h"

int
XpmCreateImageFromData(display, data, image_return,
		       shapeimage_return, attributes)
    Display *display;
    char **data;
    XImage **image_return;
    XImage **shapeimage_return;
    XpmAttributes *attributes;
{
    XpmImage image;
    XpmInfo info;
    int ErrorStatus;

    /* create an XpmImage from the file */
    if (attributes) {
	xpmInitAttributes(attributes);
	xpmSetInfoMask(&info, attributes);
	ErrorStatus = XpmCreateXpmImageFromData(data, &image, &info);
    } else
	ErrorStatus = XpmCreateXpmImageFromData(data, &image, NULL);

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

    XpmFreeXpmImage(&image);

    return (ErrorStatus);
}

int
XpmCreateXpmImageFromData(data, image, info)
    char **data;
    XpmImage *image;
    XpmInfo *info;
{
    xpmData mdata;
    int ErrorStatus;

    /* init returned values */
    xpmInitXpmImage(image);
    xpmInitXpmInfo(info);

    /* open data */
    xpmOpenArray(data, &mdata);

    /* create the XpmImage from the XpmData */
    ErrorStatus = xpmParseData(&mdata, image, info);

    xpmDataClose(&mdata);

    return (ErrorStatus);
}

