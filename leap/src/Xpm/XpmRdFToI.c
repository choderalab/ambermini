/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmRdFToI.c:                                                                *
*                                                                             *
*  XPM library                                                                *
*  Parse an XPM file and create the image and possibly its mask               *
*                                                                             *
*  Developed by Arnaud Le Hors                                                *
\*****************************************************************************/

#include "xpmP.h"

int
XpmReadFileToImage(display, filename,
		   image_return, shapeimage_return, attributes)
    Display *display;
    char *filename;
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
	ErrorStatus = XpmReadFileToXpmImage(filename, &image, &info);
    } else
	ErrorStatus = XpmReadFileToXpmImage(filename, &image, NULL);

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
XpmReadFileToXpmImage(filename, image, info)
    char *filename;
    XpmImage *image;
    XpmInfo *info;
{
    xpmData mdata;
    int ErrorStatus;

    /* init returned values */
    xpmInitXpmImage(image);
    xpmInitXpmInfo(info);

    /* open file to read */
    if ((ErrorStatus = xpmReadFile(filename, &mdata)) != XpmSuccess)
	return (ErrorStatus);

    /* create the XpmImage from the XpmData */
    ErrorStatus = xpmParseData(&mdata, image, info);

    xpmDataClose(&mdata);

    return (ErrorStatus);
}
