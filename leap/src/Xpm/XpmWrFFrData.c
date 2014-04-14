/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmWrFFrData.c:                                                             *
*                                                                             *
*  XPM library                                                                *
*  Parse an Xpm array and write a file that corresponds to it.                *
*                                                                             *
*  Developed by Dan Greening dgreen@cs.ucla.edu / dgreen@sti.com              *
\*****************************************************************************/

#include "xpmP.h"

int
XpmWriteFileFromData(filename, data)
    char *filename;
    char **data;
{
    XpmImage image;
    XpmInfo info;
    int ErrorStatus;

    info.valuemask = XpmReturnComments | XpmReturnExtensions;

    ErrorStatus = XpmCreateXpmImageFromData(data, &image, &info);

    if (ErrorStatus != XpmSuccess)
	return (ErrorStatus);

    ErrorStatus = XpmWriteFileFromXpmImage(filename, &image, &info);

    XpmFreeXpmImage(&image);
    XpmFreeXpmInfo(&info);

    return (ErrorStatus);
}
