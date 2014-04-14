/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* XpmWrFFrI.c:                                                                *
*                                                                             *
*  XPM library                                                                *
*  Write an image and possibly its mask to an XPM file                        *
*                                                                             *
*  Developed by Arnaud Le Hors                                                *
\*****************************************************************************/

#include "xpmP.h"
#ifdef VMS
#include "sys$library:string.h"
#else
#include <string.h>
#endif

LFUNC(WriteColors, void, (FILE *file, XpmColor *colors, unsigned int ncolors));

LFUNC(WritePixels, int, (FILE *file, unsigned int width, unsigned int height,
			 unsigned int cpp, unsigned int *pixels,
			 XpmColor *colors));

LFUNC(WriteExtensions, void, (FILE *file, XpmExtension *ext,
			      unsigned int num));

int
XpmWriteFileFromImage(display, filename, image, shapeimage, attributes)
    Display *display;
    char *filename;
    XImage *image;
    XImage *shapeimage;
    XpmAttributes *attributes;
{
    XpmImage xpmimage;
    XpmInfo info;
    int ErrorStatus;

    /* create an XpmImage from the image */
    ErrorStatus = XpmCreateXpmImageFromImage(display, image, shapeimage,
					     &xpmimage, attributes);
    if (ErrorStatus != XpmSuccess)
        return (ErrorStatus);

    /* write the file from the XpmImage */
    if (attributes) {
	xpmSetInfo(&info, attributes);
	ErrorStatus = XpmWriteFileFromXpmImage(filename, &xpmimage, &info);
    } else
	ErrorStatus = XpmWriteFileFromXpmImage(filename, &xpmimage, NULL);

    /* free the XpmImage */
    XpmFreeXpmImage(&xpmimage);

    return (ErrorStatus);
}

int
XpmWriteFileFromXpmImage(filename, image, info)
    char *filename;
    XpmImage *image;
    XpmInfo *info;
{
    xpmData mdata;
    char *name, *dot, *s, new_name[BUFSIZ];
    int ErrorStatus;

    /* open file to write */
    if ((ErrorStatus = xpmWriteFile(filename, &mdata)) != XpmSuccess)
	return (ErrorStatus);

    /* figure out a name */
    if (filename) {
#ifdef VMS
	name = filename;
#else
	if (!(name = rindex(filename, '/')))
	    name = filename;
	else
	    name++;
#endif
	if ((dot = index(name, '.'))) {
	    strcpy(new_name, name);
	    /* change '.' to '_' to get a valid C syntax name */
	    name = s = new_name;
	    while ((dot = index(s, '.'))) {
		*dot = '_';
		s = dot;
	    }
	}
    } else
	name = "image_name";

    /* write the XpmData from the XpmImage */
    if (ErrorStatus == XpmSuccess)
	ErrorStatus = xpmWriteData(&mdata, image, name, info);

    xpmDataClose(&mdata);

    return (ErrorStatus);
}

int
xpmWriteData(mdata, image, name, info)
    xpmData *mdata;
    XpmImage *image;
    char *name;
    XpmInfo *info;
{
    /* calculation variables */
    unsigned int cmts, extensions;
    FILE *file;
    int ErrorStatus;

    /* store this to speed up */
    file = mdata->stream.file;

    cmts = info && (info->valuemask & XpmComments);
    extensions = info && (info->valuemask & XpmExtensions)
	&& info->nextensions;

    /* print the header line */
    fprintf(file, "/* XPM */\nstatic char * %s[] = {\n", name);

    /* print the hints line */
    if (cmts && info->hints_cmt)
	fprintf(file, "/*%s*/\n", info->hints_cmt);

    fprintf(file, "\"%d %d %d %d", image->width, image->height,
	    image->ncolors, image->cpp);

    if (info && (info->valuemask & XpmHotspot))
	fprintf(file, " %d %d", info->x_hotspot, info->y_hotspot);

    if (extensions)
	fprintf(file, " XPMEXT");

    fprintf(file, "\",\n");

    /* print colors */
    if (cmts && info->colors_cmt)
	fprintf(file, "/*%s*/\n", info->colors_cmt);

    WriteColors(file, image->colorTable, image->ncolors);

    /* print pixels */
    if (cmts && info->pixels_cmt)
	fprintf(file, "/*%s*/\n", info->pixels_cmt);

    ErrorStatus = WritePixels(file, image->width, image->height, image->cpp,
			      image->data, image->colorTable);
    if (ErrorStatus != XpmSuccess)
	return (ErrorStatus);

    /* print extensions */
    if (extensions)
	WriteExtensions(file, info->extensions, info->nextensions);

    /* close the array */
    fprintf(file, "};\n");

    return (XpmSuccess);
}

static void
WriteColors(file, colors, ncolors)
    FILE *file;
    XpmColor *colors;
    unsigned int ncolors;
{
    unsigned int a, key;
    char *s;
    char **defaults;

    for (a = 0; a < ncolors; a++, colors++) {

	defaults = (char **) colors;
	fprintf(file, "\"%s", *defaults++);

	for (key = 1; key <= NKEYS; key++, defaults++) {
	    if ((s = *defaults))
		fprintf(file, "\t%s %s", xpmColorKeys[key - 1], s);
	}
	fprintf(file, "\",\n");
    }
}


static int
WritePixels(file, width, height, cpp, pixels, colors)
    FILE *file;
    unsigned int width;
    unsigned int height;
    unsigned int cpp;
    unsigned int *pixels;
    XpmColor *colors;
{
    char *s, *p, *buf;
    unsigned int x, y, h;

    h = height - 1;
    p = buf = (char *) XpmMalloc(width * cpp + 3);
    if (!buf)
	return (XpmNoMemory);
    *buf = '"';
    p++;
    for (y = 0; y < h; y++) {
	s = p;
	for (x = 0; x < width; x++, pixels++) {
	    strncpy(s, colors[*pixels].string, cpp);
	    s += cpp;
	}
	*s++ = '"';
	*s = '\0';
	fprintf(file, "%s,\n", buf);
    }
    /* duplicate some code to avoid a test in the loop */
    s = p;
    for (x = 0; x < width; x++, pixels++) {
	strncpy(s, colors[*pixels].string, cpp);
	s += cpp;
    }
    *s++ = '"';
    *s = '\0';
    fprintf(file, "%s", buf);

    XpmFree(buf);
    return (XpmSuccess);
}

static void
WriteExtensions(file, ext, num)
    FILE *file;
    XpmExtension *ext;
    unsigned int num;
{
    unsigned int x, y, n;
    char **line;

    for (x = 0; x < num; x++, ext++) {
	fprintf(file, ",\n\"XPMEXT %s\"", ext->name);
	n = ext->nlines;
	for (y = 0, line = ext->lines; y < n; y++, line++)
	    fprintf(file, ",\n\"%s\"", *line);
    }
    fprintf(file, ",\n\"XPMENDEXT\"");
}
