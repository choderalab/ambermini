/* Copyright 1989-94 GROUPE BULL -- See license conditions in file COPYRIGHT */
/*****************************************************************************\
* misc.c:                                                                     *
*                                                                             *
*  XPM library                                                               *
*  Miscellaneous utilities                                                    *
*                                                                             *
*  Developed by Arnaud Le Hors                                                *
\*****************************************************************************/

#include "xpmP.h"
#ifdef VMS
#include "sys$library:stat.h"
#include "sys$library:fcntl.h"
#include "sys$library:string.h"
#else
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#endif

/* 3.2 backward compatibility code */
LFUNC(CreateOldColorTable, int, (XpmColor *ct, int ncolors,
				 XpmColor ***oldct));

LFUNC(FreeOldColorTable, void, (XpmColor **colorTable, int ncolors));

/*
 * Create a colortable compatible with the old style colortable
 */
static int
CreateOldColorTable(ct, ncolors, oldct)
    XpmColor *ct;
    int ncolors;
    XpmColor ***oldct;
{
    XpmColor **colorTable, **color;
    int a;

    colorTable = (XpmColor **) XpmMalloc(ncolors * sizeof(XpmColor *));
    if (!colorTable) {
	*oldct = NULL;
	return (XpmNoMemory);
    }
    for (a = 0, color = colorTable; a < ncolors; a++, color++, ct++)
	*color = ct;
    *oldct = colorTable;
    return (XpmSuccess);
}

static void
FreeOldColorTable(colorTable, ncolors)
    XpmColor **colorTable;
    int ncolors;
{
    int a, b;
    XpmColor **color;
    char **sptr;

    if (colorTable) {
	for (a = 0, color = colorTable; a < ncolors; a++, color++) {
	    for (b = 0, sptr = (char **) *color; b <= NKEYS; b++, sptr++)
		if (*sptr)
		    XpmFree(*sptr);
	}
	XpmFree(*colorTable);
	XpmFree(colorTable);
    }
}

/* end 3.2 bc */


/*
 * Free the computed color table
 */
void
xpmFreeColorTable(colorTable, ncolors)
    XpmColor *colorTable;
    int ncolors;
{
    int a, b;
    XpmColor *color;
    char **sptr;

    if (colorTable) {
	for (a = 0, color = colorTable; a < ncolors; a++, color++) {
	    for (b = 0, sptr = (char **) color; b <= NKEYS; b++, sptr++)
		if (*sptr)
		    XpmFree(*sptr);
	}
	XpmFree(colorTable);
    }
}

/*
 * Free array of extensions
 */
void
XpmFreeExtensions(extensions, nextensions)
    XpmExtension *extensions;
    int nextensions;
{
    unsigned int i, j, nlines;
    XpmExtension *ext;
    char **sptr;

    if (extensions) {
	for (i = 0, ext = extensions; i < nextensions; i++, ext++) {
	    if (ext->name)
		XpmFree(ext->name);
	    nlines = ext->nlines;
	    for (j = 0, sptr = ext->lines; j < nlines; j++, sptr++)
		if (*sptr)
		    XpmFree(*sptr);
	    if (ext->lines)
		XpmFree(ext->lines);
	}
	XpmFree(extensions);
    }
}


/*
 * Return the XpmAttributes structure size
 */

int 
XpmAttributesSize()
{
    return sizeof(XpmAttributes);
}

/*
 * Init returned data to free safely later on
 */
void
xpmInitAttributes(attributes)
    XpmAttributes *attributes;
{
    if (attributes) {
	attributes->pixels = NULL;
	attributes->npixels = 0;
	attributes->colorTable = NULL;
	attributes->ncolors = 0;
/* 3.2 backward compatibility code */
	attributes->hints_cmt = NULL;
	attributes->colors_cmt = NULL;
	attributes->pixels_cmt = NULL;
/* end 3.2 bc */
	attributes->extensions = NULL;
	attributes->nextensions = 0;
    }
}

/*
 * Fill in the XpmAttributes with the XpmImage and the XpmInfo
 */
void
xpmSetAttributes(attributes, image, info)
    XpmAttributes *attributes;
    XpmImage *image;
    XpmInfo *info;
{
    if (attributes->valuemask & XpmReturnColorTable) {
	attributes->colorTable = image->colorTable;
	attributes->ncolors = image->ncolors;

	/* avoid deletion of copied data */
	image->ncolors = 0;
	image->colorTable = NULL;
    }
/* 3.2 backward compatibility code */
    else if (attributes->valuemask & XpmReturnInfos) {
	int ErrorStatus;

	ErrorStatus = CreateOldColorTable(image->colorTable, image->ncolors,
					  (XpmColor ***)
					  &attributes->colorTable);

	/* if error just say we can't return requested data */
	if (ErrorStatus != XpmSuccess) {
	    attributes->valuemask &= ~XpmReturnInfos;
	    if (!(attributes->valuemask & XpmReturnPixels)) {
		XpmFree(attributes->pixels);
		attributes->pixels = NULL;
		attributes->npixels = 0;
	    }
	    attributes->ncolors = 0;
	} else {
	    attributes->ncolors = image->ncolors;
	    attributes->hints_cmt = info->hints_cmt;
	    attributes->colors_cmt = info->colors_cmt;
	    attributes->pixels_cmt = info->pixels_cmt;

	    /* avoid deletion of copied data */
	    image->ncolors = 0;
	    image->colorTable = NULL;
	    info->hints_cmt = NULL;
	    info->colors_cmt = NULL;
	    info->pixels_cmt = NULL;
	}
    }
/* end 3.2 bc */
    if (attributes->valuemask & XpmReturnExtensions) {
	attributes->extensions = info->extensions;
	attributes->nextensions = info->nextensions;

	/* avoid deletion of copied data */
	info->extensions = NULL;
	info->nextensions = 0;
    }
    if (info->valuemask & XpmHotspot) {
	attributes->valuemask |= XpmHotspot;
	attributes->x_hotspot = info->x_hotspot;
	attributes->y_hotspot = info->y_hotspot;
    }
    attributes->valuemask |= XpmCharsPerPixel;
    attributes->cpp = image->cpp;
    attributes->valuemask |= XpmSize;
    attributes->width = image->width;
    attributes->height = image->height;
}

/*
 * Free the XpmAttributes structure members
 * but the structure itself
 */
void
XpmFreeAttributes(attributes)
    XpmAttributes *attributes;
{
    if (attributes->valuemask & XpmReturnPixels && attributes->npixels) {
	XpmFree(attributes->pixels);
	attributes->pixels = NULL;
	attributes->npixels = 0;
    }
    if (attributes->valuemask & XpmReturnColorTable) {
	xpmFreeColorTable(attributes->colorTable, attributes->ncolors);
	attributes->colorTable = NULL;
	attributes->ncolors = 0;
    }
/* 3.2 backward compatibility code */
    else if (attributes->valuemask & XpmInfos) {
	if (attributes->colorTable) {
	    FreeOldColorTable((XpmColor **) attributes->colorTable,
			      attributes->ncolors);
	    attributes->colorTable = NULL;
	    attributes->ncolors = 0;
	}
	if (attributes->hints_cmt) {
	    XpmFree(attributes->hints_cmt);
	    attributes->hints_cmt = NULL;
	}
	if (attributes->colors_cmt) {
	    XpmFree(attributes->colors_cmt);
	    attributes->colors_cmt = NULL;
	}
	if (attributes->pixels_cmt) {
	    XpmFree(attributes->pixels_cmt);
	    attributes->pixels_cmt = NULL;
	}
	if (attributes->pixels) {
	    XpmFree(attributes->pixels);
	    attributes->pixels = NULL;
	    attributes->npixels = 0;
	}
    }
/* end 3.2 bc */
    if (attributes->valuemask & XpmReturnExtensions
	&& attributes->nextensions) {
	XpmFreeExtensions(attributes->extensions, attributes->nextensions);
	attributes->extensions = NULL;
	attributes->nextensions = 0;
    }
    attributes->valuemask = 0;
}

/*
 * Init returned data to free safely later on
 */
void
xpmInitXpmImage(image)
    XpmImage *image;
{
    image->ncolors = 0;
    image->colorTable = NULL;
    image->data = NULL;
}

/*
 * Free the XpmImage data which have been allocated
 */
void
XpmFreeXpmImage(image)
    XpmImage *image;
{
    if (image->colorTable)
	xpmFreeColorTable(image->colorTable, image->ncolors);
    XpmFree(image->data);
    image->data = NULL;
}

/*
 * Init returned data to free safely later on
 */
void
xpmInitXpmInfo(info)
    XpmInfo *info;
{
    if (info) {
	info->hints_cmt = NULL;
	info->colors_cmt = NULL;
	info->pixels_cmt = NULL;
	info->extensions = NULL;
	info->nextensions = 0;
    }
}

/*
 * Free the XpmInfo data which have been allocated
 */
void
XpmFreeXpmInfo(info)
    XpmInfo *info;
{
    if (info) {
	if (info->valuemask & XpmComments) {
	    if (info->hints_cmt) {
		XpmFree(info->hints_cmt);
		info->hints_cmt = NULL;
	    }
	    if (info->colors_cmt) {
		XpmFree(info->colors_cmt);
		info->colors_cmt = NULL;
	    }
	    if (info->pixels_cmt) {
		XpmFree(info->pixels_cmt);
		info->pixels_cmt = NULL;
	    }
	}
	if (info->valuemask & XpmReturnExtensions && info->nextensions) {
	    XpmFreeExtensions(info->extensions, info->nextensions);
	    info->extensions = NULL;
	    info->nextensions = 0;
	}
	info->valuemask = 0;
    }
}

/*
 * Set the XpmInfo valuemask to retrieve required info
 */
void
xpmSetInfoMask(info, attributes)
    XpmInfo *info;
    XpmAttributes *attributes;
{
    info->valuemask = 0;
    if (attributes->valuemask & XpmReturnInfos)
	info->valuemask |= XpmReturnComments;
    if (attributes->valuemask & XpmReturnExtensions)
	info->valuemask |= XpmReturnExtensions;
}

/*
 * Fill in the XpmInfo with the XpmAttributes
 */
void
xpmSetInfo(info, attributes)
    XpmInfo *info;
    XpmAttributes *attributes;
{
    info->valuemask = 0;
    if (attributes->valuemask & XpmInfos) {
	info->valuemask |= XpmComments | XpmColorTable;
	info->hints_cmt = attributes->hints_cmt;
	info->colors_cmt = attributes->colors_cmt;
	info->pixels_cmt = attributes->pixels_cmt;
    }
    if (attributes->valuemask & XpmExtensions) {
	info->valuemask |= XpmExtensions;
	info->extensions = attributes->extensions;
	info->nextensions = attributes->nextensions;
    }
    if (attributes->valuemask & XpmHotspot) {
	info->valuemask |= XpmHotspot;
	info->x_hotspot = attributes->x_hotspot;
	info->y_hotspot = attributes->y_hotspot;
    }
}


#ifdef NEED_STRDUP
/*
 * in case strdup is not provided by the system here is one
 * which does the trick
 */
char *
strdup(s1)
    char *s1;
{
    char *s2;
    int l = strlen(s1) + 1;

    if (s2 = (char *) XpmMalloc(l))
	strncpy(s2, s1, l);
    return s2;
}

#endif

/*
 *  File / Buffer utilities
 */
int
XpmReadFileToBuffer(filename, buffer_return)
    char *filename;
    char **buffer_return;
{
    int fd, fcheck, len;
    char *ptr;
    struct stat stats;
    FILE *fp;

    *buffer_return = NULL;

    fd = open(filename, O_RDONLY);
    if (fd < 0)
	return XpmOpenFailed;

    if (fstat(fd, &stats)) {
	close(fd);
	return XpmOpenFailed;
    }
    fp = fdopen(fd, "r");
    if (!fp) {
	close(fd);
	return XpmOpenFailed;
    }
    len = (int) stats.st_size;
    ptr = (char *) XpmMalloc(len + 1);
    if (!ptr) {
	fclose(fp);
	return XpmNoMemory;
    }
    fcheck = fread(ptr, len, 1, fp);
    fclose(fp);
    if (fcheck != 1) {
	XpmFree(ptr);
	return XpmOpenFailed;
    }
    ptr[len] = '\0';
    *buffer_return = ptr;
    return XpmSuccess;
}

int
XpmWriteFileFromBuffer(filename, buffer)
    char *filename;
    char *buffer;
{
    int fcheck, len;
    FILE *fp = fopen(filename, "w");

    if (!fp)
	return XpmOpenFailed;

    len = strlen(buffer);
    fcheck = fwrite(buffer, len, 1, fp);
    fclose(fp);
    if (fcheck != 1)
	return XpmOpenFailed;

    return XpmSuccess;
}


/*
 * Small utility function
 */
char *
XpmGetErrorString(errcode)
    int errcode;
{
    switch (errcode) {
    case XpmColorError:
	return ("XpmColorError");
    case XpmSuccess:
	return ("XpmSuccess");
    case XpmOpenFailed:
	return ("XpmOpenFailed");
    case XpmFileInvalid:
	return ("XpmFileInvalid");
    case XpmNoMemory:
	return ("XpmNoMemory");
    case XpmColorFailed:
	return ("XpmColorFailed");
    default:
	return ("Invalid XpmError");
    }
}

/*
 * The following function provides a way to figure out if the linked library is
 * newer or older than the one with which a program has been first compiled.
 */
int
XpmLibraryVersion()
{
    return XpmIncludeVersion;
}



void
xpmCreatePixmapFromImage(display, d, ximage, pixmap_return)
    Display *display;
    Drawable d;
    XImage *ximage;
    Pixmap *pixmap_return;
{
    GC gc;

    *pixmap_return = XCreatePixmap(display, d, ximage->width,
				   ximage->height, ximage->depth);
    gc = XCreateGC(display, *pixmap_return, 0, NULL);

    XPutImage(display, *pixmap_return, gc, ximage, 0, 0, 0, 0,
	      ximage->width, ximage->height);

    XFreeGC(display, gc);
}

void
xpmCreateImageFromPixmap(display, pixmap, ximage_return, width, height)
    Display *display;
    Pixmap pixmap;
    XImage **ximage_return;
    unsigned int *width;
    unsigned int *height;
{
    unsigned int dum;
    int dummy;
    Window win;

    if (*width == 0 && *height == 0)
	XGetGeometry(display, pixmap, &win, &dummy, &dummy,
		     width, height, &dum, &dum);

    *ximage_return = XGetImage(display, pixmap, 0, 0, *width, *height,
			       AllPlanes, ZPixmap);
}
