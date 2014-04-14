#include	<X11/Intrinsic.h>
#include	<X11/IntrinsicP.h>
#include	<X11/StringDefs.h>
#include	"../Xpm/xpm.h"


/*
 * XmuConvertStringToPixmap:
 *
 * creates a depth-n Pixmap suitable for window manager icons.
 * "string" represents a pixmap(1) filename which may be absolute,
 * or relative to the global resource pixmapFilePath, class
 * PixmapFilePath.  If the resource is not defined, the default
 * value is the build symbol PIXMAPDIR.
 *
 * chares lots of code with XmuConvertStringToCursor.  
 *
 * To use, include the following in your ClassInitialize procedure:

static XtConvertArgRec screenConvertArg[] = {
    {XtBaseOffset, (caddr_t) XtOffset(Widget, core.screen), sizeof(Screen *)}
};

    XtAddConverter("String", "Pixmap", XmuCvtStringToPixmap,
		   screenConvertArg, XtNumber(screenConvertArg));
 *
 */

#define	done(address, type) \
	{ (*toVal).size = sizeof(type); (*toVal).addr = (caddr_t) address; \
	    printf("got it\n");}


/*ARGSUSED*/
void 
XmuCvtStringToPixmap(XrmValuePtr args, Cardinal *num_args, 
	XrmValuePtr fromVal, XrmValuePtr toVal)
{
  static Pixmap pixmap;		/* static for cvt magic */
  static Pixmap shape_pixmap;
  Display      *dpy  = XDisplayOfScreen(*((Screen **) args[0].addr));
  Window        d    = XDefaultRootWindow(dpy);
  char         *name = (char *)fromVal->addr;
  int           status;
  char*         system_env;
  char*         path;
/*  char          suf[5] = ".xpm";
    char          sepr[2] = ":";
*/
  char          file_name[200];
  
  if (*num_args != 2)
    XtErrorMsg("wrongParameters","cvtStringToPixmap","XtToolkitError",
	       "String to pixmap conversion needs screen argument",
	       (String *)NULL, (Cardinal *)NULL);
  
  system_env = (char *) getenv("PIXMAP_PATH");
  XtStringConversionWarning ("PIXMAP_PATH -- ", system_env);
  printf("PIXMAP_PATH = %s\n",system_env);
  
  if (system_env != NULL ) {

    for(; (path = strtok(system_env, ":"))!=NULL; )  {
      file_name[0] = '\0';
      strcat(file_name, path);
      strcat(file_name, name);
      strcat(file_name, ".xpm");
      
      status = XpmReadFileToPixmap(dpy, d,
				   file_name, &pixmap, &shape_pixmap, NULL);

      if ((status != XpmOpenFailed) &&
	  (status != XpmFileInvalid)&&
	  (status != XpmNoMemory)	)
	break;
      
    }
  }else {
    
    file_name[0] = '\0';
    strcat(file_name, name);
    strcat(file_name, ".xpm");
    status = XpmReadFileToPixmap(dpy, d, 
				 file_name, &pixmap, &shape_pixmap, NULL);
  }
  
  if (pixmap != None) {
    done (&pixmap, Pixmap);
  } else {
    XtStringConversionWarning (name, "Pixmap");
    printf("failed\n");
    done (&pixmap, Pixmap);
  }
}

