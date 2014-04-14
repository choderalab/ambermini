#include "COPY.h"

/*
* SCCS_data: @(#) XtMacros.c 1.3 93/01/28 14:38:33
*
* Some old implementations of Xt provide some of the Xt functions only as
* macros, in violation of the Xt spec.  Here are implementations of the ones
* not provided by SCO Open Desk Top (Motif 1.0).
*/

#include <X11/IntrinsicP.h>
#include "WcCreateP.h"

#ifdef XtMapWidget
#undef XtMapWidget
#endif
#ifdef XtUnmapWidget
#undef XtUnmapWidget
#endif

void XtMapWidget(widget)
    Widget widget;
{
    XMapWindow(XtDisplay(widget), XtWindow(widget));
}

void XtUnmapWidget(widget)
    Widget widget;
{
    XUnmapWindow(XtDisplay(widget), XtWindow(widget));
}

#ifndef ISC301
char* XtName( w )
    Widget w;
{
    if (XtIsWidget(w))
	return w->core.name;
    else
	return XrmQuarkToString(w->core.xrm_name);
}
#endif /*ISC301*/
