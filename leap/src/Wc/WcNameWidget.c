#include "COPY.h"

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* Widget Creation Library - WcName.c
*
* Implements several name-to-widget and widget-to-name functions as well as
* other general purpose parsing functions.
*
*******************************************************************************
                       -- WcChildNameToWidget

    If we have a recent version of Xt with an XtNameToWidget() which does
    not barf on Gadgets, then we can just use XtNameToWidget().
    Otherwise, we use the WcChildNameToWidget function implemented in the file
    XtName.c which is simply a fixed version of XtNameToWidget.

*/

#include <X11/IntrinsicP.h>

#ifdef sun
#include <X11/ObjectP.h>	/* why don't they just use X from mit!?! */
#include <X11/RectObjP.h>
#endif

#include <X11/StringDefs.h>

#include "WcCreateP.h"



#ifndef XtSpecificationRelease

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

char* XtName( w )
    Widget w;
{
    if (XtIsWidget(w))
	return w->core.name;
    else
	return XrmQuarkToString(w->core.xrm_name);
}

static Boolean ClassIsSubclassOf(class, superclass)
    WidgetClass class, superclass;
{
    for (; class != NULL; class = class->core_class.superclass) {
	if (class == superclass) return True;
    }
    return False;
}

void XtGetConstraintResourceList(widget_class, resources, num_resources)
	WidgetClass widget_class;
	XtResourceList *resources;
	Cardinal *num_resources;
{
	int size;
	register int i, dest = 0;
	register XtResourceList *list, dlist;
	ConstraintWidgetClass class = (ConstraintWidgetClass)widget_class;

	if (   (class->core_class.class_inited &&
		!(class->core_class.class_inited & _XtConstraintBit)) /* DES */
	    || (!class->core_class.class_inited &&
		!ClassIsSubclassOf(widget_class, constraintWidgetClass))
	    || class->constraint_class.num_resources == 0) {

	    *resources = NULL;
	    *num_resources = 0;
	    return;
	}

	size = class->constraint_class.num_resources * sizeof(XtResource);
	*resources = (XtResourceList) XtMalloc((unsigned) size);

	if (!class->core_class.class_inited) {
	    /* Easy case */

	    bcopy((char *)class->constraint_class.resources,
		    (char *) *resources, size);
	    *num_resources = class->constraint_class.num_resources;
	    return;
	}

	/* Nope, it's the hard case */

	list = (XtResourceList *) class->constraint_class.resources;
	dlist = *resources;
	for (i = 0; i < class->constraint_class.num_resources; i++) {
	    if (list[i] != NULL) {
		dlist[dest].resource_name = (String)
			XrmQuarkToString((XrmQuark) list[i]->resource_name);
		dlist[dest].resource_class = (String) 
			XrmQuarkToString((XrmQuark) list[i]->resource_class);
		dlist[dest].resource_type = (String)
			XrmQuarkToString((XrmQuark) list[i]->resource_type);
		dlist[dest].resource_size = list[i]->resource_size;
		dlist[dest].resource_offset = -(list[i]->resource_offset + 1);
		dlist[dest].default_type = (String)
			XrmQuarkToString((XrmQuark) list[i]->default_type);
		dlist[dest].default_addr = list[i]->default_addr;
		dest++;
	    }
	}
	*num_resources = dest;
}

#ifndef ALLOCATE_LOCAL
#define ALLOCATE_LOCAL(size) XtMalloc((unsigned long)(size))
#define DEALLOCATE_LOCAL(ptr) XtFree((caddr_t)(ptr))
#endif /* ALLOCATE_LOCAL */

#define _XtAllocError		XtError

static Widget NameListToWidget();

typedef Widget (*NameMatchProc)();

static Widget MatchExactChildren(names, bindings, children, num,
	in_depth, out_depth, found_depth)
    XrmNameList     names;
    XrmBindingList  bindings;
    register WidgetList children;
    register int num;
    int in_depth, *out_depth, *found_depth;
{
    register Cardinal   i;
    register XrmName    name = *names;
    Widget w, result = NULL;
    int d, min = 10000;

    for (i = 0; i < num; i++) {
	if (name == children[i]->core.xrm_name) {
	    w = NameListToWidget(children[i], &names[1], &bindings[1],
		    in_depth+1, &d, found_depth);
	    if (w != NULL && d < min) {result = w; min = d;}
	}
    }
    *out_depth = min;
    return result;
}

static Widget MatchWildChildren(names, bindings, children, num,
	in_depth, out_depth, found_depth)
    XrmNameList     names;
    XrmBindingList  bindings;
    register WidgetList children;
    register int num;
    int in_depth, *out_depth, *found_depth;
{
    register Cardinal   i;
    Widget w, result = NULL;
    int d, min = 10000;

    for (i = 0; i < num; i++) {
	w = NameListToWidget(children[i], names, bindings,
		in_depth+1, &d, found_depth);
	if (w != NULL && d < min) {result = w; min = d;}
    }
    *out_depth = min;
    return result;
}

static Widget SearchChildren(root, names, bindings, matchproc,
			     in_depth, out_depth, found_depth)
     Widget root;
     XrmNameList     names;
     XrmBindingList  bindings;
     NameMatchProc matchproc;
     int in_depth, *out_depth, *found_depth;
{
  Widget w1, w2;
  int d1, d2;
  
  if (XtIsComposite(root)) {
    w1 = (*matchproc)(names, bindings,
		      ((CompositeWidget) root)->composite.children,
		      ((CompositeWidget) root)->composite.num_children,
		      in_depth, &d1, found_depth);
  } else d1 = 10000;
  w2 = (*matchproc)(names, bindings, root->core.popup_list,
		    root->core.num_popups, in_depth, &d2, found_depth);
  *out_depth = (d1 < d2 ? d1 : d2);
  return (d1 < d2 ? w1 : w2);
}

static Widget NameListToWidget(root, names, bindings,
			       in_depth, out_depth, found_depth)
     register Widget root;
     XrmNameList     names;
     XrmBindingList  bindings;
     int in_depth, *out_depth, *found_depth;
{
  Widget w1, w2;
  int d1, d2;
  
  if (in_depth >= *found_depth) {
    *out_depth = 10000;
    return NULL;
  }
  
  if (names[0] == NULLQUARK) {
    *out_depth = *found_depth = in_depth;
    return root;
  }
  
  if (! XtIsWidget(root)) {
    *out_depth = 10000;
    return NULL;
  }
  
  if (*bindings == XrmBindTightly) {
    return SearchChildren(root, names, bindings, MatchExactChildren,
			  in_depth, out_depth, found_depth);
    
  } else {	/* XrmBindLoosely */
    w1 = SearchChildren(root, names, bindings, MatchExactChildren,
			in_depth, &d1, found_depth);
    w2 = SearchChildren(root, names, bindings, MatchWildChildren,
			in_depth, &d2, found_depth);
    *out_depth = (d1 < d2 ? d1 : d2);
    return (d1 < d2 ? w1 : w2);
  }
}

Widget WcChildNameToWidget(root, name)
    Widget root;
     String name;
{
  XrmName *names;
  XrmBinding *bindings;
  int len, depth, found = 10000;
  Widget result;
  
  len = strlen(name);
  if (len == 0) return NULL;
  
  names = (XrmName *) ALLOCATE_LOCAL((unsigned) (len+1) * sizeof(XrmName));
  bindings = (XrmBinding *)
    ALLOCATE_LOCAL((unsigned) (len+1) * sizeof(XrmBinding));
  if (names == NULL || bindings == NULL) _XtAllocError(NULL);
  
  XrmStringToBindingQuarkList(name, bindings, names);
  if (names[0] == NULLQUARK) {
    DEALLOCATE_LOCAL((char *) bindings);
    DEALLOCATE_LOCAL((char *) names);
    return NULL;
  }
  
  result = NameListToWidget(root, names, bindings, 0, &depth, &found);
  
  DEALLOCATE_LOCAL((char *) bindings);
  DEALLOCATE_LOCAL((char *) names);
  return result;
}



#endif /* undefine XtSpecificationRelease */





#if XtSpecificationRelease == 4




#ifndef ALLOCATE_LOCAL
#define ALLOCATE_LOCAL(size) XtMalloc((unsigned long)(size))
#define DEALLOCATE_LOCAL(ptr) XtFree((caddr_t)(ptr))
#endif /* ALLOCATE_LOCAL */

#define _XtAllocError		XtError

static Widget NameListToWidget();

typedef Widget (*NameMatchProc)();

static Widget MatchExactChildren(names, bindings, children, num,
	in_depth, out_depth, found_depth)
    XrmNameList     names;
    XrmBindingList  bindings;
    register WidgetList children;
    register int num;
    int in_depth, *out_depth, *found_depth;
{
    register Cardinal   i;
    register XrmName    name = *names;
    Widget w, result = NULL;
    int d, min = 10000;

    for (i = 0; i < num; i++) {
	if (name == children[i]->core.xrm_name) {
	    w = NameListToWidget(children[i], &names[1], &bindings[1],
		    in_depth+1, &d, found_depth);
	    if (w != NULL && d < min) {result = w; min = d;}
	}
    }
    *out_depth = min;
    return result;
}

static Widget MatchWildChildren(names, bindings, children, num,
	in_depth, out_depth, found_depth)
    XrmNameList     names;
    XrmBindingList  bindings;
    register WidgetList children;
    register int num;
    int in_depth, *out_depth, *found_depth;
{
    register Cardinal   i;
    Widget w, result = NULL;
    int d, min = 10000;

    for (i = 0; i < num; i++) {
	w = NameListToWidget(children[i], names, bindings,
		in_depth+1, &d, found_depth);
	if (w != NULL && d < min) {result = w; min = d;}
    }
    *out_depth = min;
    return result;
}

static Widget SearchChildren(root, names, bindings, matchproc,
			     in_depth, out_depth, found_depth)
     Widget root;
     XrmNameList     names;
     XrmBindingList  bindings;
     NameMatchProc matchproc;
     int in_depth, *out_depth, *found_depth;
{
  Widget w1, w2;
  int d1, d2;
  
  if (XtIsComposite(root)) {
    w1 = (*matchproc)(names, bindings,
		      ((CompositeWidget) root)->composite.children,
		      ((CompositeWidget) root)->composite.num_children,
		      in_depth, &d1, found_depth);
  } else d1 = 10000;
  w2 = (*matchproc)(names, bindings, root->core.popup_list,
		    root->core.num_popups, in_depth, &d2, found_depth);
  *out_depth = (d1 < d2 ? d1 : d2);
  return (d1 < d2 ? w1 : w2);
}

static Widget NameListToWidget(root, names, bindings,
			       in_depth, out_depth, found_depth)
     register Widget root;
     XrmNameList     names;
     XrmBindingList  bindings;
     int in_depth, *out_depth, *found_depth;
{
  Widget w1, w2;
  int d1, d2;
  
  if (in_depth >= *found_depth) {
    *out_depth = 10000;
    return NULL;
  }
  
  if (names[0] == NULLQUARK) {
    *out_depth = *found_depth = in_depth;
    return root;
  }
  
  if (! XtIsWidget(root)) {
    *out_depth = 10000;
    return NULL;
  }
  
  if (*bindings == XrmBindTightly) {
    return SearchChildren(root, names, bindings, MatchExactChildren,
			  in_depth, out_depth, found_depth);
    
  } else {	/* XrmBindLoosely */
    w1 = SearchChildren(root, names, bindings, MatchExactChildren,
			in_depth, &d1, found_depth);
    w2 = SearchChildren(root, names, bindings, MatchWildChildren,
			in_depth, &d2, found_depth);
    *out_depth = (d1 < d2 ? d1 : d2);
    return (d1 < d2 ? w1 : w2);
  }
}

Widget WcChildNameToWidget(root, name)
    Widget root;
     String name;
{
  XrmName *names;
  XrmBinding *bindings;
  int len, depth, found = 10000;
  Widget result;
  
  len = strlen(name);
  if (len == 0) return NULL;
  
  names = (XrmName *) ALLOCATE_LOCAL((unsigned) (len+1) * sizeof(XrmName));
  bindings = (XrmBinding *)
    ALLOCATE_LOCAL((unsigned) (len+1) * sizeof(XrmBinding));
  if (names == NULL || bindings == NULL) _XtAllocError(NULL);
  
  XrmStringToBindingQuarkList(name, bindings, names);
  if (names[0] == NULLQUARK) {
    DEALLOCATE_LOCAL((char *) bindings);
    DEALLOCATE_LOCAL((char *) names);
    return NULL;
  }
  
  result = NameListToWidget(root, names, bindings, 0, &depth, &found);
  
  DEALLOCATE_LOCAL((char *) bindings);
  DEALLOCATE_LOCAL((char *) names);
  return result;
}



#else /* XtSpecificationRelease */



Widget WcChildNameToWidget( ref, childName )
     Widget ref;
     char*  childName;
{
  return XtNameToWidget( ref, childName );
}


#endif /* XtSpecificationRelease */

