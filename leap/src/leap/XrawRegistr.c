#include <X11/Xatom.h>
#include <X11/StringDefs.h>

#include <X11/IntrinsicP.h>
#include <X11/CoreP.h>
#include <X11/ObjectP.h>
#include <X11/RectObjP.h>
#include <X11/CompositeP.h>
#include <X11/ConstrainP.h>

#include "../Wc/WcCreateP.h"

#include "../Xraw/ArrowP.h"
#include "../Xraw/AsciiSinkP.h"
#include "../Xraw/AsciiSrcP.h"
#include "../Xraw/AsciiTextP.h"
#include "../Xraw/BoxP.h"
#include "../Xraw/ClockP.h"
#include "../Xraw/CommandP.h"
#include "../Xraw/DialogP.h"
#include "../Xraw/FormP.h"
#include "../Xraw/FrameP.h"
#include "../Xraw/GripP.h"
#include "../Xraw/LabelP.h"
#include "../Xraw/ListP.h"
#include "../Xraw/LogoP.h"
#include "../Xraw/MailboxP.h"
#include "../Xraw/MenuButtoP.h"
#include "../Xraw/PanedP.h"
#include "../Xraw/RepeaterP.h"
#include "../Xraw/SeparatorP.h"
#include "../Xraw/ScrollbarP.h"
#include "../Xraw/ScrolledTableP.h"
#include "../Xraw/SimpleMenP.h"
#include "../Xraw/SimpleP.h"
#include "../Xraw/SmeBSBP.h"
#include "../Xraw/SmeLineP.h"
#include "../Xraw/SmeP.h"
#include "../Xraw/StripCharP.h"
#include "../Xraw/TableP.h"
#include "../Xraw/TextP.h"
#include "../Xraw/TextSinkP.h"
#include "../Xraw/TextSrcP.h"
#include "../Xraw/ToggleP.h"
#include "../Xraw/ViewportP.h"

#include "XrawRegistr.h"

/* --  Register all Athena widgets, converters, callbacks, actions.
*******************************************************************************
*/

void XrawRegisterAll ( XtAppContext app )
{
  ONCE_PER_XtAppContext( app );

  /* -- Force Athena to load XmuNewCvtStringToWidget so
   *    WcCvtStringToWidget stays in effect when loaded by WcInitialize()
   */

  XtInitializeWidgetClass( toggleWidgetClass );
  XtInitializeWidgetClass( formWidgetClass );

#define RCP( name, class ) WcRegisterClassPtr  ( app, name, class );
#define RCO( name, func )  WcRegisterConstructor(app, name, func  )
  
    RCP("Arrow",			arrowWidgetClass);
    RCP("arrowWidgetClass",		arrowWidgetClass);
    RCP("AsciiSink",			asciiSinkObjectClass);
    RCP("asciiSinkObjectClass",		asciiSinkObjectClass);
    RCP("AsciiSrc",			asciiSrcObjectClass);
    RCP("asciiSrcObjectClass",		asciiSrcObjectClass);
    RCP("AsciiText",			asciiTextWidgetClass);
    RCP("asciiTextWidgetClass",		asciiTextWidgetClass);
    RCP("Box",				boxWidgetClass);
    RCP("boxWidgetClass",		boxWidgetClass);
    RCP("Clock",			clockWidgetClass);
    RCP("clockWidgetClass",		clockWidgetClass);
    RCP("Command",			commandWidgetClass);
    RCP("commandWidgetClass",		commandWidgetClass);
    RCP("Dialog",			dialogWidgetClass);
    RCP("dialogWidgetClass",		dialogWidgetClass);
    RCP("Form",				formWidgetClass);
    RCP("formWidgetClass",		formWidgetClass);
    RCP("Frame",			frameWidgetClass);
    RCP("frameWidgetClass",		frameWidgetClass);
    RCP("Grip",				gripWidgetClass);
    RCP("gripWidgetClass",		gripWidgetClass);
    RCP("Label",			labelWidgetClass);
    RCP("labelWidgetClass",		labelWidgetClass);
    RCP("List",				listWidgetClass);
    RCP("listWidgetClass",		listWidgetClass);
    RCP("Logo",				logoWidgetClass);
    RCP("logoWidgetClass",		logoWidgetClass);
    RCP("MenuButton",			menuButtonWidgetClass);
    RCP("menuButtonWidgetClass",	menuButtonWidgetClass);
    RCP("Paned",			panedWidgetClass);
    RCP("panedWidgetClass",		panedWidgetClass);
#if 0
    RCP("Panner",			pannerWidgetClass);
    RCP("pannerWidgetClass",		pannerWidgetClass);
    RCP("Porthole",			portholeWidgetClass);
    RCP("portholeWidgetClass",		portholeWidgetClass);
#endif
    RCP("Repeater",			repeaterWidgetClass);
    RCP("repeaterWidgetClass",		repeaterWidgetClass);
    RCP("Scrollbar",			scrollbarWidgetClass);
    RCP("scrollbarWidgetClass",		scrollbarWidgetClass);
    RCP("ScrolledTable",                scrolledTableWidgetClass);
    RCP("scrolledTableWidgetClass",     scrolledTableWidgetClass);
    RCP("Separator",			separatorWidgetClass);
    RCP("separatorWidgetClass",		separatorWidgetClass);
    RCP("Simple",			simpleWidgetClass);
    RCP("SimpleMenu",			simpleMenuWidgetClass);
    RCP("simpleMenuWidgetClass",	simpleMenuWidgetClass);
    RCP("simpleWidgetClass",		simpleWidgetClass);
    RCP("Sme",				smeObjectClass);
    RCP("SmeBSB",			smeBSBObjectClass);
    RCP("smeBSBObjectClass",		smeBSBObjectClass);
    RCP("SmeLine",			smeLineObjectClass);
    RCP("smeLineObjectClass",		smeLineObjectClass);
    RCP("smeObjectClass",		smeObjectClass);
    RCP("StripChart",			stripChartWidgetClass);
    RCP("stripChartWidgetClass",	stripChartWidgetClass);
    RCP("Table",                        tableWidgetClass);  
    RCP("Text",				textWidgetClass);
    RCP("TextSink",			textSinkObjectClass);
    RCP("textSinkObjectClass",		textSinkObjectClass);
    RCP("TextSrc",			textSrcObjectClass);
    RCP("textSrcObjectClass",		textSrcObjectClass);
    RCP("textWidgetClass",		textWidgetClass);
    RCP("Toggle",			toggleWidgetClass);
    RCP("toggleWidgetClass",		toggleWidgetClass);
    RCP("Viewport",			viewportWidgetClass);
    RCP("viewportWidgetClass",		viewportWidgetClass);

}

