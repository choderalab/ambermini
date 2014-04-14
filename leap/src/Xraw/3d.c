#include "3d.h"
#include "color.h"

#include "XrawDebug.h"


#define WHITE WhitePixelOfScreen
#define BLACK BlackPixelOfScreen
#define SWITCH(a,b) (top_or_bottom == TOP ? a : b)

#define DEPTH_SCREEN(w) DefaultDepthOfScreen(XtScreen(w))

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define TOP_VALUE_RATIO (1.3)
#define BOT_VALUE_RATIO (0.6)
#define ARM_VALUE_RATIO (0.85)

#define Top_Cash (1<<0L)
#define Bot_Cash (1<<1L)
#define Arm_Cash (1<<2L)

typedef struct _ColorCashRec {
  struct _ColorCashRec* nxt;
  XColor         col;
  RGB            rgb[1];
  HSV            hsv[1];
  RGB            top[1];
  RGB            bot[1];
  RGB            arm[1];
  unsigned int   set; 
}ColorCashRec, *ColorCash;

static ColorCashRec* color_cash = NULL;
static void GetTopShadow();


static void GetTopShadow(trom ,to)
     XColor* trom;
     XColor* to;
{
  ColorCashRec* cash;
  float         save;
  
  for (cash = color_cash; cash != NULL; cash = cash->nxt)
  {
    if (cash->col.red   == trom->red &&
	cash->col.green == trom->green &&
	cash->col.blue  == trom->blue)
    {
      if (cash->set & Top_Cash)
      {
	to->red   = cash->top[0].r;
	to->green = cash->top[0].g;
	to->blue  = cash->top[0].b;
      }
      else
      {
	save         = cash->hsv[0].v;
	cash->hsv[0].v *= TOP_VALUE_RATIO;
	cash->hsv[0].v  = MIN(cash->hsv[0].v, 1.0);
	
	HSVToRGB((HSV*)cash->hsv, (RGB*)cash->top);

	cash->hsv[0].v  = save;
	
	to->red   = cash->top[0].r;
	to->green = cash->top[0].g;
	to->blue  = cash->top[0].b;

	cash->set |= Top_Cash;
      }

      return;
    }
  }

  cash =  XtNew (ColorCashRec);

  cash->col.red   = cash->rgb[0].r = trom->red;
  cash->col.green = cash->rgb[0].g = trom->green;
  cash->col.blue  = cash->rgb[0].b = trom->blue;
  
  RGBToHSV((RGB*)cash->rgb, (HSV*)cash->hsv);

  save         = cash->hsv[0].v;
  cash->hsv[0].v *= TOP_VALUE_RATIO;
  cash->hsv[0].v  = MIN(cash->hsv[0].v, 1.0);

  HSVToRGB((HSV*)cash->hsv, (RGB*)cash->top);
  
  cash->hsv[0].v  = save;
	
  to->red   = cash->top[0].r;
  to->green = cash->top[0].g;
  to->blue  = cash->top[0].b;
  
  cash->set = Top_Cash;

  cash->nxt  = color_cash;
  color_cash = cash;
}

static void GetBotShadow(trom ,to)
     XColor* trom;
     XColor* to;
{
  ColorCashRec* cash;
  float         save;
  
  for (cash = color_cash; cash != NULL; cash = cash->nxt)
  {
    if (cash->col.red   == trom->red &&
	cash->col.green == trom->green &&
	cash->col.blue  == trom->blue)
    {
      if (cash->set & Bot_Cash)
      {
	to->red   = cash->bot[0].r;
	to->green = cash->bot[0].g;
	to->blue  = cash->bot[0].b;
      }
      else
      {
	save         = cash->hsv[0].v;
	cash->hsv[0].v *= BOT_VALUE_RATIO;
	
	HSVToRGB((HSV*)cash->hsv, (RGB*)cash->bot);

	cash->hsv[0].v  = save;
	
	to->red   = cash->bot[0].r;
	to->green = cash->bot[0].g;
	to->blue  = cash->bot[0].b;

	cash->set |= Bot_Cash;
      }

      return;
    }
  }

  cash =  XtNew (ColorCashRec);

  cash->col.red   = cash->rgb[0].r = trom->red;
  cash->col.green = cash->rgb[0].g = trom->green;
  cash->col.blue  = cash->rgb[0].b = trom->blue;
  
  RGBToHSV((RGB*)cash->rgb, (HSV*)cash->hsv);

  save         = cash->hsv[0].v;
  cash->hsv[0].v *= BOT_VALUE_RATIO;

  HSVToRGB((HSV*)cash->hsv, (RGB*)cash->bot);
  
  cash->hsv[0].v  = save;
	
  to->red   = cash->bot[0].r;
  to->green = cash->bot[0].g;
  to->blue  = cash->bot[0].b;
  
  cash->set = Bot_Cash;

  cash->nxt  = color_cash;
  color_cash = cash;
}

static void GetArmShadow(trom ,to)
     XColor* trom;
     XColor* to;
{
  ColorCashRec* cash;
  float         save;
  
  for (cash = color_cash; cash != NULL; cash = cash->nxt)
  {
    if (cash->col.red   == trom->red &&
	cash->col.green == trom->green &&
	cash->col.blue  == trom->blue)
    {
      if (cash->set & Arm_Cash)
      {
	to->red   = cash->arm[0].r;
	to->green = cash->arm[0].g;
	to->blue  = cash->arm[0].b;
      }
      else
      {
	save         = cash->hsv[0].v;
	cash->hsv[0].v *= ARM_VALUE_RATIO;
	
	HSVToRGB((HSV*)cash->hsv, (RGB*)cash->arm);

	cash->hsv[0].v  = save;
	
	to->red   = cash->arm[0].r;
	to->green = cash->arm[0].g;
	to->blue  = cash->arm[0].b;

	cash->set |= Arm_Cash;
      }

      return;
    }
  }

  cash =  XtNew (ColorCashRec);

  cash->col.red   = cash->rgb[0].r = trom->red;
  cash->col.green = cash->rgb[0].g = trom->green;
  cash->col.blue  = cash->rgb[0].b = trom->blue;
  
  RGBToHSV((RGB*)cash->rgb, (HSV*)cash->hsv);

  save         = cash->hsv[0].v;
  cash->hsv[0].v *= ARM_VALUE_RATIO;

  HSVToRGB((HSV*)cash->hsv, (RGB*)cash->arm);
  
  cash->hsv[0].v  = save;
	
  to->red   = cash->arm[0].r;
  to->green = cash->arm[0].g;
  to->blue  = cash->arm[0].b;
  
  cash->set = Arm_Cash;

  cash->nxt  = color_cash;
  color_cash = cash;
}

unsigned int shadowpm_width = 8;
unsigned int shadowpm_height= 8;

static char shadow_bits[] = {
    0xaa, 0x55, 0xaa, 0x55, 0xaa, 0x55, 0xaa, 0x55};

static char mtshadowpm_bits[] = {
    0x92, 0x24, 0x49, 0x92, 0x24, 0x49, 0x92, 0x24};

static char mbshadowpm_bits[] = {
    0x6d, 0xdb, 0xb6, 0x6d, 0xdb, 0xb6, 0x6d, 0xdb};


GC 
AllocGCFromPixmap (Widget w, Pixmap pixmap)
{
  XGCValues values;

  if (pixmap != None) {
    values.tile       = pixmap;
    values.fill_style = FillTiled;
    return XtGetGC(w, (XtGCMask)(GCTile | GCFillStyle), &values);
  }

  return XtGetGC(w, (XtGCMask)0, (XGCValues *)NULL);
}

GC 
AllocGCFromPixel (Widget w, Pixel fore)
{
  XGCValues values;
  
  values.foreground = fore;
  return XtGetGC(w, GCForeground, &values);
}

static Pixmap Depth_1_ShadowPixmap (w, top_or_bottom)
     Widget w;
     int top_or_bottom;
{
  Screen  *scn = XtScreen (w);

  if (DEPTH_SCREEN(w) == 1) 
    return XCreatePixmapFromBitmapData (XtDisplay (w),
					RootWindowOfScreen (scn),
					shadow_bits,
					shadowpm_width,
					shadowpm_height,
					SWITCH(BLACK (scn), WHITE (scn)),
					SWITCH(WHITE (scn), BLACK (scn)),
					1);
  else
    return None;
}

static Pixmap Depth_NOT_1_ShadowPixmap (w, colour, top_or_bottom)
     Widget w;
     Pixel colour;
     int top_or_bottom;
{
  Display *dpy = XtDisplay (w);
  Screen  *scn = XtScreen (w);
  unsigned long fore;
  unsigned long back;
  char         *pm_data;
  int           depth = DEPTH_SCREEN(w);
  
  if (depth == 1)
    return None;
  
  if (colour == WHITE (scn) || colour == BLACK (scn)) {
    fore    = WHITE (scn);
    back    = BLACK (scn);
    pm_data = SWITCH(mtshadowpm_bits, mbshadowpm_bits); 
  } else {
    fore    = colour;
    back    = SWITCH (WHITE (scn), BLACK (scn));
    pm_data = shadow_bits;
  }

  return XCreatePixmapFromBitmapData (dpy,
				      RootWindowOfScreen (scn),
				      pm_data,
				      shadowpm_width,
				      shadowpm_height,
				      fore,
				      back,
				      depth);
}

Pixmap 
CreateShadowPixmap (Widget w, Pixel colour, int top_or_bottom)
{
  if (DEPTH_SCREEN(w) == 1)
    return Depth_1_ShadowPixmap (w, top_or_bottom);
  else
    return Depth_NOT_1_ShadowPixmap (w, colour, top_or_bottom);
}

#define _MIN(x,y) (unsigned short) ((x) < (y)) ? (x) : (y)
#define _MAX(x,y) (unsigned short) ((x) < (y)) ? (y) : (x)


Boolean 
AllocShadowPixel (
     Widget w,
     Pixel  base,
     int    brightness,   
     Pixel *result)
{
  XColor         set;
  XColor         get;
  double         mult;
  Colormap       cmap= ((CoreWidget)w)->core.colormap;
  unsigned short red;
  unsigned short green;
  unsigned short blue; 
  
  get.pixel = base;
  XQueryColor (XtDisplay (w), cmap, &get);
  mult = (double)(100 + brightness) / ((double) 100.0);

  red   = mult * (double)get.red;
  green = mult * (double)get.green;
  blue  = mult * (double)get.blue;

  set.red   = _MAX(0,_MIN (65535, red));
  set.green = _MAX(0,_MIN (65535, green));
  set.blue  = _MAX(0,_MIN (65535, blue));

#define EQ(field) (set.field == get.field)  
  if (EQ(red) && EQ(green) && EQ(blue))
#undef EQ  
    return False;

  if (XAllocColor (XtDisplay (w), cmap, &set) != 0) {
    *result = set.pixel;
    return True;
  } else
    return False;
}

GC 
MakeGC (Widget w,
	Pixel base,
	int brightness,
	Boolean pseudo,
	int top_or_bottom)
{
  Pixel fore;
  Pixmap tile;

  if (DEPTH_SCREEN(w) > 1)  {
    if (pseudo)     {
      tile = Depth_NOT_1_ShadowPixmap (w, base, top_or_bottom);
      return AllocGCFromPixmap(w, tile);
    }     else     {
      if (AllocShadowPixel(w, base, brightness, &fore))      {
	return AllocGCFromPixel(w, fore);
      }  else       {
	tile = Depth_NOT_1_ShadowPixmap (w, base, top_or_bottom);
	return AllocGCFromPixmap(w, tile);
      }
    }
  }  else   {
    tile =  Depth_1_ShadowPixmap (w, top_or_bottom);
    return AllocGCFromPixmap(w, tile);
  }
}

GC 
MakeTopShadowGC (Widget w, Pixel base)
{
  Pixel fore;
  Pixmap tile;

  if (DEPTH_SCREEN(w) > 1)  {
    if (TopShadowColor(w, base, &fore))      {
      return AllocGCFromPixel(w, fore);
    }  else       {
      tile = Depth_NOT_1_ShadowPixmap (w, base, TOP);
      return AllocGCFromPixmap(w, tile);
    }
  }  else   {
    tile =  Depth_1_ShadowPixmap (w, TOP);
    return AllocGCFromPixmap(w, tile);
  }
}


GC 
MakeBottomShadowGC (Widget w, Pixel base)
{
  Pixel fore;
  Pixmap tile;

  if (DEPTH_SCREEN(w) > 1)  {
    if (BottomShadowColor(w, base, &fore))      {
      return AllocGCFromPixel(w, fore);
    }  else       {
      tile = Depth_NOT_1_ShadowPixmap (w, base, BOTTOM);
      return AllocGCFromPixmap(w, tile);
    }
  }  else   {
    tile =  Depth_1_ShadowPixmap (w, BOTTOM);
    return AllocGCFromPixmap(w, tile);
  }
}

GC 
MakeArmedGC (Widget w, Pixel base)
{
  Pixel fore;
  Pixmap tile;

  if (DEPTH_SCREEN(w) > 1)  {
    if (ArmedColor(w, base, &fore))      {
      return AllocGCFromPixel(w, fore);
    }  else       {
      tile = Depth_NOT_1_ShadowPixmap (w, base, BOTTOM);
      return AllocGCFromPixmap(w, tile);
    }
  }  else   {
    tile =  Depth_1_ShadowPixmap (w, BOTTOM);
    return AllocGCFromPixmap(w, tile);
  }
}

void 
XawDrawFrame (Widget       gw,
	      Position     x,
	      Position     y,
	      Dimension    w,
	      Dimension    h,
	      XawFrameType frame_type,
	      Dimension    t,
	      GC           lightgc,
	      GC           darkgc)
{
  XPoint top_polygon[6];
  XPoint bottom_polygon[6];
  XPoint points[3];

  
  if (t == 0 || w == 0 || h == 0)
    return;

  if( lightgc == (GC)NULL ){
    XtWarning("XawDrawFrame: lightgc is NULL in XawDrawFrame.");
    return;
  }
    
  if( darkgc == (GC)NULL ){
    XtWarning("XawDrawFrame: darkgc is NULL in XawDrawFrame.");
    return;
  }

  if (!XtIsRealized(gw)) {
    XtWarning("XawDrawFrame: widget is not realized!!!");
    return;
  }
  
#define topPolygon(i,xx,yy)           \
  top_polygon[i].x = (short) (xx);    \
  top_polygon[i].y = (short) (yy)

#define bottomPolygon(i,xx,yy)        \
  bottom_polygon[i].x = (short) (xx); \
  bottom_polygon[i].y = (short) (yy)
  

  if (frame_type == XawTACK && t <= 2)
    frame_type = XawLEDGED;

  switch (frame_type) {

  case XawRAISED :
  case XawSUNKEN :

    topPolygon (0,x    ,y    ); bottomPolygon (0,x+w  ,y+h  ); 
    topPolygon (1,x+w  ,y    ); bottomPolygon (1,x    ,y+h  );
    topPolygon (2,x+w-t,y+t  ); bottomPolygon (2,x+t  ,y+h-t);
    topPolygon (3,x+t  ,y+t  ); bottomPolygon (3,x+w-t,y+h-t);
    topPolygon (4,x+t  ,y+h-t); bottomPolygon (4,x+w-t,y+t  );
    topPolygon (5,x    ,y+h  ); bottomPolygon (5,x+w  ,y    );
    
    if (frame_type == XawSUNKEN) 
    {
      XFillPolygon(XtDisplayOfObject(gw), XtWindowOfObject(gw), darkgc,
		   top_polygon, 6, Nonconvex, CoordModeOrigin);
      
      XFillPolygon(XtDisplayOfObject(gw), XtWindowOfObject(gw), lightgc,
		   bottom_polygon, 6, Nonconvex, CoordModeOrigin);
    } 
    else if (frame_type == XawRAISED)
    {
      XFillPolygon(XtDisplayOfObject(gw), XtWindowOfObject(gw), lightgc,
		   top_polygon, 6, Nonconvex, CoordModeOrigin);
      XFillPolygon(XtDisplayOfObject(gw), XtWindowOfObject(gw), darkgc,
		   bottom_polygon, 6, Nonconvex, CoordModeOrigin);
    }
    
    break;

  case XawTACK :

    t -= 2;

    topPolygon (0,x    ,y    ); bottomPolygon (0,x+w  ,y+h  ); 
    topPolygon (1,x+w  ,y    ); bottomPolygon (1,x    ,y+h  );
    topPolygon (2,x+w-t,y+t  ); bottomPolygon (2,x+t  ,y+h-t);
    topPolygon (3,x+t  ,y+t  ); bottomPolygon (3,x+w-t,y+h-t);
    topPolygon (4,x+t  ,y+h-t); bottomPolygon (4,x+w-t,y+t  );
    topPolygon (5,x    ,y+h  ); bottomPolygon (5,x+w  ,y    );
    
    XFillPolygon(XtDisplayOfObject(gw), XtWindowOfObject(gw), lightgc,
		 top_polygon, 6, Nonconvex, CoordModeOrigin);
    XFillPolygon(XtDisplayOfObject(gw), XtWindowOfObject(gw), darkgc,
		 bottom_polygon, 6, Nonconvex, CoordModeOrigin);
    
    points[0].x = x + t + 1;	points[0].y = y + h - t - 2;
    points[1].x = x + t + 1;	points[1].y = y + t + 1;
    points[2].x = x + w - t - 2;	points[2].y = y + t + 1;
    
    XDrawLines (XtDisplayOfObject(gw), XtWindowOfObject(gw), darkgc,
		points, 3, CoordModeOrigin);
    
 /* points[0].x = x + t + 1;	        points[0].y = y + h -t - 1;   */
    points[1].x = x + w - t - 2;	points[1].y = y + h - t - 2;
 /* points[2].x = x + w - t - 1;	points[2].y = y + t + 1;      */
    
    XDrawLines (XtDisplayOfObject(gw), XtWindowOfObject(gw), lightgc,
		points, 3, CoordModeOrigin);

    break;

  case XawLEDGED :

    XawDrawFrame(gw, x, y, w, h, XawRAISED, t/2, lightgc, darkgc);
    XawDrawFrame(gw, (Position)(x + t/2), (Position)(y + t/2),
		 (Dimension)(w - 2 * (t/2)), (Dimension)(h - 2 * (t/2)),
		 XawSUNKEN, t/2, lightgc, darkgc);
    break;

  case XawCHISELED :

    XawDrawFrame(gw, x, y, w, h, XawSUNKEN, t/2, lightgc, darkgc);
    XawDrawFrame(gw, (Position)(x + t/2),(Position)(y + t/2),
		 (Dimension)(w - 2 * (t/2)), (Dimension)(h - 2 * (t/2)),
		 XawRAISED, t/2, lightgc, darkgc);
    break;

  default :
    break;

  }

#undef topPolygon
#undef bottomPolygon
}


Boolean
BottomShadowColor( Widget widget,
		   Pixel  base,
		   Pixel  *result)
{
  Colormap colormap;
  XColor color;
  
  if (XtIsWidget(widget))
    colormap = widget->core.colormap;
  else
    colormap = (XtParent(widget))->core.colormap;

  color.pixel = base;
  
  XQueryColor(XtDisplay(widget), colormap, &color);

  GetBotShadow(&color, &color);
  
  if (XAllocColor (XtDisplay(widget), colormap, &color) != 0)
  {
    *result = color.pixel;
    return True;
  }
  else
  {
    HSV hsv;
    RGB rgb;

    rgb.r = color.red;
    rgb.g = color.green;
    rgb.b = color.blue;

    RGBToHSV ((RGB*)&rgb, (HSV*)&hsv);

    if (hsv.v > 0.5)
      *result = WhitePixelOfScreen(XtScreen (widget));
    else
      *result = BlackPixelOfScreen(XtScreen (widget));
    return True;
  }

}

Boolean
TopShadowColor( Widget widget,
	       Pixel  base,
	       Pixel  *result)
{
  Colormap colormap;
  XColor color;
  
  if (XtIsWidget(widget))
    colormap = widget->core.colormap;
  else
    colormap = (XtParent(widget))->core.colormap;

  color.pixel = base;
  
  XQueryColor(XtDisplay(widget), colormap, &color);

  GetTopShadow(&color, &color);

  if (XAllocColor (XtDisplay(widget), colormap, &color) != 0)
  {
    *result = color.pixel;
    return True;
  }
  else
  {
    HSV hsv;
    RGB rgb;

    rgb.r = color.red;
    rgb.g = color.green;
    rgb.b = color.blue;

    RGBToHSV ((RGB*)&rgb, (HSV*)&hsv);

    if (hsv.v > 0.5)
      *result = WhitePixelOfScreen(XtScreen (widget));
    else
      *result = BlackPixelOfScreen(XtScreen (widget));
    return True;
  }

}

Boolean
ArmedColor( Widget widget,
	       Pixel  base,
	       Pixel  *result)
{
  Colormap colormap;
  XColor color;
  
  if (XtIsWidget(widget))
    colormap = widget->core.colormap;
  else
    colormap = (XtParent(widget))->core.colormap;

  color.pixel = base;
  
  XQueryColor(XtDisplay(widget), colormap, &color);

  GetArmShadow(&color, &color);

  if (XAllocColor (XtDisplay(widget), colormap, &color) != 0)
  {
    *result = color.pixel;
    return True;
  }
  else
  {
    HSV hsv;
    RGB rgb;

    rgb.r = color.red;
    rgb.g = color.green;
    rgb.b = color.blue;

    RGBToHSV ((RGB*)&rgb, (HSV*)&hsv);

    if (hsv.v > 0.5)
      *result = WhitePixelOfScreen(XtScreen (widget));
    else
      *result = BlackPixelOfScreen(XtScreen (widget));
    return True;
  }

}



#undef assign_max
#undef assign_min


void 
DrawRhombus (
	      Widget w,
	      short x,
	      short y,
	      short g,
	      short t,
	      GC top_shadow_GC,
	      GC foreground_gc,
	      GC bottom_shadow_GC,
	      Boolean state )
{
  XPoint top_shade[6];
  XPoint bot_shade[6];
  XPoint center[4];

#define topPolygon(i,a,b)    top_shade[i].x = a; top_shade[i].y = b
#define bottomPolygon(i,a,b) bot_shade[i].x = a; bot_shade[i].y = b
#define centerPolygon(i,a,b)    center[i].x = a;    center[i].y = b
  
  topPolygon(0, x-g  , y    );    bottomPolygon(0, x-g  , y    );
  topPolygon(1, x-g+t, y    );    bottomPolygon(1, x-g+t, y    );
  topPolygon(2, x    , y-g+t);    bottomPolygon(2, x    , y+g-t);
  topPolygon(3, x+g-t, y    );    bottomPolygon(3, x+g-t, y    );
  topPolygon(4, x+g  , y    );    bottomPolygon(4, x+g  , y    );
  topPolygon(5, x    , y-g  );    bottomPolygon(5, x    , y+g  );
  
  if (state)
  {
    if (foreground_gc)
    {
      centerPolygon(0, x-g+t, y    );
      centerPolygon(1, x    , y-g+t);
      centerPolygon(2, x+g-t, y    );
      centerPolygon(3, x    , y+g-t);
      
      XFillPolygon(XtDisplayOfObject(w), XtWindowOfObject(w), foreground_gc,
		   center,  XtNumber(center),
		   Convex, CoordModeOrigin);
    }

    if (bottom_shadow_GC)
      XFillPolygon(XtDisplayOfObject(w), XtWindowOfObject(w), bottom_shadow_GC,
		   top_shade, XtNumber(top_shade), Nonconvex, CoordModeOrigin);
    
    if (top_shadow_GC)
      XFillPolygon(XtDisplayOfObject(w), XtWindowOfObject(w), top_shadow_GC,
		   bot_shade, XtNumber(bot_shade), Nonconvex, CoordModeOrigin);
  }else{

    if (top_shadow_GC)
      XFillPolygon(XtDisplayOfObject(w), XtWindowOfObject(w), top_shadow_GC,
		   top_shade, XtNumber(top_shade),
		   Nonconvex, CoordModeOrigin);
    
    if (bottom_shadow_GC)
      XFillPolygon(XtDisplayOfObject(w), XtWindowOfObject(w), bottom_shadow_GC,
		   bot_shade, XtNumber(bot_shade),
		   Nonconvex, CoordModeOrigin);
  }    
  
#undef topPolygon
#undef bottomPolygon
#undef centerPolygon
  
}

Boolean
FetchPixel (Widget w, String name, Pixel* pixel)
{
  XrmValue source, dest;

  source.size = strlen(name)+1;
  source.addr = name;
  dest.size = sizeof(Pixel);
  dest.addr = (caddr_t) pixel;
  
  return XtConvertAndStore(w, XtRString, &source, XtRPixel, &dest);
}
