/*
 *      File: xTank.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              This file contains routines for maintaining a Widget
 *              called a TANK.
 *              A TANK is a window on the screen that displays the
 *              contents of a 3D space containing a UNIT.
 *
 *              TANKs do not use translation tables for mouse events
 *              because of the problem of actions not being passed
 *              client data.  There would be no way of knowing
 *              which TANK is responsible for an event without
 *              searching a list of Widgets/TANKs.
 *
 *              A word about coordinate systems:
 *
 *              World coordinates are the coordinate system of the atoms.
 *              View coordinates are the coordinate system of the viewer.
 *                      The viewer is at the origin and is looking down the
 *                      negative Z axis with up as the positive Y axis and
 *                      right is the positive X axis.
 *              Screen coordinates are 2D coordinates of the screen.
 *                      When points are converted from screen coordinates to
 *                      view coordinates, their Z value is the average of
 *                      the front and back clipping planes.
 *
 */

#define         DEBUG_TRANSFORM 1

#include        <X11/IntrinsicP.h>
#include        <X11/StringDefs.h>
#include        <X11/cursorfont.h>

#include        "basics.h"


#include        "vector.h"
#include        "threed.h"
#include        "classes.h"
#include        "varArray.h"
#include        "parser.h"
#include        "xAction.h"
#include        "xTank.h"




                /* Actions all defined in action.h */


XtActionsRec    SxtaaActions[] = {
        { "xtapActionButtonADown",        (VFUNCTION)xtapActionButtonADown },
        { "xtapActionButtonAUp",          (VFUNCTION)xtapActionButtonAUp },
        { "xtapActionButtonBDown",        (VFUNCTION)xtapActionButtonBDown },
        { "xtapActionButtonBUp",          (VFUNCTION)xtapActionButtonBUp },
        { "xtapActionButtonCDown",        (VFUNCTION)xtapActionButtonCDown },
        { "xtapActionButtonCUp",          (VFUNCTION)xtapActionButtonCUp },
        { "xtapActionButtonDDown",        (VFUNCTION)xtapActionButtonDDown },
        { "xtapActionButtonDUp",          (VFUNCTION)xtapActionButtonDUp },
        { "xtapActionButtonEDown",        (VFUNCTION)xtapActionButtonEDown },
        { "xtapActionButtonEUp",          (VFUNCTION)xtapActionButtonEUp },
        { "xtapActionCenterUnit",         (VFUNCTION)xtapActionCenterUnit },
        { "xtapActionSetCancelHit",       (VFUNCTION)xtapActionSetCancelHit },
        { "xtapActionMouseMotion",        (VFUNCTION)xtapActionMouseMotion }
};




/*
 ************************************************************
 *
 *      Static variables
 *
 */

static  BOOL    SbDefinedCursors = FALSE;

static  Cursor  ScDraw;
static  Cursor  ScErase;
static  Cursor  ScSelect;

#include        "drawSource.xbm"
#include        "drawMask.xbm"
#include        "eraseSource.xbm"
#include        "eraseMask.xbm"
#include        "selectSource.xbm"
#include        "selectMask.xbm"



/*
 *------------------------------------------------------------
 *------------------------------------------------------------
 *
 *      Stuff that pertains to the display, graphics contexts,
 *      color maps, foreground and background colors.
 */


                /* Currently only support 8 colors */

#define TOTALCOLORS     8
#define MAXPIXELS       TOTALCOLORS


#define PARTIALPLANES   2               /* 1 front 1 back */
#define PARTIALPIXELS   MAXPIXELS       /* Used for fullcolor */
#define PARTIALTOTALPIXELS      MAXPIXELS + 4


#define TANKPAGEFLIP            1
#define TANKOFFSCREEN           2


typedef struct  {
                BOOL    bDoPartialPageFlip;
                GC      gcaSlow[MAXPIXELS];
                GC      gcaSlowDashed[MAXPIXELS];
                GC      gcaFastFront[MAXPIXELS];
                GC      gcaFastBack[MAXPIXELS];
                GC      gcPreparePageFlip;
                int     iXColorEntries;
                XColor  xcaFrontMap[2*(1<<PARTIALPLANES)];
                XColor  xcaBackMap[2*(1<<PARTIALPLANES)];
                } PAGEFLIPt;

typedef struct  {
                GC              gcaSlow[MAXPIXELS];
                GC              gcaSlowDashed[MAXPIXELS];
                GC              gcaFast[MAXPIXELS];
                } OFFSCREENt;

static  BOOL            SbInitializedClass = FALSE;
static  Colormap        ScMap;
static  BOOL            SbMonochrome;
static  int             SiBackgroundColor;
static  int             SiForegroundColor;
static  int             SiSelectedColor;
static  String          SsaRgb[MAXPIXELS];
static  String          SsMonoForeground;
static  String          SsMonoBackground;
static  int             SiaElementColors[MAXELEMENTS];

static  Display         *SdPDisplay;
static  Screen          *SsPScreen;
static  Window          SwRoot;

static  int             SiRenderWay;

                /* Depending on the value of SiRenderWay, one of the */
                /* following structures is used */

static  PAGEFLIPt       SpfPageFlip;
static  OFFSCREENt      SosOffScreen;

static  GC              SgcRubberBand;


                        

/*
 *--------------------------------------------------------------
 *--------------------------------------------------------------
 *--------------------------------------------------------------
 *
 *      Private routines.
 */




/*
 *      TankClassDefineCursors
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Define the shapes of the Cursors that are used by the
 *      TANK to indicate the state of the TANK.
 */
static void
TankClassDefineCursors()
{
Pixmap          pSource, pMask;
Drawable        dRoot;
Colormap        cMap;
XColor          xcWhite, xcBlack, xcTemp;
unsigned int    iX, iY;
BOOL            bGotSize;

    SbDefinedCursors = TRUE;

    dRoot = RootWindowOfScreen(SsPScreen);
    cMap = DefaultColormapOfScreen(SsPScreen);
    XAllocNamedColor( SdPDisplay, cMap, "White", &xcWhite, &xcTemp );
    XAllocNamedColor( SdPDisplay, cMap, "Black", &xcBlack, &xcTemp );

    bGotSize = XQueryBestCursor( SdPDisplay, dRoot,
                        drawSource_width, drawSource_height,
                        &iX, &iY );
    MESSAGE(( "Queried best cursor = %s\n", sBOOL(bGotSize) ));
    MESSAGE(( "Best cursor size: %d, %d\n", iX, iY ));

                /* Swap foreground/background pixel values if the black */
                /* pixel has bits set */

    MESSAGE(( "Black pixel value: %lX\n", BlackPixelOfScreen(SsPScreen) ));
    MESSAGE(( "White pixel value: %lX\n", WhitePixelOfScreen(SsPScreen) ));

#define PIXMAP( sData, iWidth, iHeight, fg, bg )        \
        XCreateBitmapFromData( SdPDisplay, dRoot, sData,\
                iWidth, iHeight )


        /* Define the TANKDRAW state cursor */

    pSource = PIXMAP( drawSource_bits, drawSource_width, drawSource_height,
                        xcWhite.pixel, xcBlack.pixel );
    pMask   = PIXMAP( drawMask_bits, drawMask_width, drawMask_height,
                        0, 1 );
    ScDraw = XCreatePixmapCursor( SdPDisplay, pSource, pMask,
                                        &xcWhite, &xcBlack,
                                        drawMask_x_hot,
                                        drawMask_y_hot );

        /* Define the TANKERASE state cursor */

    pSource = PIXMAP( eraseSource_bits, eraseSource_width, eraseSource_height,
                        xcWhite.pixel, xcBlack.pixel );
    pMask   = PIXMAP( eraseMask_bits, eraseMask_width, eraseMask_height,
                        0, 1 );
    ScErase = XCreatePixmapCursor( SdPDisplay, pSource, pMask,
                                        &xcWhite, &xcBlack,
                                        eraseMask_x_hot,
                                        eraseMask_y_hot );



        /* Define the TANKSELECT state cursor */

    pSource = PIXMAP( selectSource_bits, selectSource_width, 
                        selectSource_height, xcWhite.pixel, xcBlack.pixel );
    pMask   = PIXMAP( selectMask_bits, selectMask_width, selectMask_height,
                        0, 1 );
    ScSelect = XCreatePixmapCursor( SdPDisplay, pSource, pMask,
                                        &xcWhite, &xcBlack,
                                        selectMask_x_hot,
                                        selectMask_y_hot );

}



/*
 *      zcTankCursor
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the cursor for the current TANK state.
 */
static Cursor
zcTankCursor( int iNewState )
{
Cursor          cCursor;

    switch ( iNewState ) {
        case TANKDRAW:
            MESSAGE(( "Changing state to DRAW\n" ));
            cCursor = ScDraw;
            break;
        case TANKERASE:
            MESSAGE(( "Changing state to ERASE\n" ));
            cCursor = ScErase;
            break;
        case TANKSELECT:
            MESSAGE(( "Changing state to SELECT\n" ));
            cCursor = ScSelect;
            break;
        case TANKDRAGROTATE:
            MESSAGE(( "Changing state to DRAG\n" ));
            cCursor = XCreateFontCursor( SdPDisplay, XC_fleur );
            break;
        case TANKTWIST:
            MESSAGE(( "Changing state to TWIST\n" ));
            cCursor = XCreateFontCursor( SdPDisplay, XC_exchange );
            break;
        default:
            DFATAL(( "Illegal state!" ));
            break;
    }
    return(cCursor);
}




/*
 *------------------------------------------------------------
 *
 *      Initialize colormaps for the various rendering
 *      techniques.
 *
 */



/*
 *      TankClassInitElementColors
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Initialize the colors of the atoms.
 */
#define SEC( e, c ) (SiaElementColors[e] = c)
static void
TankClassInitElementColors()
{
int     i;

    SiBackgroundColor = TBLACK;
    SiForegroundColor = TWHITE;
    SiSelectedColor   = TPURPLE;

    /*
     *  set everything to default element color
     */
    for (i=0; i<MAXELEMENTS; i++)
        SEC( i, TGREEN );

    /*
     *  override default with some specific colors
     */

    SEC( LONEPAIR,      TYELLOW );
    SEC( HYDROGEN,      TWHITE );
    SEC( HELIUM,        TYELLOW );
    SEC( LITHIUM,       TYELLOW );
    SEC( BERILIUM,      TYELLOW );
    SEC( BORON,         TWHITE );
    SEC( CARBON,        TGREEN );
    SEC( NITROGEN,      TCYAN );
    SEC( OXYGEN,        TRED );
    SEC( FLOURINE,      TGREEN );
    SEC( NEON,          TYELLOW );
    SEC( SODIUM,        TYELLOW );
    SEC( MAGNESIUM,     TYELLOW );
    SEC( ALUMINUM,      TYELLOW );
    SEC( SILICON,       TYELLOW );
    SEC( PHOSPHORUS,    TYELLOW );
    SEC( SULFUR,        TYELLOW );
    SEC( CHLORINE,      TGREEN );
    SEC( ARGON,         TYELLOW );
}


/*
 *========================================================
 *========================================================
 *
 *      OffScreen buffer rendering.
 */


/*
 *      TankClassInitOffScreenColors
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This routine sets up the pixels for OFFSCREEN display
 *      of the UNIT.
 */
static void
TankClassInitOffScreenColors()
{
XColor          xcColor, xcExact;
XGCValues       xgcvNew;
GC              gcNew;
int             i;


    if ( SiForegroundColor == SiSelectedColor ) {
            DFATAL(( "Foreground color cannot be selected color" ));
    }

                /* First initialize the GC's for the offscreen buffer, */
                /* used in FAST drawing */

    if ( SbMonochrome ) {
    
    
                /* setup the fast drawing GC's */

                        /* Create a GC for monochrome foreground rendering */
                        
        XAllocNamedColor( SdPDisplay, ScMap, 
                SsMonoForeground, &xcColor, &xcExact );
        xgcvNew.foreground = xcColor.pixel;
        gcNew = XCreateGC( SdPDisplay, SwRoot, GCForeground, &xgcvNew );

        for ( i=0; i<PARTIALPIXELS; i++ ) SosOffScreen.gcaFast[i] = gcNew;

                        /* Create GC for monochrome background rendering */

        XAllocNamedColor( SdPDisplay, ScMap, 
                SsMonoBackground, &xcColor, &xcExact );
        xgcvNew.foreground = xcColor.pixel;
        SosOffScreen.gcaFast[SiBackgroundColor] =
                XCreateGC( SdPDisplay, SwRoot, GCForeground, &xgcvNew );


                /* Setup the slow drawing GC's */
        
                        /* Create solid GC's */

        for ( i=0; i<PARTIALPIXELS; i++ ) 
                SosOffScreen.gcaSlow[i] = SosOffScreen.gcaFast[i];

                                /* Create the selection GC */
        XGetGCValues( SdPDisplay, SosOffScreen.gcaSlow[SiSelectedColor], 
                        GCForeground, &xgcvNew );
        xgcvNew.line_width = DEFAULTSELECTTHICKNESS;
        SosOffScreen.gcaSlow[SiSelectedColor] = XCreateGC( SdPDisplay, 
                                                SwRoot,
                                                GCForeground|GCLineWidth,
                                                &xgcvNew );

                        /* Create Dashed GC's */
                        
        for ( i=0; i<PARTIALPIXELS; i++ ) {
            XGetGCValues( SdPDisplay, 
                SosOffScreen.gcaSlow[i], GCForeground, &xgcvNew );
            xgcvNew.line_style = LineOnOffDash;
            SosOffScreen.gcaSlowDashed[i] = XCreateGC( SdPDisplay, 
                                                SwRoot,
                                                GCForeground|GCLineStyle,
                                                &xgcvNew );
        }

    } else {
    
                /* Setup the slow drawing GC's */
        
                        /* Create solid GC's */

        for ( i=0; i<PARTIALPIXELS; i++ ) {
            XAllocNamedColor( SdPDisplay, ScMap, 
                        SsaRgb[i], &xcColor, &xcExact );
            xgcvNew.foreground = xcColor.pixel;
            SosOffScreen.gcaSlow[i] = 
                XCreateGC( SdPDisplay, SwRoot, GCForeground, &xgcvNew );
        }


                        /* Create Dashed GC's */
                        
        for ( i=0; i<PARTIALPIXELS; i++ ) {
            XGetGCValues( SdPDisplay, SosOffScreen.gcaSlow[i], GCForeground, &xgcvNew );
            xgcvNew.line_style = LineOnOffDash;
            SosOffScreen.gcaSlowDashed[i] = XCreateGC( SdPDisplay, SwRoot,
                                                        GCForeground|GCLineStyle,
                                                        &xgcvNew );
        }
    
                /* setup the fast drawing GC's */

        for ( i=0; i<PARTIALPIXELS; i++ ) {
            SosOffScreen.gcaFast[i] = SosOffScreen.gcaSlow[SiForegroundColor];
        }
        SosOffScreen.gcaFast[SiBackgroundColor] = SosOffScreen.gcaSlow[SiBackgroundColor];

    }


                /* Create the GC for rubber banding, use pixel that is the */
                /* logical OR of the two foreground pixels and use a plane */
                /* mask that is the logical OR of the two plane masks. */
                /* This will garantee that the rubber band will be visible */

    xgcvNew.foreground = (~0L);
    xgcvNew.function = GXxor;
    SgcRubberBand = XCreateGC( SdPDisplay, 
                                SwRoot,
                                GCForeground | GCFunction, &xgcvNew ); 


    MESSAGE(( "Pixel# for background=%d\n", SiBackgroundColor ));
    MESSAGE(( "Pixel# for foreground=%d\n", SiForegroundColor ));
    MESSAGE(( "Pixel# for selection =%d\n", SiSelectedColor ));
}      


/*
 *      bTankClassInitOffScreen
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Initialize the TANK for OFFSCREEN drawing.
 *      This sets up the TANK to be displayed in monochrome
 *      using a one bit deep OFFSCREEN buffer.
 */
static BOOL
bTankClassInitOffScreen()
{


    SiRenderWay = TANKOFFSCREEN;

                /* Check how deep the default colormap is */
                /* if it is only one bit deep then use monochrome colors */

    if ( DefaultDepthOfScreen(SsPScreen) == 1 ) {
        SbMonochrome = TRUE;
    } else {
        SbMonochrome = FALSE;
    }
    
    MESSAGE((  "Off screen buffer with color slow draw = %s\n",
        sBOOL(SbMonochrome) ));

    TankClassInitOffScreenColors();
    return(TRUE);
}



/*
 *===================================================================
 *===================================================================
 *
 *      Page flip graphics initialization
 *
 */

/*
 *      TankClassInitPartialPageFlip
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Setup the colormaps, and GC arrays for partial page
 *      flipping.
 *      Currently this is hard wired to support only one bit deep
 *      pages.
 */
static void
TankClassInitPartialPageFlip( unsigned long ulaPlaneMasks[PARTIALPIXELS], 
        unsigned long ulaPixels[PARTIALPIXELS] )
{
unsigned long   ulFront, ulBack;
int             i;
XGCValues       xgcvTemplate;
unsigned long   ulValues;
XColor          xcTemp;
XColor          xcaFullColors[PARTIALPIXELS];

#define MAKEMAPENTRY( array, i, pix, color ) {\
        array[i].pixel = pix;\
        XParseColor( SdPDisplay, \
                        ScMap, \
                        SsaRgb[color], \
                        &xcTemp ); \
        array[i].red   = xcTemp.red;\
        array[i].green = xcTemp.green;\
        array[i].blue  = xcTemp.blue;\
        array[i].flags = DoRed | DoGreen | DoBlue;}


    SpfPageFlip.bDoPartialPageFlip = TRUE;

        /* First create the XColor array used to modify the colormap */

    SpfPageFlip.iXColorEntries = 4;
    


MESSAGE(( "Back plane background pixel = %lX\n", ulaPixels[0] ));
MESSAGE(( "Back plane foreground pixel = %lX\n", 
                ulaPixels[0] | ulaPlaneMasks[1] ));
MESSAGE(( "Front plane background pixel = %lX\n", ulaPixels[0] ));
MESSAGE(( "Front plane foreground pixel = %lX\n", 
                ulaPixels[0] | ulaPlaneMasks[0] ));

                /* Do the front plane color map */
    i = 0;
    MAKEMAPENTRY( SpfPageFlip.xcaFrontMap, i, 
                  ulaPixels[0], 
                  SiBackgroundColor );
    i++;
    MAKEMAPENTRY( SpfPageFlip.xcaFrontMap, i, 
                  ulaPixels[0] | ulaPlaneMasks[0], 
                  SiForegroundColor );
    i++;
    MAKEMAPENTRY( SpfPageFlip.xcaFrontMap, i, 
                  ulaPixels[0] | ulaPlaneMasks[1], 
                  SiBackgroundColor );
    i++;
    MAKEMAPENTRY( SpfPageFlip.xcaFrontMap, i, 
                  ulaPixels[0] | ulaPlaneMasks[0] | ulaPlaneMasks[1], 
                  SiForegroundColor );
    i++;

                /* Do the back plane color map */
    i = 0;
    MAKEMAPENTRY( SpfPageFlip.xcaBackMap, i, 
                  ulaPixels[0], 
                  SiBackgroundColor );
    i++;
    MAKEMAPENTRY( SpfPageFlip.xcaBackMap, i, 
                  ulaPixels[0] | ulaPlaneMasks[0], 
                  SiBackgroundColor );
    i++;
    MAKEMAPENTRY( SpfPageFlip.xcaBackMap, i, 
                  ulaPixels[0] | ulaPlaneMasks[1], 
                  SiForegroundColor );
    i++;
    MAKEMAPENTRY( SpfPageFlip.xcaBackMap, i, 
                  ulaPixels[0] | ulaPlaneMasks[0] | ulaPlaneMasks[1], 
                  SiForegroundColor );
    i++;

        /* Then create two GCs for each plane, one for the background */
        /* color and one for the foreground */
        /* The background GCs will be the same */

    ulFront = ulaPixels[0] | ulaPlaneMasks[0];
    ulBack  = ulaPixels[0] | ulaPlaneMasks[1];

    ulValues = GCFunction | GCForeground | GCBackground | GCPlaneMask;
    xgcvTemplate.function   = GXcopy;
    xgcvTemplate.background = ulaPixels[0];

                /* Set the plane_mask for the foreground so that drawing */
                /* will only effect this plane */
    xgcvTemplate.plane_mask = ulaPlaneMasks[0];
    xgcvTemplate.foreground = ulaPixels[0];
    SpfPageFlip.gcaFastFront[SiBackgroundColor] =
        XCreateGC( SdPDisplay, SwRoot,
                        ulValues, &xgcvTemplate );
    xgcvTemplate.foreground = ulFront;
    SpfPageFlip.gcaFastFront[SiForegroundColor] =
        XCreateGC( SdPDisplay, SwRoot,
                        ulValues, &xgcvTemplate );

                /* Set the plane_mask for the foreground so that drawing */
                /* will only effect this plane */
    xgcvTemplate.plane_mask = ulaPlaneMasks[1];
    xgcvTemplate.foreground = ulaPixels[0];
    SpfPageFlip.gcaFastBack[SiBackgroundColor] =
        XCreateGC( SdPDisplay, SwRoot,
                        ulValues, &xgcvTemplate );
    xgcvTemplate.foreground = ulBack;
    SpfPageFlip.gcaFastBack[SiForegroundColor] =
        XCreateGC( SdPDisplay, SwRoot,
                        ulValues, &xgcvTemplate );

        /* Create a GC that is used to prepare the TANK for page flipping */
        /* It fills the TANK with the base pixel value for the two planes */
        /* It does not have a restricted plane_mask */

    xgcvTemplate.foreground = ulaPixels[0];
    SpfPageFlip.gcPreparePageFlip =
        XCreateGC( SdPDisplay, SwRoot,
                        GCForeground, &xgcvTemplate );
    
        /* Create the GC for rubber banding, use pixel that is the */
        /* logical OR of the two foreground pixels and use a plane */
        /* mask that is the logical OR of the two plane masks. */
        /* This will garantee that the rubber band will be visible */

    xgcvTemplate.foreground = ulaPixels[0] | ulaPixels[1];
    xgcvTemplate.plane_mask = ulaPixels[0] | ulaPixels[1];
    xgcvTemplate.function = GXxor;
    SgcRubberBand = XCreateGC( SdPDisplay, 
                                   SwRoot,
                                   GCPlaneMask | GCForeground | GCFunction, 
                                   &xgcvTemplate );

        /* Now create a full set of GCs for the full color displays */

    for ( i=0; i<PARTIALPIXELS; i++ ) {
        MAKEMAPENTRY( xcaFullColors, i, ulaPixels[i], i );
        xgcvTemplate.foreground = ulaPixels[i];
        SpfPageFlip.gcaSlow[i] =
                XCreateGC( SdPDisplay, 
                                   SwRoot,
                                   GCForeground,
                                   &xgcvTemplate );
    }

        /* Now create a full set of GCs for the full color displays */
        /* and dashed lines */

    for ( i=0; i<PARTIALPIXELS; i++ ) {
        xgcvTemplate.foreground = ulaPixels[i];
        xgcvTemplate.line_style = LineOnOffDash;
        SpfPageFlip.gcaSlowDashed[i] =
                XCreateGC( SdPDisplay, 
                                   SwRoot,
                                   GCForeground | GCLineStyle,
                                   &xgcvTemplate );
    }

        /* Define the colormap for the full color displays */

    XStoreColors( SdPDisplay, ScMap,
                     xcaFullColors, PARTIALPIXELS );

}


   



/*
 *      bTankClassInitPageFlip
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Setup the colors for page flipping.
 *      Determine what kind of page flipping can be done.
 *      If not enough bit planes can be allocated for full color
 *      in both the front and back planes then just allocate
 *      one plane each, and use monochrome for fast redraws and
 *      color for slow redraws.
 *      Return FALSE only if two bit planes could not be allocated.
 *      If two bit planes could not be allocated then we
 *      have to fall back to another graphics method.
 */
static BOOL
bTankClassInitPageFlip()
{
BOOL            bSuccess;
unsigned long   ulaPlaneMasks[PARTIALPLANES];
unsigned long   ulaPixels[PARTIALPIXELS];
XVisualInfo     xviNeed;
XVisualInfo     *xviPList;
int             iGot;


    SiRenderWay = TANKPAGEFLIP;

                /* Look up what kinds of visuals are supported */

    bSuccess = FALSE;
    iGot = 0;

                /* I can't get the colormaps too work, */
                /* for now parasitize the */
                /* default colormap */

    ScMap = DefaultColormapOfScreen(SsPScreen);
    xviNeed.class = PseudoColor;
    xviNeed.depth = 8;
    xviPList = XGetVisualInfo( SdPDisplay, 
                                VisualDepthMask | VisualClassMask,
                                &xviNeed, &iGot );

    if ( iGot != 0 ) {
        bSuccess = XAllocColorCells( SdPDisplay,
                                ScMap, FALSE,
                                ulaPlaneMasks, PARTIALPLANES,
                                ulaPixels, PARTIALPIXELS );
    } else {
        bSuccess = FALSE;
    }

    if ( iGot != 0 ) {
        XFree(xviPList);
    }
    if ( !bSuccess ) {
        MESSAGE(( "Could not implement partial page flipping\n" ));
    } else {
        MESSAGE(( "Implementing partial page flipping\n" ));
        TankClassInitPartialPageFlip( ulaPlaneMasks, ulaPixels );
    }
    return(bSuccess);
}




/* 
 *===================================================================
 *===================================================================
 *
 *      Initialise the TANK class.
 *
 *      Read Xdefaults for colors, and rendering method.
 *
 *      Initialise the rendering method specified using the
 *      'graphicMode' resource which can have the values:
 *              'determine', 'pageflip', 'offscreen'.
 */


/*
 *      TankClassInitialize
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Try to initialize the different drawing techniques
 *      in order of best technique to worst.
 *      The different tests also determine whether or not
 *      monochrome colors should be used for the NormalColors.
 */
static void
TankClassInitialize( Display *dPDisplay, Screen *sPScreen )
{
char    *cPMode, *cPTemp;

#define XGETDEF( sDest, sName, sDefault ) {\
    cPTemp = XGetDefault( dPDisplay, GsProgramName, sName );\
    if ( cPTemp == NULL ) sDest = sDefault;\
    else                  sDest = cPTemp;}

                /* Set up static variables */


                        /* Find out what graphics mode to use */
                        /* Can be one of: pageflip, offscreen, */
                        /* (anything else ) */

    XGETDEF( cPMode, "graphicMode", "determine" );
    XGETDEF( SsaRgb[0], "color0", "Black" );
    XGETDEF( SsaRgb[1], "color1", "Red" );
    XGETDEF( SsaRgb[2], "color2", "Green" );
    XGETDEF( SsaRgb[3], "color3", "Blue" );
    XGETDEF( SsaRgb[4], "color4", "magenta" );
    XGETDEF( SsaRgb[5], "color5", "cyan" );
    XGETDEF( SsaRgb[6], "color6", "yellow" );
    XGETDEF( SsaRgb[7], "color7", "White" );
    XGETDEF( SsMonoForeground, "monoForeground", "Black" );
    XGETDEF( SsMonoBackground, "monoBackground", "White" );

    SdPDisplay = dPDisplay;
    SsPScreen = sPScreen;
    SwRoot = RootWindowOfScreen(SsPScreen);

                /* Set up the cursors */

    TankClassDefineCursors();

                /* Initialize the drawing method */

    SbMonochrome = FALSE;

    TankClassInitElementColors();


                /* Now determine the rendering method */
                
    if ( strcmp( cPMode, "determine" ) == 0 ) {
        if ( !bTankClassInitPageFlip() ) 
            if ( !bTankClassInitOffScreen() ) 
                DFATAL(( "The display will not support any graphics!\n" ));
        return;
                
    } else if ( strcmp( cPMode, "pageflip" )==0 ) {
        if ( bTankClassInitPageFlip() ) return;
        DFATAL(( "Cannot initialise pageflip graphicMode" ));

    } else if ( strcmp( cPMode, "offscreen" )==0 ) {
        if ( bTankClassInitOffScreen() ) return;
        DFATAL(( "Cannot initialise offscreen graphicMode" ));
    }
    DFATAL(( 
     "Unknown graphicMode type:  'determine', 'pageflip', or 'offscreen'" ));

}






/*
 *------------------------------------------------------------------
 *
 *      Get graphics contexts for the different rendering
 *      methods.
 */


/*
 *      gcTankColorIndex
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the gc for the color represented by iColor.
 */
static GC
gcTankColorIndex( TANK tTank, int iColor, BOOL bFast, BOOL bSolidLine )
{
GC              gc;

    switch ( SiRenderWay ) {
        case TANKOFFSCREEN:
            if ( bFast ) {
                gc = SosOffScreen.gcaFast[iColor];
            } else {
                if ( bSolidLine ) gc = SosOffScreen.gcaSlow[iColor];
                else gc = SosOffScreen.gcaSlowDashed[iColor];
            }
            break;
        case TANKPAGEFLIP:
            if ( SpfPageFlip.bDoPartialPageFlip && bFast ) {
                if ( iColor != SiBackgroundColor )
                    iColor = SiForegroundColor;
            }
            if ( !SpfPageFlip.bDoPartialPageFlip ||
                 (SpfPageFlip.bDoPartialPageFlip && !bFast) ) {
                if ( bSolidLine ) {
                    gc = SpfPageFlip.gcaSlow[iColor];
                } else {
                    gc = SpfPageFlip.gcaSlowDashed[iColor];
                }
            } else if ( tTank->tank.bPageFlipFront ) {
                gc = SpfPageFlip.gcaFastFront[iColor];
            } else {
                gc = SpfPageFlip.gcaFastBack[iColor];
            }
            break;
    }
    return(gc);
}





/*
 *      gcTankColor
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the the graphics context representing the color for the atom.
 *      If the atom is NULL then return the background color.
 */
static GC
gcTankColor( TANK tTank, ATOM aAtom, BOOL bFast, BOOL bSolid )
{
int             iColor, iElement;
GC              gcNew;

    if ( aAtom == NULL ) {
        iColor = SiBackgroundColor;
    } else if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) { 
        iColor = SiSelectedColor;
    } else if ( SbMonochrome ) {
        /* If monochrome then limit the selection of colors */
        iColor = SiForegroundColor;
    } else {
        iElement = iAtomElement(aAtom);
        if ( iElement == NOELEMENT ) {
                iColor = TGREEN;
        } else if ( iElement < 0 ) {
                fprintf(stderr, "Bad element; atomname %s iAtomElement %d\n", 
                        sAtomName( aAtom ), iAtomElement(aAtom));
                exit(1);
        } else {
                iColor = SiaElementColors[iAtomElement(aAtom)];
        }
    }
   
    gcNew = gcTankColorIndex( tTank, iColor, bFast, bSolid );
    return(gcNew);
}




/*
 *      TankPageFlip
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Flip pages for the TANK.
 */
static void
TankPageFlip( TANK tTank )
{
int             iEntries;

    iEntries = SpfPageFlip.iXColorEntries;

        /* Fill the color table with the correct color values */

    if ( !tTank->tank.bPageFlipFront ) {
        XStoreColors( SdPDisplay, ScMap,
                SpfPageFlip.xcaBackMap,
                iEntries );
    } else {
        XStoreColors( SdPDisplay, ScMap,
                SpfPageFlip.xcaFrontMap,
                iEntries );
   }
    tTank->tank.bPageFlipFront = !tTank->tank.bPageFlipFront;
    XFlush( SdPDisplay );
}






/*
 *      TankBuildArrays
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Build arrays for speeding up graphics.
 */
static void
TankBuildArrays( TANK tTank )
{
LOOP            lAtoms;
ATOM            aAtom;
int             iAtoms;
ATOM            *aPAtom;


    if ( tTank->tank.uUnit == NULL ) return;

    if ( tTank->tank.vaAtomPtrs != NULL ) 
        VarArrayDestroy(&(tTank->tank.vaAtomPtrs) );
    iAtoms = 0;
    lAtoms = lLoop( (OBJEKT)tTank->tank.uUnit, ATOMS );
    while ( oNext(&lAtoms) ) iAtoms++;
    tTank->tank.vaAtomPtrs = vaVarArrayCreate(sizeof(ATOM));
    VarArraySetSize( (tTank->tank.vaAtomPtrs), iAtoms );
    
    if ( iAtoms ) {
        aPAtom = PVAI( tTank->tank.vaAtomPtrs, ATOM, 0 );
        lAtoms = lLoop( (OBJEKT)tTank->tank.uUnit, ATOMS );
        while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
                *aPAtom = aAtom;
                aPAtom++;
        }
    }
}




/*
 *      TankBuildString
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Build a string by concatenating
 *      the different strings, atom name, type, charge.
 *      Depending on flags set within the TANK structure.
 *
 *TODO: Add support for displaying CHARGES, and RESIDUE names
 */
static void
TankBuildString( TANK tTank, ATOM aAtom, char *sStr )
{
STRING          sTemp;

    strcpy( sStr, "" );
    if ( bTankFlagsSet( tTank, TANKSHOWRESIDUES ) ) {
        if ( iContainerSequence(aAtom) == 1 ) {
            sprintf( sTemp, "<%s:%d>", 
                sContainerName(cContainerWithin(aAtom)),
                iContainerSequence(cContainerWithin(aAtom)) );
            strcat( sStr, sTemp );
        }
    }
    if ( bTankFlagsSet( tTank, TANKSHOWNAMES ) ) {
        strcat( sStr, sContainerName(aAtom) );
    }
    if ( bTankFlagsSet( tTank, TANKSHOWPERTNAMES ) ) {
        strcat( sStr, "/" );
        strcat( sStr, sAtomPertName(aAtom) );
    }
    if ( bTankFlagsSet( tTank, TANKSHOWTYPES ) ||
         bTankFlagsSet( tTank, TANKSHOWPERTTYPES ) ) {
        strcat( sStr, "(" );
        if ( bTankFlagsSet( tTank, TANKSHOWTYPES ) ) {
            strcat( sStr, sAtomType(aAtom) );
        }
        if ( bTankFlagsSet( tTank, TANKSHOWPERTTYPES ) ) {
            strcat( sStr, "/" );
            strcat( sStr, sAtomPertType(aAtom) );
        }
        strcat( sStr, ")" );
    }
    if ( bTankFlagsSet( tTank, TANKSHOWCHARGES ) ) {
        sprintf(sTemp, "%7.4f", dAtomCharge(aAtom));
        strcat( sStr, sTemp);
    }
}



/*
 *-----------------------------------------------------------------
 *
 *      Handle the virtual sphere
 *
 */
  

/*
 *      TankDrawCircle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw a circle on the screen that is used for
 *      the virtual sphere.
 *      The circle is always drawn in the foreground color
 *      and in bFast=TRUE mode.
 */
static void
TankDrawCircle( TANK tTank, Drawable dDest )
{
int             iTLX, iTLY;
int             iWidth, iHeight;

    iWidth = iX3dScreenUseWidth(tTank->tank.x3dEngine) * CIRCLERADIUS;
    iHeight = iX3dScreenUseHeight(tTank->tank.x3dEngine) * CIRCLERADIUS;

    iTLX = (iX3dScreenWidth(tTank->tank.x3dEngine) >> 1) - (iWidth >> 1);
    iTLY = (iX3dScreenHeight(tTank->tank.x3dEngine) >> 1) - (iHeight >> 1);

    XDrawArc( XtDisplay(tTank), dDest, 
        gcTankColorIndex( tTank, SiForegroundColor, TRUE, TRUE ),
        iTLX, iTLY,
        iWidth, iHeight,
        0, 23040 );
}


/*
 *      bTankPointInCircle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if the point is within the circle.
 */
BOOL
bTankPointInCircle( TANK tTank, int iX, int iY )
{
double          dRealX, dRealY;
double          dDist2, dSize2;
BOOL            bInside;
double          dWidth, dHeight;


    dWidth = (double)iX3dScreenUseWidth(tTank->tank.x3dEngine);
    dHeight= (double)iX3dScreenUseHeight(tTank->tank.x3dEngine);

    dRealX = iX-(iX3dScreenWidth(tTank->tank.x3dEngine)>>1);
    dRealY = (iY-(iX3dScreenHeight(tTank->tank.x3dEngine)>>1)) * 
                dWidth/dHeight;
    dDist2 = dRealX*dRealX + dRealY*dRealY;
    dSize2 = dWidth / 2.0 * CIRCLERADIUS;
    dSize2 = dSize2*dSize2;

                /* If the point is outside the circle return FALSE */
    if ( dDist2 < dSize2 ) bInside = TRUE;
    else bInside = FALSE;

    return(bInside);
}


/*
 *------------------------------------------------------------------------
 *
 *      Draw rubber band and rubber box.
 */

        
/*
 *      TankToggleDrawRubberBand
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw the rubber band from iStartX, iStartY 
 *      to iStopX, iStopY.
 */
void
TankToggleDrawRubberBand( TANK tTank )
{

    XDrawLine( XtDisplay(tTank), XtWindow(tTank),
                SgcRubberBand,
                tTank->tank.iXStart, tTank->tank.iYStart,
                tTank->tank.iXStop, tTank->tank.iYStop );
    tTank->tank.bRubberOn = !tTank->tank.bRubberOn;
}





        
/*
 *      TankToggleDrawRubberBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw the rubber box using iStartX,SiStartY 
 *      to iStopX, iStopY.
 */
void
TankToggleDrawRubberBox( TANK tTank )
{
    XDrawLine( XtDisplay(tTank), XtWindow(tTank),
                SgcRubberBand,
                tTank->tank.iXStart, tTank->tank.iYStart,
                tTank->tank.iXStop, tTank->tank.iYStart );
    XDrawLine( XtDisplay(tTank), XtWindow(tTank),
                SgcRubberBand,
                tTank->tank.iXStop, tTank->tank.iYStart,
                tTank->tank.iXStop, tTank->tank.iYStop );
    XDrawLine( XtDisplay(tTank), XtWindow(tTank),
                SgcRubberBand,
                tTank->tank.iXStop, tTank->tank.iYStop,
                tTank->tank.iXStart, tTank->tank.iYStop );
    XDrawLine( XtDisplay(tTank), XtWindow(tTank),
                SgcRubberBand,
                tTank->tank.iXStart, tTank->tank.iYStop,
                tTank->tank.iXStart, tTank->tank.iYStart );
    tTank->tank.bRubberOn = !tTank->tank.bRubberOn;
}



/*
 *-----------------------------------------------------------------
 *
 *      Rendering
 */
 


/*
 *      TankClearBackground
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Clear the BACKGROUND of the display
 *      prior to drawing the molecule.
 */
static void
TankClearBackground( TANK tTank, Drawable *dPDest, BOOL bFast )
{
GC              gc;
Drawable        dDest;

        /* First set up the GC and the Drawable */


    switch ( SiRenderWay ) {
        case TANKOFFSCREEN:
            dDest = tTank->tank.pOffScreen;
            break;
        case TANKPAGEFLIP:
            dDest = XtWindow(tTank);
            break;
    };
    *dPDest = dDest;


                /* Clear the background */

    switch ( SiRenderWay ) {
        case TANKOFFSCREEN:
            XFillRectangle( XtDisplay(tTank), dDest, 
                        gcTankColor( tTank, NULL, bFast, TRUE ),
                        0, 0, 
                        iX3dScreenWidth(tTank->tank.x3dEngine),
                        iX3dScreenHeight(tTank->tank.x3dEngine) );
            break;
        case TANKPAGEFLIP:

            if ( bFast != tTank->tank.bFastPageFlipLastTime ) {

                        /* If last time it was slow and now it must be */
                        /* fast then clear the TANK using the PreparePageFlip */
                        /* GC */

                if ( bFast ) {
                    XFillRectangle( XtDisplay(tTank), dDest,
                                SpfPageFlip.gcPreparePageFlip,
                                0, 0, 
                                iX3dScreenWidth(tTank->tank.x3dEngine),
                                iX3dScreenHeight(tTank->tank.x3dEngine) );
                } else {
                    XFillRectangle( XtDisplay(tTank), dDest,
                                gcTankColor( tTank, NULL, FALSE, TRUE ),
                                0, 0, 
                                iX3dScreenWidth(tTank->tank.x3dEngine),
                                iX3dScreenHeight(tTank->tank.x3dEngine) );
                }
            } else {
                if ( bFast ) {
                    gc = gcTankColor( tTank, NULL, bFast, TRUE );
                    XFillRectangle( XtDisplay(tTank), dDest, gc,
                                0, 0, 
                                iX3dScreenWidth(tTank->tank.x3dEngine),
                                iX3dScreenHeight(tTank->tank.x3dEngine) );
                    XFlush(XtDisplay(tTank));
                } else {

                            /* If last time it was slow and now it is slow */
                            /* again then just clear the color screen */
                    gc = gcTankColor( tTank, NULL, bFast, TRUE );
                    XFillRectangle( XtDisplay(tTank), dDest, gc,
                                0, 0, 
                                iX3dScreenWidth(tTank->tank.x3dEngine),
                                iX3dScreenHeight(tTank->tank.x3dEngine) );
                    XFlush(XtDisplay(tTank));
                }
            }
            break;
    }
    tTank->tank.bNeedToRedraw = FALSE;
}



/*
 *      TankDrawLoneAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw a diamond representing a lone atom.
 */
static void
TankDrawLoneAtom( TANK tTank, Drawable dDest, BOOL bFast, 
        int iX, int iY, ATOM aA )
{
int             iXd, iYd;
int             iX1, iY1;
int             iX2, iY2;
int             iX3, iY3;
int             iX4, iY4;
GC              gc;

    iXd = 3;
    iYd = 3;
    iX1 = iX;           iY1 = iY-iYd;
    iX2 = iX+iXd;       iY2 = iY;
    iX3 = iX;           iY3 = iY+iYd;
    iX4 = iX-iXd;       iY4 = iY;

    gc = gcTankColor( tTank, aA, bFast, TRUE );

    XDrawLine( XtDisplay(tTank), dDest, gc, iX1, iY1, iX2, iY2 );
    XDrawLine( XtDisplay(tTank), dDest, gc, iX2, iY2, iX3, iY3 );
    XDrawLine( XtDisplay(tTank), dDest, gc, iX3, iY3, iX4, iY4 );
    XDrawLine( XtDisplay(tTank), dDest, gc, iX4, iY4, iX1, iY1 );
}


/*
 *      TankDrawUnBuiltAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw a cross representing an unbuilt atom.
 */
static void
TankDrawUnBuiltAtom( TANK tTank, Drawable dDest, BOOL bFast, 
        int iX, int iY, ATOM aA )
{
int             iXd, iYd;
int             iX1, iY1;
int             iX2, iY2;
int             iX3, iY3;
int             iX4, iY4;
GC              gc;


    iXd = 3;
    iYd = 3;
    iX1 = iX-iXd;       iY1 = iY-iYd;
    iX2 = iX+iXd;       iY2 = iY+iYd;
    iX3 = iX-iXd;       iY3 = iY+iYd;
    iX4 = iX+iXd;       iY4 = iY-iYd;

    gc = gcTankColor( tTank, aA, bFast, TRUE );

    XDrawLine( XtDisplay(tTank), dDest, gc, iX1, iY1, iX2, iY2 );
    XDrawLine( XtDisplay(tTank), dDest, gc, iX3, iY3, iX4, iY4 );
}






/*
 *      TankDrawSingleBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw a single bond.
 */
static void
TankDrawSingleBond( TANK tTank, Drawable dDest, BOOL bFast, 
        int iXA, int iYA, int iXB, int iYB, ATOM aA, ATOM aB )
{
int             iXC, iYC;

    iXC = ( iXA + iXB ) >> 1;
    iYC = ( iYA + iYB ) >> 1;
    XDrawLine( XtDisplay(tTank), dDest,
        gcTankColor( tTank, aA, bFast, TRUE ),
        iXA, iYA, iXC, iYC );
    XDrawLine( XtDisplay(tTank), dDest,
        gcTankColor( tTank, aB, bFast, TRUE ),
        iXC, iYC, iXB, iYB );
}




/*
 *      TankDrawDoubleBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw a double bond.
 *      Currently this routine assumes an aspect ratio of 1.0
 *      this shouldn't make too much difference, but maybe
 *      it should be fixed.
 *
 */
static void
TankDrawDoubleBond( TANK tTank, Drawable dDest, BOOL bFast,
        int iXStart, int iYStart, int iXStop, int iYStop,
        ATOM aStart, ATOM aStop )
{
int             iX, iY, iX2, iY2, iDX, iDY;
GC              gcStart, gcStop;

    iX = iXStop - iXStart;
    iY = iYStart - iYStop;

    MESSAGE(( "iX, iY = %d, %d\n", iX, iY ));

    iX2 = iX * 2;
    iY2 = iY * 2;

                /* Get the offsets from the atom position for the */
                /* starting points of the bond */
    if ( abs(iX) >= abs(iY2) ) {
        iDX = 0;
        iDY = 3;
        MESSAGE(( "Horizontal\n" ));
    } else if ( abs(iY) >= abs(iX2) ) {
        iDX = 3;
        iDY = 0;
        MESSAGE(( "Vertical\n" ));
    } else if ( iY*iX > 0 ) {
        iDX = 2;
        iDY = 2;
        MESSAGE(( "Diagonal negative\n" ));
    } else {
        iDX = 2;
        iDY = -2;
        MESSAGE(( "Diagonal positive\n" ));
    }

                /* Find the point halfway between the atoms */

    iX = ( iXStart + iXStop ) / 2;
    iY = ( iYStart + iYStop ) / 2;

                /* Now draw the bonds */

    gcStart = gcTankColor( tTank, aStart, bFast, TRUE );
    gcStop  = gcTankColor( tTank, aStop, bFast, TRUE );

    XDrawLine( XtDisplay(tTank), dDest, gcStart, 
                iXStart+iDX, iYStart+iDY, iX+iDX, iY+iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStop,
                iX+iDX, iY+iDY, iXStop+iDX, iYStop+iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStart, 
                iXStart-iDX, iYStart-iDY, iX-iDX, iY-iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStop,
                iX-iDX, iY-iDY, iXStop-iDX, iYStop-iDY );
}






/*
 *      TankDrawTripleBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw a triple bond.
 *      Currently this routine assumes an aspect ratio of 1.0
 *      this shouldn't make too much difference, but maybe
 *      it should be fixed.
 *
 *      TODO: Fix this for aspect ratio's other than 1.0
 */
static void
TankDrawTripleBond( TANK tTank, Drawable dDest, BOOL bFast,
        int iXStart, int iYStart, int iXStop, int iYStop,
        ATOM aStart, ATOM aStop )
{
int             iX, iY, iX2, iY2, iDX, iDY;
GC              gcStart, gcStop;

    iX = iXStop - iXStart;
    iY = iYStart - iYStop;

    MESSAGE(( "iX, iY = %d, %d\n", iX, iY ));

    iX2 = iX * 2;
    iY2 = iY * 2;

                /* Get the offsets from the atom position for the */
                /* starting points of the bond */
    if ( abs(iX) >= abs(iY2) ) {
        iDX = 0;
        iDY = 4;
        MESSAGE(( "Horizontal\n" ));
    } else if ( abs(iY) >= abs(iX2) ) {
        iDX = 4;
        iDY = 0;
        MESSAGE(( "Vertical\n" ));
    } else if ( iY*iX > 0 ) {
        iDX = 3;
        iDY = 3;
        MESSAGE(( "Diagonal negative\n" ));
    } else {
        iDX = 3;
        iDY = -3;
        MESSAGE(( "Diagonal positive\n" ));
    }

                /* Find the point halfway between the atoms */

    iX = ( iXStart + iXStop ) / 2;
    iY = ( iYStart + iYStop ) / 2;

                /* Now draw the bonds */

    gcStart = gcTankColor( tTank, aStart, bFast, TRUE );
    gcStop  = gcTankColor( tTank, aStop, bFast, TRUE );

        /* Draw outer line 1, middle line, outer line 2 */
    XDrawLine( XtDisplay(tTank), dDest, gcStart, 
                iXStart+iDX, iYStart+iDY, iX+iDX, iY+iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStop,
                iX+iDX, iY+iDY, iXStop+iDX, iYStop+iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStart, 
                iXStart, iYStart, iX, iY );
    XDrawLine( XtDisplay(tTank), dDest, gcStop,
                iX, iY, iXStop, iYStop );
    XDrawLine( XtDisplay(tTank), dDest, gcStart, 
                iXStart-iDX, iYStart-iDY, iX-iDX, iY-iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStop,
                iX-iDX, iY-iDY, iXStop-iDX, iYStop-iDY );
}



/*
 *      TankDrawAromaticBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw an aromatic bond.
 *      Currently this routine assumes an aspect ratio of 1.0
 *      this shouldn't make too much difference, but maybe
 *      it should be fixed.
 *
 *      The side that the dashed line is drawn on is arbitrary.
 *      This should be changed so that it is on the inside of
 *      aromatic rings.
 */
static void
TankDrawAromaticBond( TANK tTank, Drawable dDest, BOOL bFast,
        int iXStart, int iYStart, int iXStop, int iYStop,
        ATOM aStart, ATOM aStop )
{
int             iX, iY, iX2, iY2, iDX, iDY;
GC              gcStart, gcStop;
GC              gcStartDashed, gcStopDashed;

    iX = iXStop - iXStart;
    iY = iYStart - iYStop;

    MESSAGE(( "iX, iY = %d, %d\n", iX, iY ));

    iX2 = iX * 2;
    iY2 = iY * 2;

                /* Get the offsets from the atom position for the */
                /* starting points of the bond */
    if ( abs(iX) >= abs(iY2) ) {
        iDX = 0;
        iDY = 3;
        MESSAGE(( "Horizontal\n" ));
    } else if ( abs(iY) >= abs(iX2) ) {
        iDX = 3;
        iDY = 0;
        MESSAGE(( "Vertical\n" ));
    } else if ( iY*iX > 0 ) {
        iDX = 2;
        iDY = 2;
        MESSAGE(( "Diagonal negative\n" ));
    } else {
        iDX = 2;
        iDY = -2;
        MESSAGE(( "Diagonal positive\n" ));
    }

                /* Find the point halfway between the atoms */

    iX = ( iXStart + iXStop ) / 2;
    iY = ( iYStart + iYStop ) / 2;

                /* Now draw the bonds */

    gcStart = gcTankColor( tTank, aStart, FALSE, TRUE );
    gcStop  = gcTankColor( tTank, aStop, FALSE, TRUE );
    gcStartDashed = gcTankColor( tTank, aStart, FALSE, FALSE );
    gcStopDashed = gcTankColor( tTank, aStop, FALSE, FALSE );

    XDrawLine( XtDisplay(tTank), dDest, gcStart, 
                iXStart+iDX, iYStart+iDY, iX+iDX, iY+iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStop,
                iX+iDX, iY+iDY, iXStop+iDX, iYStop+iDY );

                /* Draw the aromatic part, the dashed part */

    XDrawLine( XtDisplay(tTank), dDest, gcStartDashed, 
                iXStart-iDX, iYStart-iDY, iX-iDX, iY-iDY );
    XDrawLine( XtDisplay(tTank), dDest, gcStopDashed,
                iX-iDX, iY-iDY, iXStop-iDX, iYStop-iDY );
}



/*
 *      zTankDrawLine
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw a line in the tank of a certian color.
 */
static void
zTankDrawLine( TANK tTank, Drawable dDest, GVECTOR *gvPA, GVECTOR *gvPB, 
        int iColor, BOOL bFast )
{
GVECTOR         gvA, gvB;
GC              gc;
int             iXA, iYA, iXB, iYB;


    X3dWorldToScreen( tTank->tank.x3dEngine, gvPA, &gvA );
    X3dWorldToScreen( tTank->tank.x3dEngine, gvPB, &gvB );
    if ( bX3dVisibleLineClip( tTank->tank.x3dEngine,
                        &iXA, &iYA, &iXB, &iYB,
                        &gvA, &gvB ) ) {
        gc = gcTankColorIndex( tTank, iColor, bFast, TRUE );
        XDrawLine( XtDisplay(tTank), dDest, gc, 
                   iXA, iYA, iXB, iYB );
    }
}



/*
 *      zTankDrawBoundingBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw the bounding box.
 */
static void
zTankDrawBoundingBox( TANK tTank, Drawable dDest, BOOL bFast )
{
double          dX, dY, dZ;
GVECTOR         gvA, gvB;

    UnitGetBox( tTank->tank.uUnit, &dX, &dY, &dZ );

    dX /= 2.0;   dY /= 2.0;   dZ /= 2.0;
    GVectorDef( gvA, -dX, -dY, -dZ );
    GVectorDef( gvB, dX, -dY, -dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, -dX, dY, -dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, -dX, -dY, dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    
    GVectorDef( gvA, dX, dY, -dZ );
    GVectorDef( gvB, -dX, dY, -dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, dX, -dY, -dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, dX, dY, dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );

    GVectorDef( gvA, -dX, dY, dZ );
    GVectorDef( gvB, dX, dY, dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, -dX, -dY, dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, -dX, dY, -dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );

    GVectorDef( gvA, dX, -dY, dZ );
    GVectorDef( gvB, dX, dY, dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, -dX, -dY, dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
    GVectorDef( gvB, dX, -dY, -dZ );
    zTankDrawLine( tTank, dDest, &gvA, &gvB,
                    SiForegroundColor, bFast );
}



/*
 *      zTankDrawCoordinateAxis
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw the coordinate axis.
 */
static void
zTankDrawCoordinateAxis( TANK tTank, Drawable dDest, BOOL bFast )
{
int             iColor;
GC              gc;
int             iXA, iYA, iXB, iYB;
GVECTOR         gvA, gvB, gvTA, gvTB;


    iColor = SiForegroundColor;
    gc = gcTankColorIndex( tTank, iColor, bFast, TRUE );
    GVectorDef( gvA, 0.0, 0.0, 0.0 );
    GVectorDef( gvB, AXISLENGTH, 0.0, 0.0 );
    X3dWorldToScreen( tTank->tank.x3dEngine, &gvA, &gvTA );
    X3dWorldToScreen( tTank->tank.x3dEngine, &gvB, &gvTB );
    if ( bX3dVisibleLineClip( tTank->tank.x3dEngine,
                    &iXA, &iYA, &iXB, &iYB,
                    &gvTA, &gvTB ) ) {
        XDrawLine( XtDisplay(tTank), dDest, gc, 
                iXA, iYA, iXB, iYB );
        XDrawString( XtDisplay(tTank), dDest, gc, iXB, iYB, "X", 1 );
    }
    GVectorDef( gvB, 0.0, AXISLENGTH, 0.0 );
    X3dWorldToScreen( tTank->tank.x3dEngine, &gvB, &gvTB );
    if ( bX3dVisibleLineClip( tTank->tank.x3dEngine,
                    &iXA, &iYA, &iXB, &iYB,
                    &gvTA, &gvTB ) ) {
        XDrawLine( XtDisplay(tTank), dDest, gc, 
                iXA, iYA, iXB, iYB );
        XDrawString( XtDisplay(tTank), dDest, gc, iXB, iYB, "Y", 1 );
    }
    GVectorDef( gvB, 0.0, 0.0, AXISLENGTH );
    X3dWorldToScreen( tTank->tank.x3dEngine, &gvB, &gvTB );
    if ( bX3dVisibleLineClip( tTank->tank.x3dEngine,
                    &iXA, &iYA, &iXB, &iYB,
                    &gvTA, &gvTB ) ) {
        XDrawLine( XtDisplay(tTank), dDest, gc, 
                iXA, iYA, iXB, iYB );
        XDrawString( XtDisplay(tTank), dDest, gc, iXB, iYB, "Z", 1 );
    }

}




/*
 *      TankDisplay
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display the contents of the TANK within the
 *      Widget.  
 *      If bFast is TRUE then draw the image
 *      as quickly as possible, this is only used in TANKPAGEFLIP
 *      mode when there are not enough bit planes to do page flipping
 *      in color.  Instead do it in monochrome.
 *
 */
void
TankDisplay( TANK tTank, BOOL bFast )
{
GVECTOR         gvA;
ATOM            aA, aB;
int             iXA, iYA, iXB, iYB;
Display         *d;
Drawable        dDest;
ATOM            *aPA;
GRAPHATOMt      *gaPA, *gaPB;
GVECTOR         *gvPA, *gvPB;
int             iAtom, iAtomCount, iBond;
STRING          sTemp;
BOOL            bDrewBond;



    if ( !tTank->tank.bGotInitialExpose ) return;
    if ( tTank->tank.uUnit == NULL ) return;
    if ( tTank->tank.vaAtomPtrs == NULL ) 
                DFATAL(( "No vaAtomPtrs array defined" ));

                /* If we are doing page-flipping and we are changing */
                /* how we draw (fast/slow) then we must redraw even if */
                /* the transform matrix has not been changed */

    if ( SiRenderWay == TANKPAGEFLIP ) {
        if ( bFast != tTank->tank.bFastPageFlipLastTime ) {
            tTank->tank.bNeedToRedraw = TRUE;
        }
    }

/*    if ( !tTank->tank.bNeedToRedraw ) return; */

    MESSAGE(( "About to render TANK fast: %s\n", sBOOL(bFast) ));
    
    d = XtDisplay(tTank);

                /* Clear the background and get the drawable onto */
                /* which we can draw */

    TankClearBackground( tTank, &dDest, bFast );

                /* Draw the coordinate axis if they are required */

    if ( bTankFlagsSet( tTank, TANKSHOWAXIS ) ) 
        zTankDrawCoordinateAxis( tTank, dDest, bFast );
        
                /* Draw the circle if virtual sphere is on */

    if ( bTankFlagsSet( tTank, TANKDRAWCIRCLE ) ) {
        TankDrawCircle( tTank, dDest );
    }

                /* Draw the bounding box if it is required */

    if ( bTankFlagsSet( tTank, TANKSHOWBOX ) && 
                bUnitUseBox(tTank->tank.uUnit) ) 
        zTankDrawBoundingBox( tTank, dDest, bFast );

    if ( !iVarArrayElementCount(tTank->tank.vaAtomPtrs) ) 
        return;


                /* Transform the coordinates */

    iAtomCount = iVarArrayElementCount( tTank->tank.vaAtomPtrs );
    aPA = PVAI( tTank->tank.vaAtomPtrs, ATOM, 0 );
    for ( iAtom=0; iAtom<iAtomCount; iAtom++, aPA++ ) {
        if ( bAtomHasPosition(*aPA) ) {
            gaPA = (GRAPHATOMt*)PAtomGraphicsPointer(*aPA);
            GVectorDef( gvA, dVX(&vAtomPosition(*aPA)),
                      dVY(&vAtomPosition(*aPA)),
                      dVZ(&vAtomPosition(*aPA)) );
            X3dWorldToScreen( tTank->tank.x3dEngine, 
                                &gvA, &(gaPA->gvScreen) );
        }
    }



                /* Now render the bonds */

    aPA = PVAI( tTank->tank.vaAtomPtrs, ATOM, 0 );
    for ( iAtom=0; iAtom<iAtomCount; iAtom++, aPA++ ) {
        aA = *aPA;
        gaPA = (GRAPHATOMt*)PAtomGraphicsPointer(aA);
        gvPA = &(gaPA->gvScreen);
        if ( !bAtomVisible(aA) ) continue;
        bDrewBond = FALSE;
        for ( iBond=0; iBond<iAtomCoordination(aA); iBond++ ) {
            aB = aAtomBondedNeighbor( aA, iBond );
            if ( !bAtomVisible(aB) ) continue;
            bDrewBond = TRUE;
            if ( iAtomId(aA) < iAtomId(aB) ) {
                MESSAGE(( "Skipping bond\n" ));
                continue;
            }

                /* Get the coordinates */

            gaPB = (GRAPHATOMt*)PAtomGraphicsPointer(aB);
            gvPB = &(gaPB->gvScreen);
                
            if ( !bX3dVisibleLineClip( tTank->tank.x3dEngine, 
                                 &iXA, &iYA, &iXB, &iYB, 
                                 gvPA, gvPB ) ) continue;


            if ( bFast ) {
                XDrawLine( XtDisplay(tTank), dDest,
                        gcTankColor( tTank, aA, bFast, TRUE ),
                        iXA, iYA, iXB, iYB );

            } else {
                switch( iAtomBondOrder( aA, iBond ) ) {
                    case BONDSINGLE:
                        TankDrawSingleBond( tTank, dDest, bFast,
                                                iXA, iYA, iXB, iYB,
                                                aA, aB );
                        break;
                    case BONDDOUBLE:
                        TankDrawDoubleBond( tTank, dDest, bFast,
                                                iXA, iYA, iXB, iYB,
                                                aA, aB );
                        break;
                    case BONDTRIPLE:
                        TankDrawTripleBond( tTank, dDest, bFast,
                                                iXA, iYA, iXB, iYB,
                                                aA, aB );
                        break;
                    case BONDAROMATIC:
                        TankDrawAromaticBond( tTank, dDest, bFast,
                                                iXA, iYA, iXB, iYB,
                                                aA, aB );
                        break;
                    default:
                        DFATAL(( "Unknown bond type: %d\n",
                                iAtomBondOrder( aA, iBond ) ));
                        break;
                }
            }
MESSAGE(( "Drew line between: %d, %d - %d, %d\n", 
                        iXA, iYA, iXB, iYB ));
        }
        gvPA = &(gaPA->gvScreen);
        if ( cGVectorClipStatus(*gvPA) == CLIPMIDDLE ) {
            iXA = (int)nGVX(*gvPA);
            iYA = (int)nGVY(*gvPA);
            if ( bAtomFlagsSet( aA, ATOMPOSITIONDRAWN ) ) 
                TankDrawUnBuiltAtom( tTank, dDest, bFast, 
                                        iXA, iYA, aA );
            if ( !bDrewBond ) TankDrawLoneAtom( tTank, dDest, 
                                bFast, iXA, iYA, aA );
        }
    }
    XFlush(d);


        /* Display labels attached to atoms */

    if ( bTankFlagsSet( tTank, TANKSHOWNAMES ) || 
         bTankFlagsSet( tTank, TANKSHOWPERTNAMES ) ||
         bTankFlagsSet( tTank, TANKSHOWTYPES ) ||
         bTankFlagsSet( tTank, TANKSHOWPERTTYPES ) ||
         bTankFlagsSet( tTank, TANKSHOWCHARGES ) || 
         bTankFlagsSet( tTank, TANKSHOWRESIDUES ) ) {
        aPA = PVAI(tTank->tank.vaAtomPtrs, ATOM, 0 );
        for ( iAtom=0; iAtom<iAtomCount; iAtom++, aPA++ ) {
            if ( !bAtomVisible(*aPA) ) continue;
            gaPA = (GRAPHATOMt*)PAtomGraphicsPointer(*aPA);
            if ( bGVectorMiddle(gaPA->gvScreen) ) {
                iXA = nGVX(gaPA->gvScreen);
                iYA = nGVY(gaPA->gvScreen);
                TankBuildString( tTank, *aPA, sTemp );
                XDrawString( XtDisplay(tTank), dDest,
                             gcTankColor( tTank, *aPA, bFast, TRUE ),
                             iXA, iYA,
                             sTemp, strlen(sTemp) );
            }
        }
    }

        /* Now do what must be done to bring the image */
        /* onto the screen */

    switch ( SiRenderWay ) {
        case TANKOFFSCREEN:
            XCopyArea( XtDisplay(tTank), dDest, XtWindow(tTank),
                        gcTankColor( tTank, NULL, bFast, TRUE ),
                        0, 0, 
                        iX3dScreenWidth(tTank->tank.x3dEngine),
                        iX3dScreenHeight(tTank->tank.x3dEngine),
                        0, 0 );
            break;
        case TANKPAGEFLIP:
            if ( bFast || !SpfPageFlip.bDoPartialPageFlip ) {
                TankPageFlip(tTank);
            }
            tTank->tank.bFastPageFlipLastTime = bFast;
            break;
    }
}





/*
 *================================================================
 *
 *      Callback functions
 */




/*
 *      TankHandleExpose
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
static XtCallbackProc
TankHandleExpose( Widget wWidget, TANK PData, caddr_t PTemp )
{
TANK                    tTank;
int                     iX, iY;
unsigned int            uiWidth, uiHeight, uiBorder, uiDepth;
double                  dXmm, dXPixels;
double                  dYmm, dYPixels;
double                  dAspectRatio;
Cursor                  cCursor;
Window                  wWin;
XSetWindowAttributes    xswaWin;



    tTank = (TANK)wWidget;

                /* Handle the very first Expose event gotten by the TANK */
                /* Widget */

    if ( !SbInitializedClass ) {
        SbInitializedClass = TRUE;
        TankClassInitialize( XtDisplay(tTank), XtScreen(tTank) );
    }

                /* Handle the INITIAL Expose event for the Widget */

    if ( !tTank->tank.bGotInitialExpose ) {
        wWin = XtWindow(tTank);
        if ( wWin != (Window)NULL ) {
            MESSAGE(( "Initializing TANK Window properties\n" ));
            tTank->tank.bGotInitialExpose = TRUE;
            cCursor = zcTankCursor( tTank->tank.iTankState );

            xswaWin.cursor = cCursor;
            xswaWin.colormap = ScMap;
            XChangeWindowAttributes( SdPDisplay, wWin, 
                                        CWColormap | CWCursor,
                                        &xswaWin );

                        /* Tell the Window manager to acknowlege this */
                        /* windows colormap */

            XSetWMColormapWindows( SdPDisplay, wWin, &wWin, 1 );

            MESSAGE(( "Window installed with colormap: 0x%lX\n", ScMap ));
            MESSAGE(( "Default colormap: 0x%lX\n", 
                        DefaultColormapOfScreen(SsPScreen) ));
        }
    }



                /* Check if the window has been resized */

    XGetGeometry( SdPDisplay, XtWindow(tTank),
                  &SwRoot, &iX, &iY, &uiWidth, &uiHeight,
                  &uiBorder, &uiDepth );

                /* Handle a Resize */

    if ( iX3dScreenWidth(tTank->tank.x3dEngine) != uiWidth ||
         iX3dScreenHeight(tTank->tank.x3dEngine)!= uiHeight ) {
        X3dSetScreenWidth( tTank->tank.x3dEngine, uiWidth );
        X3dSetScreenHeight( tTank->tank.x3dEngine, uiHeight );

                        /* Recalculate the height and width to use to keep */
                        /* the drawing from being affected by the dimensions */
                        /* of the window. */
                        /* Do this by calculating the height and width */
                        /* of a window that would contain the drawing */
                        /* without distorsion */

        dYPixels = (double)HeightOfScreen(XtScreen(tTank));
        dYmm     = (double)HeightMMOfScreen(XtScreen(tTank));
        dXPixels = (double)WidthOfScreen(XtScreen(tTank));
        dXmm     = (double)WidthMMOfScreen(XtScreen(tTank));

        dAspectRatio = dYPixels/dYmm*dXmm/dXPixels;

                        /* If the width of the screen times the aspect */
                        /* ratio is less than the height then use the */
                        /* width of the screen and width*aspectRatio */
                        /* otherwise use height and height/aspectRatio */
        if ( dAspectRatio*uiWidth< uiHeight ) {
            X3dSetScreenUseWidth( tTank->tank.x3dEngine, uiWidth );
            X3dSetScreenUseHeight( tTank->tank.x3dEngine, 
                                        uiWidth*dAspectRatio );
        } else {
            X3dSetScreenUseWidth( tTank->tank.x3dEngine, 
                                uiWidth/dAspectRatio );
            X3dSetScreenUseHeight( tTank->tank.x3dEngine, uiWidth );
        }
        switch ( SiRenderWay ) {
            case TANKPAGEFLIP:
                break;
            case TANKOFFSCREEN:

                        /* Allocate a new PIXMAP for the offscreen */
                        /* buffer */
                if ( tTank->tank.pOffScreen != None )
                    XFreePixmap( SdPDisplay, tTank->tank.pOffScreen );
                tTank->tank.pOffScreen =
                        XCreatePixmap( SdPDisplay,
                                       XtWindow(tTank),
                                iX3dScreenWidth(tTank->tank.x3dEngine),
                                iX3dScreenHeight(tTank->tank.x3dEngine),
                                       uiDepth );
                break;
        }
    }


                /* Now display the stuff */

    X3dBuildTransform( tTank->tank.x3dEngine );
    TankDisplay( tTank , FALSE );
    return NULL;
}


static XtCallbackProc   TankCreate( Widget wRequest, TANK tTank );
static XtWidgetProc     TankDestroy( TANK tTank );


TankClassRec tankClassRec = {
  { /* core fields */
    /* superclass               */      (WidgetClass) &widgetClassRec,
    /* class_name               */      "Tank",
    /* widget_size              */      sizeof(TankRec),
    /* class_initialize         */      NULL,
    /* class_part_initialize    */      NULL,
    /* class_inited             */      FALSE,
    /* initialize               */      (VFUNCTION)TankCreate,
    /* initialize_hook          */      NULL,
    /* realize                  */      XtInheritRealize,
    /* actions                  */      SxtaaActions,
    /* num_actions              */      XtNumber(SxtaaActions),
    /* resources                */      NULL,
    /* num_resources            */      0,
    /* xrm_class                */      NULLQUARK,
    /* compress_motion          */      TRUE,
    /* compress_exposure        */      TRUE,
    /* compress_enterleave      */      TRUE,
    /* visible_interest         */      FALSE,
    /* destroy                  */      (VFUNCTION)TankDestroy,
    /* resize                   */      NULL,
    /* expose                   */      (VFUNCTION)TankHandleExpose,
    /* set_values               */      NULL,
    /* set_values_hook          */      NULL,
    /* set_values_almost        */      XtInheritSetValuesAlmost,
    /* get_values_hook          */      NULL,
    /* accept_focus             */      NULL,
    /* version                  */      XtVersion,
    /* callback_private         */      NULL,
    /* tm_table                 */      NULL,
    /* query_geometry           */      XtInheritQueryGeometry,
    /* display_accelerator      */      XtInheritDisplayAccelerator,
    /* extension                */      NULL
  },
  { /* template fields */
    /* empty                    */      0
  }
};


WidgetClass tankWidgetClass = (WidgetClass)&tankClassRec;



/*
 *      TankCreate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Creator for the TANK widget.
 */
static XtCallbackProc
TankCreate( Widget wRequest, TANK tTank )
{

    MESSAGE(( "------Creating TANK\n" ));

    tTank->tank.x3dEngine = x3dX3dCreate();

    X3dSetScreenWidth( tTank->tank.x3dEngine, -1 );
    X3dSetScreenHeight( tTank->tank.x3dEngine, -1 );
    X3dSetPerspectiveOn( tTank->tank.x3dEngine, TRUE );

    tTank->tank.iButtonState = BNONE;
    tTank->tank.iRawButtonState = 0;
    tTank->tank.uUnit         = NULL;
    tTank->tank.vaAtomPtrs       = NULL;
    tTank->tank.bGotInitialExpose = FALSE;
    tTank->tank.fFlags = TANKDEFAULTFLAGS;
    tTank->tank.iCurrentDrawingElement = CARBON;
    tTank->tank.iNextAtomNumber = 1;
    tTank->tank.vaTwistTorsions = NULL;
    tTank->tank.iTankState = TANKNOSTATE;
    tTank->tank.vaRotateAtoms = NULL;
    tTank->tank.pOffScreen = (Pixmap) NULL;

    tTank->tank.bPageFlipFront = TRUE;
    tTank->tank.bFastPageFlipLastTime = FALSE;  /* guess */
    return NULL;
}



/*
 *      TankDestroy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Destroy the TANK widget.
 */
static XtWidgetProc
TankDestroy( TANK tTank )
{

    MESSAGE(( "------Destroying TANK\n" ));

                /* Destroy the X3DENGINE */
                
    if ( tTank->tank.x3dEngine != NULL ) {
        X3dDestroy( &(tTank->tank.x3dEngine) );
    }
    
                /* Destroy the ATOMs VARARRAY */
                
    if ( tTank->tank.vaAtomPtrs != NULL ) {
        VarArrayDestroy(&(tTank->tank.vaAtomPtrs) );
    }

    return NULL;
}









/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *      Maintain lists of torsions to twist.
 *      Maintain lists of ATOMs to rotate.
 */


/*
 *      zTankStartTwistTorsions
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Find all the torsions that we are going to rotate
 *      around and store them in a list so that we can
 *      do the rotations quickly.
 */
static void
zTankStartTwistTorsions( TANK tTank )
{
LOOP            lTorsions, lSpan;
ATOM            aA, aB, aC, aD, aAtom;
int             iCountB, iCountC;
TWISTTORSIONt   rtSave;
BOOL            bInRing;


    if ( tTank->tank.vaTwistTorsions != NULL ) {
        VarArrayDestroy(&(tTank->tank.vaTwistTorsions) );
        tTank->tank.vaTwistTorsions = NULL;
    }

                /* Create a VARARRAY to store the information required */
                /* for each torsion to rotate the ATOMs around that TORSION */

    tTank->tank.vaTwistTorsions = 
                vaVarArrayCreate( sizeof(TWISTTORSIONt) );

                /* Look over all the selected PROPERS and */
                /* find the side of the torsion that has the shortest */
                /* chain, (or if it is in a ring) and store the information */
                /* required to rotate only that side */

    bInRing = FALSE;
    lTorsions = lLoop( (OBJEKT)tTank->tank.uUnit, PROPERS );
    while ( oNext(&lTorsions) ) {
        LoopGetTorsion( &lTorsions, &aA, &aB, &aC, &aD );
        if ( bAtomFlagsSet( aA, ATOMSELECTED ) &&
             bAtomFlagsSet( aB, ATOMSELECTED ) &&
             bAtomFlagsSet( aC, ATOMSELECTED ) &&
             bAtomFlagsSet( aD, ATOMSELECTED ) ) {

                        /* First check the size of the one side of the */
                        /* torsion */
            iCountB = 0;
            lSpan = lLoop( (OBJEKT)aB, SPANNINGTREE );
            LoopDefineInvisibleAtom( &lSpan, aC );
            while ( (aAtom = (ATOM)oNext(&lSpan)) ) {
                if ( aAtom == aC ) {
                    bInRing = TRUE;
                    break;
                }
                iCountB++;
            }   
            if ( bInRing ) break;

                        /* Check the size of the other side of the torsion */
            iCountC = 0;
            lSpan = lLoop( (OBJEKT)aC, SPANNINGTREE );
            LoopDefineInvisibleAtom( &lSpan, aB );
            while ( (aAtom = (ATOM)oNext(&lSpan)) ) iCountC++;

                /* Set up the information for rotating the torsion */

            if ( bInRing ) {
                rtSave.aAtomStart = aB;
                rtSave.aAtomInvisible = aC;
                rtSave.bInRing = TRUE;
            } else {
                if ( iCountB > iCountC ) {
                    rtSave.aAtomStart = aC;
                    rtSave.aAtomInvisible = aB;
                    rtSave.bInRing = FALSE;
                } else {
                    rtSave.aAtomStart = aB;
                    rtSave.aAtomInvisible = aC;
                    rtSave.bInRing = FALSE;
                }
            }

            MESSAGE(( "Setting up rotate torsion: %s -> %s   : inRing: %d\n",
                        sAtomName(rtSave.aAtomInvisible),
                        sAtomName(rtSave.aAtomStart),
                        rtSave.bInRing ));

            VarArrayAdd( tTank->tank.vaTwistTorsions, (GENP)&rtSave );

        }
    }
}




        
            
            

/*
 *      zTankStopTwistTorsions
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Find all the torsions that we are going to rotate
 *      around and store them in a list so that we can
 *      do the rotations quickly.
 */
static void
zTankStopTwistTorsions( TANK tTank )
{
    VarArrayDestroy(&(tTank->tank.vaTwistTorsions) );
    tTank->tank.vaTwistTorsions = NULL;
}



/*
 *      zTankStartRotateAtoms
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Find all the ATOMs that are selected and put them in
 *      the vaRotateAtoms array.  Then find the ATOM closest to
 *      the geometric center of the ATOMs to rotate.
 */
static void
zTankStartRotateAtoms( TANK tTank )
{
LOOP            lAtoms;
ATOM            aAtom, aMin;
VECTOR          vCenter, vPos, vTrans;
int             iCount;
double          dMin, dLen;


    MESSAGE(( "++++++Start rotate atoms\n" ));
    iCount = 0;
    VectorDef( &vCenter, 0.0, 0.0, 0.0 );
    tTank->tank.vaRotateAtoms = vaVarArrayCreate(sizeof(ATOM));
    lAtoms = lLoop( (OBJEKT)tTank->tank.uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
            vCenter = vVectorAdd( &vCenter, &vAtomPosition(aAtom) );
            iCount++;
            VarArrayAdd( tTank->tank.vaRotateAtoms, (GENP)&aAtom );
        }
    }
    if ( iCount > 0 ) {
        vCenter = vVectorTimesScalar( &vCenter, 1/(double)iCount );
    }
    dMin = 9999999.0;
    lAtoms = lLoop( (OBJEKT)tTank->tank.uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
            vPos = vAtomPosition(aAtom);
            vTrans = vVectorSub( &vPos, &vCenter );
            dLen = dVectorLen(&vTrans);
            if ( dLen < dMin ) {
                dMin = dLen;
                aMin = aAtom;
            }
        }
    }
    tTank->tank.aRotateCenter = aMin;
}


        
            
            

/*
 *      zTankStopRotateAtoms
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
static void
zTankStopRotateAtoms( TANK tTank )
{
    MESSAGE(( "++++++Stop rotate ATOMS\n" ));
    VarArrayDestroy(&(tTank->tank.vaRotateAtoms) );
    tTank->tank.vaRotateAtoms = NULL;
}





/*
 *====================================================
 *====================================================
 *====================================================
 *
 *      Public routines.
 */






/*
 *      TankUseUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Inform the TANK which UNIT it should display.
 *      This routine should be called everytime the UNIT
 *      is changed outside of the TANK routines.
 */
void
TankUseUnit( TANK tTank, UNIT uUnit )
{
VECTOR          vCenter;
GVECTOR         gvCenter;


    tTank->tank.uUnit = uUnit;

    if ( uUnit == NULL ) return;

                /* Build a matrix that will cause the geometric center */
                /* of the UNIT to lie on the origin */

    vCenter = vContainerGeometricCenter( (CONTAINER)uUnit ); 
    
    GVectorDef( gvCenter,
                        dVX(&vCenter),
                        dVY(&vCenter),
                        dVZ(&vCenter) );
    cGVectorClipStatus( gvCenter ) = 0;                 /* for Purify */

    X3dSetCenterOfRotation( tTank->tank.x3dEngine, gvCenter );

    MESSAGE(( "Using UNIT: %s\n", sContainerName(uUnit) ));
    MESSAGE(( "Geometric center @ %f, %f, %f\n",
                nGVX(gvX3dCenterOfRotation(tTank->tank.x3dEngine)),
                nGVY(gvX3dCenterOfRotation(tTank->tank.x3dEngine)),
                nGVZ(gvX3dCenterOfRotation(tTank->tank.x3dEngine)) ));

    TankBuildArrays( tTank );
    X3dInit( tTank->tank.x3dEngine );
    X3dBuildTransform( tTank->tank.x3dEngine );
    X3dBuildCenteredScaledTransform( tTank->tank.x3dEngine, 
                                        tTank->tank.uUnit );

    tTank->tank.iTankSelectLastState = TANKSELECTNOTHING;
    tTank->tank.tLastClickMSec = 0;

    TankDisplay( tTank , FALSE );
}



/*
 *      TankRedisplayUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Inform the TANK to recalculate everything for the UNIT
 *      and then re-display it.
 */
void
TankRedisplayUnit( TANK tTank )
{
    if ( tTank->tank.uUnit == NULL ) return;
    TankBuildArrays( tTank );
    X3dBuildTransform( tTank->tank.x3dEngine );
    TankDisplay( tTank, FALSE );
}



/*
 *      TankFastRedisplayUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Inform the TANK to recalculate everything for the UNIT
 *      and then re-display it.
 */
void
TankFastRedisplayUnit( TANK tTank )
{
    if ( tTank->tank.uUnit == NULL ) return;
    tTank->tank.bNeedToRedraw = TRUE;
    TankDisplay( tTank, TRUE );
}





/*
 *      uTankUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the UNIT that the TANK is currently displaying.
 */
UNIT
uTankUnit( TANK tTank )
{
    return(tTank->tank.uUnit);
}



/*
 *      TankBeep
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Beep the speaker, or flash the TANK to indicate an error.
 */
void
TankBeep( TANK tTank )
{
    XBell( SdPDisplay, 100 );
}



/*
 *      TankSetState
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the state of the TANK.  This determines whether
 *      the user draws/erases/selects using the pointer button.
 */
void
TankSetState( TANK tTank, int iNewState )
{
Cursor          cCursor;
Window          wWin;

                /* Setup and breakdown rotate torsions information if */
                /* we are entering or leaving the TANKTWIST state */


    switch ( tTank->tank.iTankState ) {
        case TANKTWIST:
            zTankStopTwistTorsions(tTank);
            break;
        case TANKDRAGROTATE:
            zTankStopRotateAtoms(tTank);
            break;
        default:
            break;
    }

    switch ( iNewState ) {
        case TANKTWIST:
            zTankStartTwistTorsions(tTank);
            break;
        case TANKDRAGROTATE:
            zTankStartRotateAtoms(tTank);
            break;
        default:
            break;
    }

    tTank->tank.iTankState = iNewState;

    cCursor = zcTankCursor( iNewState );
    wWin = XtWindow(tTank);
    if ( wWin != (Window)NULL ) {
        XDefineCursor( SdPDisplay, wWin, cCursor );
    }

}



/*
 *      iTankState
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the current state of the TANK.
 */
int
iTankState( TANK tTank )
{
    return(tTank->tank.iTankState);
}




/*
 *      TankSetFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the flags controlling different properties of the TANK.
 */
void
TankSetFlags( TANK tTank, FLAGS fFlags )
{
    tTank->tank.fFlags |= fFlags;
}


/*
 *      TankResetFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Reset the flags controlling different properties of the TANK.
 */
void
TankResetFlags( TANK tTank, FLAGS fFlags )
{
    tTank->tank.fFlags &= ~fFlags;
}


/*
 *      bTankFlagsSet
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if the TANK flags are set.
 */
BOOL
bTankFlagsSet( TANK tTank, FLAGS fFlags )
{
    return((( tTank->tank.fFlags & fFlags ) == fFlags ));
}


/*
 *      TankSetDrawingElement
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the current drawing element for the TANK.
 */
void
TankSetDrawingElement( TANK tTank, int iElement )
{
    tTank->tank.iCurrentDrawingElement = iElement;
}


/*
 *      iTankDrawingElement
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the current drawing element.
 */
int
iTankDrawingElement( TANK tTank )
{
    return(tTank->tank.iCurrentDrawingElement);
}



/*
 *      TankDefinePrintSink
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Define where output goes from the TANK
 */
void
TankDefinePrintSink( TANK tTank, int iSink )
{
    tTank->tank.iPrintSink = iSink;
}




/*
 *      TankAtomGraphicsDataCreator
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Allocate memory for a structure that
 *      the ATOM can point to that will contain graphics
 *      rendering information.
 */
void
TankAtomGraphicsDataCreator( ATOM aAtom )
{
GRAPHATOMt      *gaPScreen;

    MALLOC( gaPScreen, GRAPHATOMt*, sizeof(GRAPHATOMt) );
    memset( gaPScreen, 0, sizeof(GRAPHATOMt) );         /* for Purify */
    AtomSetGraphicsPointer( aAtom, gaPScreen );

}



/*
 *      TankAtomGraphicsDataDestructor
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Destroy the Graphics data.
 */
void
TankAtomGraphicsDataDestructor( ATOM aAtom )
{
GRAPHATOMt      *gaPScreen;

    gaPScreen = (GRAPHATOMt*)PAtomGraphicsPointer(aAtom);
    FREE( gaPScreen );
    AtomSetGraphicsPointer( aAtom, NULL );

}


/*
 *      x3dTankX3dEngine
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the X3dEngine for this TANK.
 */
X3DENGINE
x3dTankX3dEngine( TANK tTank )
{
    return( tTank->tank.x3dEngine );
}




/*
 *      TankRedisplayRecenteredRescaledUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Recenter and re-scale the UNIT and
 *      then redraw it.
 */
void
TankRedisplayRecenteredRescaledUnit( TANK tTank )
{
    X3dBuildCenteredScaledTransform(
                        (tTank)->tank.x3dEngine,
                        (tTank)->tank.uUnit );
    TankDisplay( tTank, FALSE );
}

