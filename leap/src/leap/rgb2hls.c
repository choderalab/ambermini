#include <stdio.h>

#define NeedShadowColor


#ifdef NeedShadowColor
#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#endif /* NeedShadowColor */

#define top_color_light          0.25   /* top+bottom < 0.5 */
#define bottom_color_light       0.2    /* top > bottom    */
#define top_color_name           ("gray75")
#define bottom_color_name        ("gray50")


#define assign_max(a,b) (a) = (a) > (b) ?  (a) : (b)
#define assign_min(a,b) (a) = (a) < (b) ?  (a) : (b)


void 
#if NeedFunctionPrototypes
RGBtoHLS (double red, double green, double blue,
	  double *hue, double *lightness, double*saturation)
#else
RGBtoHLS (red, green, blue, hue, lightness, saturation)
     double red, green, blue;
     double *hue, *lightness, *saturation;
#endif
{
  double brgbr[5], s1, s2, max, min;
  int imax, i;
  static double base[3] = {0.0, 2.0, 4.0};

#define set(colour,i)     \
  {                       \
    double val = colour;  \
    assign_max(val,0.0);  \
    assign_min(val,1.0);  \
    brgbr[i + 1] = val;   \
  }

  set(red,  0);
  set(green,1);
  set(blue, 2);

#undef set
  
  brgbr[0] = brgbr[3];
  brgbr[4] = brgbr[1];
      
  max = brgbr[1];
  min = brgbr[1];
  imax = 1;
  
  for(i = 2; i <= 3; i++)
    {
      if(brgbr[i] > max) {     max = brgbr[i]; imax = i; }
      else if(brgbr[i] < min)  min = brgbr[i];
    }
  
  (*lightness) = (brgbr[1] + brgbr[2] + brgbr[3])/3.0;
  
  if(max == min)
    {
      (*saturation) = 0.0;
      (*hue)        = 0.0;
    }
  else
    {
      s1 = 1.0 - min/(*lightness);
      s2 = 1.0 - (1.0 - max)/(1.0 - (*lightness));
		      
      (*saturation) = s1 > s2 ? s1 : s2;
		      
      assign_max(*saturation, 0.0);
      assign_min(*saturation, 1.0);
		      
      (*hue) = base[imax - 1] +
	(brgbr[imax + 1] - brgbr[imax - 1])/(max - min);
      
      if( (*hue) < 0.0) (*hue) += 6.0;
    }
  
  (*hue) = (*hue) * 60.0;
  
}


void 
#if NeedFunctionPrototypes
HLStoRGB (double *red, double *green, double *blue,
	  double hue, double lightness, double saturation)
#else
HLStoRGB (red, green, blue, hue, lightness, saturation)
     double *red, *green, *blue;
     double hue, lightness, saturation;
#endif
{
  double rgb[3];
  double min, max, Mid, h, l, s;
  int i;
  static struct Index
  {
    double sign, base;
    int   min, mid, max;
  } index[6] =  {
    /* sign base min mid max */
     1.0, 0.0,  2,  1,  0,
    -1.0, 2.0,  2,  0,  1,
     1.0, 2.0,  0,  2,  1,
    -1.0, 4.0,  0,  1,  2,
     1.0, 4.0,  1,  0,  2,
    -1.0, 6.0,  1,  2,  0,
  };


  rgb[0] = *red;
  rgb[1] = *green;
  rgb[2] = *blue;
  
  hue /= 60.0;

  while(hue >= 6.0) hue -= 6.0;
  while(h < 0.0)  hue += 6.0;

  lightness *= 3.0;
  if(lightness >= 3.0)
    {
      rgb[0] = rgb[1] = rgb[2] = 1.0;
      return;
    }
  else if(lightness <= 0.0)
    {
      rgb[0] = rgb[1] = rgb[2] = 0.0;
      return;
  }


  if(saturation > 1.0)
    saturation = 1.0;
  else if(saturation <= 0.0)
    {
      rgb[0] = rgb[1] = rgb[2] = lightness/3.0;
      return;
    }

  i = hue;
  Mid = index[i].sign * (hue - index[i].base);
  
  if(Mid >= (lightness - 1.0))
    {
      rgb[index[i].min] = min = lightness/3.0 * (1.0 - saturation);
      rgb[index[i].max] = max = min + (lightness - 3.0*min)/(1.0 + Mid);
      rgb[index[i].mid] = lightness - min - max;
    }
  else
    {
      rgb[index[i].max] = max = 1.0 - (1.0 - saturation)*(1.0 - lightness/3.0);
      rgb[index[i].min] = min = (lightness - max*(1.0 + Mid))/(2.0 - Mid);
      rgb[index[i].mid] = lightness - max - min;
    }

  *red   = rgb[0];
  *green = rgb[1];
  *blue  = rgb[2];

}

#ifdef NeedShadowColor

Pixel
#if NeedFunctionPrototypes
  TopShadowColor( Widget self,
		 Pixel  base)
#else
TopShadowColor(self, base)
     Widget self;
     Pixel  base;
#endif
{
  Colormap colormap;
  XColor color, dummy;
  
  if (XtIsRealized(self))
    colormap = self->core.colormap;
  else
    colormap = DefaultColormapOfScreen(XtScreen(self));

  color.pixel = base;
  
  XQueryColor(XtDisplay(self), colormap, &color);

  {
     double red, green, blue;
     double hue, lightness, saturation;

     red   = (float)color.red  /65535.0;
     green = (float)color.green/65535.0,
     blue  = (float)color.blue /65535.0;

     RGBtoHLS (red, green, blue, &hue, &lightness, &saturation);

     if ( lightness > (1.0 - top_color_light) )
       lightness -= top_color_light;
     else
       lightness += top_color_light;

     HLStoRGB (&red, &green, &blue, hue, lightness, saturation);

     color.red   = red   * 65535.0;
     color.green = green * 65535.0;
     color.blue  = blue  * 65535.0;

   }

  if (! XAllocColor(XtDisplay(self), colormap, &color))
    (void)XAllocNamedColor(XtDisplay(self), colormap, top_color_name,
			   &dummy,&color);
  
  return color.pixel;
  
}


Pixel
#if NeedFunctionPrototypes
  BottomShadowColor( Widget self,
		    Pixel  base)
#else
BottomShadowColor(self, base)
     Widget self;
     Pixel  base;
#endif
{
  Colormap colormap;
  XColor color, dummy;
  
  if (XtIsRealized(self))
    colormap = self->core.colormap;
  else
    colormap = DefaultColormapOfScreen(XtScreen(self));

  color.pixel = base;
  
  XQueryColor(XtDisplay(self), colormap, &color);

  {
     double red, green, blue;
     double hue, lightness, saturation;

     red   = (float)color.red  /65535.0;
     green = (float)color.green/65535.0,
     blue  = (float)color.blue /65535.0;

     RGBtoHLS (red, green, blue, &hue, &lightness, &saturation);

     if ( lightness < bottom_color_light )
       lightness += bottom_color_light;
     else
       lightness -= bottom_color_light;

     HLStoRGB (&red, &green, &blue, hue, lightness, saturation);

     color.red   = red   * 65535.0;
     color.green = green * 65535.0;
     color.blue  = blue  * 65535.0;

   }

  if (! XAllocColor(XtDisplay(self), colormap, &color))
    (void)XAllocNamedColor(XtDisplay(self), colormap, bottom_color_name,
			   &dummy,&color);
  
  return color.pixel;
  
}



#undef assign_max
#undef assign_min

#endif /* NeedShadowColor */
