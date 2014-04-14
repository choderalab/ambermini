#ifndef __RGBtoHLS__H
#define __RGBtoHLS__H

#ifndef _
#define _(args) args
#endif

extern void RGBtoHLS _( (double,
			 double,
			 double,
			 double*,
			 double*,
			 double*)
		       );
		       
extern void HLStoRGB _( (double*,
			 double*,
			 double*,
			 double,
			 double,
			 double)
		       );

#ifdef NeedShadowColor

extern Pixel  TopShadowColor(
 Widget  /* self */,
 Pixel   /* base */
);

extern Pixel  BottomShadowColor(
 Widget  /* self */,
 Pixel   /* base */
);

#endif /* NeedShadowColor */
#endif /* __RGBtoHLS__H */
