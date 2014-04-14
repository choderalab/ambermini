# include "common.h"
# define PAI 3.14159265359
# define HALFPAI  1.570796326795
double distance(ATOM at1, ATOM at2)
{
	double dist;
	dist = (at1.x - at2.x) * (at1.x - at2.x);
	dist += (at1.y - at2.y) * (at1.y - at2.y);
	dist += (at1.z - at2.z) * (at1.z - at2.z);
	dist = sqrt(dist);
	return dist;
}

double rotate(ATOM obj1, ATOM obj2, ATOM obj3, ATOM * obj4)
{

/*
This function is used to calculate the coordinate of the fourth atom in
A-B-C-D system when we perform rotation bond B-C a certain degree
The function can calculate the twist angle of four atoms before rotating
and if you want twist angle be changed, please give you value(twist), the
varible AngleChange should larger than 0
*/

	double twist = 0.0;
	double w1, w2, w3, t;
	double ya5, xc5, yd5, zd5;
	double xa, ya, za, xc, yc, zc, xd, yd, zd;
	double xa3, ya3, za3, xc3, zc3, xd3, yd3, zd3;
	double ya4, za4, xc4, yd4, zd4;
	/* double xd4 ; */

	twist = twist / 180 * PAI;
	xa = obj1.x - obj2.x;
	ya = obj1.y - obj2.y;
	za = obj1.z - obj2.z;
	xc = obj3.x - obj2.x;
	yc = obj3.y - obj2.y;
	zc = obj3.z - obj2.z;
	xd = (*obj4).x - obj2.x;
	yd = (*obj4).y - obj2.y;
	zd = (*obj4).z - obj2.z;
	if (yc == 0)
		w1 = 0.0;
	else if (xc >= 0 && yc >= 0)
		w1 = acos(xc / sqrt(xc * xc + yc * yc));
	else if (xc <= 0 && yc >= 0)
		w1 = -acos(-xc / sqrt(xc * xc + yc * yc)) + PAI;
	else if (xc <= 0 && yc <= 0)
		w1 = acos(-xc / sqrt(xc * xc + yc * yc)) - PAI;
	else if (xc >= 0 && yc <= 0)
		w1 = -acos(xc / sqrt(xc * xc + yc * yc));

	/* rotating all the point clockwise about Z axis a w1 degree
	   then the atom C is lying in the Z-X plane */

	xa3 = xa * cos(w1) + ya * sin(w1);
	ya3 = -xa * sin(w1) + ya * cos(w1);
	za3 = za;
	xc3 = xc * cos(w1) + yc * sin(w1);

	/*  yc3=-xc*sin(w1)+yc*cos(w1);
	   this value is not use anymore */

	zc3 = zc;
	xd3 = xd * cos(w1) + yd * sin(w1);
	yd3 = -xd * sin(w1) + yd * cos(w1);
	zd3 = zd;

	/* rotating all  point clockwise about Y axis a w2 degree then
	   atom C is lying the X axis */

	if (zc3 >= 0 && xc3 >= 0)
		w2 = -acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 >= 0)
		w2 = acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 <= 0)
		w2 = -acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) + PAI;
	else if (zc3 >= 0 && xc3 <= 0)
		w2 = acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) - PAI;
	za4 = za3 * cos(w2) + xa3 * sin(w2);
	ya4 = ya3;
	xc4 = -zc3 * sin(w2) + xc3 * cos(w2);
	zd4 = zd3 * cos(w2) + xd3 * sin(w2);
	yd4 = yd3;
	/*
	   xd4 = -zd3 * sin(w2) + xd3 * cos(w2);
	 */
	/* The following program rotates A,D clockwise about X axis by w3
	   degrees, then the atom A is lying in the X-Y plane */

	if (ya4 >= 0 && za4 <= 0)
		w3 = -acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 >= 0 && za4 >= 0)
		w3 = acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 <= 0 && za4 >= 0)
		w3 = -acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) + PAI;
	else if (ya4 <= 0 && za4 <= 0)
		w3 = acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) - PAI;

	/*  xa5=xa4; */

	ya5 = ya4 * cos(w3) + za4 * sin(w3);
	xc5 = xc4;
	yd5 = yd4 * cos(w3) + zd4 * sin(w3);
	zd5 = -yd4 * sin(w3) + zd4 * cos(w3);

	if (xc5 > 0) {
		if (ya5 == 0)
			ya5 += 0.0000000001;
		if (yd5 == 0)
			yd5 += 0.0000000001;
		if ((ya5 > 0 && yd5 > 0) || (ya5 < 0 && yd5 < 0)) {
			/* A and D are on the same side of B-C */
			if (zd5 == 0)
				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = PAI - acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
		} else {

			/* A, D are not on the same side of B-C bonds, so the initial 
			   degree is PAI or -3.1415926 */

			if (zd5 == 0)
				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
			if (zd5 > 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = -acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) + PAI;
		}
	} else {
		if (xc5 == 0)
			xc5 += 0.0000000001;
		if (yd5 == 0)
			yd5 += 0.0000000001;
		if (ya5 == 0)
			ya5 += 0.0000000001;
		if ((ya5 > 0 && yd5 > 0) || (ya5 < 0 && yd5 < 0)) {
			/* A and D are on the same side of B-C */
			if (zd5 == 0)
				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = PAI - acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
		} else {

			/* A, D are not on the same side of bond B-C */

			if (zd5 == 0)
				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = -acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) + PAI;
			if (zd5 > 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
		}
	}
	t = t * 180 / PAI;
	/* printf("\n The twist angle is %8.3f\n",t); */
	return t;
}


double anglecal(ATOM obj1, ATOM obj2, ATOM obj3)
{

/*
This function is used to calculate the coordinate of the fourth atom in
A-B-C-D system when we perform rotation bond B-C a certain degree
The function can calculate the twist angle of four atoms before rotating
and if you want twist angle be changed, please give you value(twist), the
variable AngleChange should larger than 0
*/

	double angle;
	double d1, d2, d3;
	d1 = (obj1.x - obj2.x) * (obj1.x - obj2.x);
	d1 += (obj1.y - obj2.y) * (obj1.y - obj2.y);
	d1 += (obj1.z - obj2.z) * (obj1.z - obj2.z);
	d2 = (obj1.x - obj3.x) * (obj1.x - obj3.x);
	d2 += (obj1.y - obj3.y) * (obj1.y - obj3.y);
	d2 += (obj1.z - obj3.z) * (obj1.z - obj3.z);
	d3 = (obj3.x - obj2.x) * (obj3.x - obj2.x);
	d3 += (obj3.y - obj2.y) * (obj3.y - obj2.y);
	d3 += (obj3.z - obj2.z) * (obj3.z - obj2.z);
	angle = d1 + d3 - d2;
	angle = angle / (2.0 * sqrt(d1) * sqrt(d3));
	if (angle < -1)
		angle = -1;
	if (angle > 1)
		angle = 1;
	if (angle >= 0)
		angle = acos(angle);
	else
		angle = PAI - acos(-angle);
	angle = angle * 180.0 / PAI;
	return angle;
}

void euler(double ang1, double ang2, double ang3, double R[3][3])
{
	double s1, s2, s3, c1, c2, c3;
	c1 = cos(ang1);
	c2 = cos(ang2);
	c3 = cos(ang3);
	s1 = sin(ang1);
	s2 = sin(ang2);
	s3 = sin(ang3);
	R[0][0] = c1 * c2 * c3 - s1 * s3;
	R[0][1] = -c1 * c2 * s3 - s1 * c3;
	R[0][2] = c1 * s2;
	R[1][0] = s1 * c2 * c3 + c1 * s3;
	R[1][1] = -s1 * c2 * s3 + c1 * c3;
	R[1][2] = s1 * s2;
	R[2][0] = -s2 * c3;
	R[2][1] = s2 * s3;
	R[2][2] = c2;
}

void rotate0(ATOM obj1, ATOM obj2, ATOM obj3, ATOM * obj4,
			 double bond, double angle, double twist)
{

/*
This function is used to calculate the coordinate of the fourth atom in
A-B-C-D system when we perform rotation bond B-C a certain degree
The function can calculate the twist angle of four atoms before rotating
and if you want twist angle be changed, please give you value(twist), the
varible AngleChange should larger than 0
*/


	/*  
	   double w1,w2,w3,t;
	   double ya5,xc5,xd5, yd5,zd5;
	   double xa,ya,za,xc,yc,zc,xd,yd,zd;
	   double xa3,ya3,za3,xc3,zc3,xd3,yd3,zd3;
	   double ya4,za4,xc4,xd4,yd4,zd4;
	   double xx,yy,zz,xxx,yyy,zzz;
	 */
	/*
	   double w1,w2,w3,t;
	   double ya5,xc5,xd5, yd5,zd5;
	   double xa,ya,za,xc,yc,zc,xd,yd,zd;
	   double xa3,ya3,za3,xc3,zc3,xd3,yd3,zd3;
	   double ya4,za4,xc4,xd4,yd4,zd4;
	   double xx,yy,zz,xxx,yyy,zzz;
	 */
	double w1, w2, w3;
	double ya5, xc5, xd5, yd5, zd5;
	double xa, ya, za, xc, yc, zc;
	double xa3, ya3, za3, xc3, zc3;
	double ya4, za4, xc4;
	double xx, yy, zz, xxx, yyy, zzz;
	twist = twist / 180.0 * PAI;
	angle = angle / 180.0 * PAI;

	xa = obj1.x - obj2.x;
	ya = obj1.y - obj2.y;
	za = obj1.z - obj2.z;
	xc = obj3.x - obj2.x;
	yc = obj3.y - obj2.y;
	zc = obj3.z - obj2.z;
	if ((yc * yc + xc * xc) == 0)
		w1 = 0.0;
	else if (xc >= 0 && yc >= 0)
		w1 = acos(xc / sqrt(xc * xc + yc * yc));
	else if (xc <= 0 && yc >= 0)
		w1 = -acos(-xc / sqrt(xc * xc + yc * yc)) + PAI;
	else if (xc <= 0 && yc <= 0)
		w1 = acos(-xc / sqrt(xc * xc + yc * yc)) - PAI;
	else if (xc >= 0 && yc <= 0)
		w1 = -acos(xc / sqrt(xc * xc + yc * yc));

	/* rotating all the point clockwise about Z axis a w1 degree
	   then the atom C is lying in the Z-X plane */

	xa3 = xa * cos(w1) + ya * sin(w1);
	ya3 = -xa * sin(w1) + ya * cos(w1);
	za3 = za;
	xc3 = xc * cos(w1) + yc * sin(w1);

	/*  yc3=-xc*sin(w1)+yc*cos(w1);
	   this value is not use anymore */

	zc3 = zc;

	/* rotating all  point clockwise about Y axis a w2 degree then
	   atom C is lying the X axis */

	if ((xc3 * xc3 + zc3 * zc3) == 0)
		w2 = 0.0;
	else if (zc3 >= 0 && xc3 >= 0)
		w2 = -acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 >= 0)
		w2 = acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 <= 0)
		w2 = -acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) + PAI;
	else if (zc3 >= 0 && xc3 <= 0)
		w2 = acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) - PAI;
	za4 = za3 * cos(w2) + xa3 * sin(w2);
	ya4 = ya3;
	xc4 = -zc3 * sin(w2) + xc3 * cos(w2);

	/* The following program  rotating A,D clockwise about X axis a w3 degree
	   then the atom A is lying in the X-Y plane */

	if ((ya4 * ya4 + za4 * za4) == 0)
		w3 = 0.0;

	else if (ya4 >= 0 && za4 <= 0)
		w3 = -acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 >= 0 && za4 >= 0)
		w3 = acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 <= 0 && za4 >= 0)
		w3 = -acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) + PAI;
	else if (ya4 <= 0 && za4 <= 0)
		w3 = acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) - PAI;

	/*  xa5=xa4; */

	ya5 = ya4 * cos(w3) + za4 * sin(w3);
	xc5 = xc4;
	if (xc5 > 0) {
		xd5 = xc5 + bond * cos(PAI - angle);
		if (ya5 > 0)
			yd5 = bond * sin(PAI - angle);
		else
			yd5 = -bond * sin(PAI - angle);
	} else {
		xd5 = xc5 - bond * cos(PAI - angle);
		if (ya5 > 0)
			yd5 = bond * sin(PAI - angle);
		else
			yd5 = -bond * sin(PAI - angle);
	}
	zd5 = 0.0;

/* The following program restore the coordinate of atom D to
the original coordinate system which mean the first three atoms'
coordinates is unchangable */

	xx = xd5;
	yy = yd5 * cos(twist) - zd5 * sin(twist);
	zz = zd5 * cos(twist) + yd5 * sin(twist);
	yyy = yy * cos(w3) - zz * sin(w3);
	zzz = yy * sin(w3) + zz * cos(w3);
	xxx = xx;
	zz = zzz * cos(w2) - xxx * sin(w2);
	yy = yyy;
	xx = zzz * sin(w2) + xxx * cos(w2);
	xxx = xx * cos(w1) + yy * sin(-w1) + obj2.x;
	yyy = xx * sin(w1) + yy * cos(w1) + obj2.y;
	zzz = zz + obj2.z;
	(*obj4).x = xxx;
	(*obj4).y = yyy;
	(*obj4).z = zzz;
}

double angle0(ATOM obj1, ATOM obj2, ATOM obj3, double *ww1,
			  double *ww2, double *ww3)
{
	double w1, w2, w3, w4;
	double xa, ya, za, xc, yc, zc;
	double xa3, ya3, za3, xc3, zc3;
	/* double yc3; */
	double xa4, ya4, za4, xc4;
	/* double yc4, zc4; */
	double xa5, ya5, xc5;
	/* double za5,yc5, zc5; */

	xa = obj1.x - obj2.x;
	ya = obj1.y - obj2.y;
	za = obj1.z - obj2.z;
	xc = obj3.x - obj2.x;
	yc = obj3.y - obj2.y;
	zc = obj3.z - obj2.z;

	if (yc == 0)
		w1 = 0.0;
	else if (xc >= 0 && yc >= 0)
		w1 = acos(xc / sqrt(xc * xc + yc * yc));
	else if (xc <= 0 && yc >= 0)
		w1 = -acos(-xc / sqrt(xc * xc + yc * yc)) + PAI;
	else if (xc <= 0 && yc <= 0)
		w1 = acos(-xc / sqrt(xc * xc + yc * yc)) - PAI;
	else if (xc >= 0 && yc <= 0)
		w1 = -acos(xc / sqrt(xc * xc + yc * yc));

	/* rotating all the point clockwise about Z axis a w1 degree
	   then the atom C is lying in the Z-X plane */

	xa3 = xa * cos(w1) + ya * sin(w1);
	ya3 = -xa * sin(w1) + ya * cos(w1);
	za3 = za;
	xc3 = xc * cos(w1) + yc * sin(w1);
	/*
	   yc3 = -xc * sin(w1) + yc * cos(w1);
	 */
	zc3 = zc;

	/* rotating all  point clockwise about Y axis a w2 degree then
	   atom C is lying the X axis */

	if (zc3 >= 0 && xc3 >= 0)
		w2 = -acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 >= 0)
		w2 = acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 <= 0)
		w2 = -acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) + PAI;
	else if (zc3 >= 0 && xc3 <= 0)
		w2 = acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) - PAI;
	xa4 = -za3 * sin(w2) + xa3 * cos(w2);
	za4 = za3 * cos(w2) + xa3 * sin(w2);
	ya4 = ya3;
	xc4 = -zc3 * sin(w2) + xc3 * cos(w2);
	/*
	   yc4 = yc3;
	   zc4 = zc3 * cos(w2) + xc3 * sin(w2);
	 */
	/* The following program  rotating A,D clockwise about X axis a w3 degree
	   then the atom A is lying in the X-Y plane */

	if (ya4 >= 0 && za4 <= 0)
		w3 = -acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 >= 0 && za4 >= 0)
		w3 = acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 <= 0 && za4 >= 0)
		w3 = -acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) + PAI;
	else if (ya4 <= 0 && za4 <= 0)
		w3 = acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) - PAI;

	xa5 = xa4;
	ya5 = ya4 * cos(w3) + za4 * sin(w3);

	/*
	   za5 = -ya4 * sin(w3) + za4 * cos(w3);
	 */

	xc5 = xc4;

	/*
	   yc5 = yc4 * cos(w3) + zc4 * sin(w3);
	   zc5 = -yc4 * sin(w3) + zc4 * cos(w3);
	 */

	/* calculate the angle C-B-A */

	w4 = 0.0;
	if (xa5 * xc5 >= 0) {
		if (xa5 >= 0)
			w4 = acos(xa5 / sqrt(ya5 * ya5 + xa5 * xa5));
		else
			w4 = -acos(-xa5 / sqrt(ya5 * ya5 + xa5 * xa5));
	} else {
		if (xa5 >= 0) {
			w4 = PAI - acos(xa5 / sqrt(ya5 * ya5 + xa5 * xa5));
		} else
			w4 = -PAI + acos(-xa5 / sqrt(ya5 * ya5 + xa5 * xa5));
	}
	*ww1 = w1;
	*ww2 = w2;
	*ww3 = w3;
	return w4 * 180.000 / PAI;
}


double distance0(ATOM obj1, ATOM obj2, double *ww1, double *ww2)
{
	double w1, w2;
	double xa, ya, za;
	double xa3, za3;
	/* double ya3; */
	double xa4;
	/* double ya4, za4; */

	xa = obj1.x - obj2.x;
	ya = obj1.y - obj2.y;
	za = obj1.z - obj2.z;

	if (ya == 0)
		w1 = 0.0;
	else if (xa >= 0 && ya >= 0)
		w1 = acos(xa / sqrt(xa * xa + ya * ya));
	else if (xa <= 0 && ya >= 0)
		w1 = -acos(-xa / sqrt(xa * xa + ya * ya)) + PAI;
	else if (xa <= 0 && ya <= 0)
		w1 = acos(-xa / sqrt(xa * xa + ya * ya)) - PAI;
	else if (xa >= 0 && ya <= 0)
		w1 = -acos(xa / sqrt(xa * xa + ya * ya));

	/* rotating all the point clockwise about Z axis a w1 degree
	   then the atom A is lying in the Z-X plane */

	xa3 = xa * cos(w1) + ya * sin(w1);
	/*
	   ya3 = -xa * sin(w1) + ya * cos(w1);
	 */
	za3 = za;

	/* rotating all  point clockwise about Y axis a w2 degree then
	   atom A is lying the X axis */

	if (za3 >= 0 && xa3 >= 0)
		w2 = -acos(xa3 / sqrt(xa3 * xa3 + za3 * za3));
	else if (za3 <= 0 && xa3 >= 0)
		w2 = acos(xa3 / sqrt(xa3 * xa3 + za3 * za3));
	else if (za3 <= 0 && xa3 <= 0)
		w2 = -acos(-xa3 / sqrt(xa3 * xa3 + za3 * za3)) + PAI;
	else if (za3 >= 0 && xa3 <= 0)
		w2 = acos(-xa3 / sqrt(xa3 * xa3 + za3 * za3)) - PAI;
	xa4 = -za3 * sin(w2) + xa3 * cos(w2);
	/*
	   za4 = za3 * cos(w2) + xa3 * sin(w2);
	   ya4 = ya3;
	 */

	/* calculate the bondlength */
	*ww1 = w1;
	*ww2 = w2;
	return xa4;
}




double rotate00(ATOM obj1, ATOM obj2, ATOM obj3, ATOM * obj4,
				double twist, int AngleChange)
{

/*
This function is used to calculate the coordinate of the fourth atom in
A-B-C-D system when we perform rotation bond B-C a certain degree
The function can calculate the twist angle of four atoms before rotating
and if you want twist angle be changed, please give you value(twist), the
varible AngleChange should larger than 0
*/

	double w1, w2, w3, t;
	double ya5, xc5, yd5, zd5;
	double xa, ya, za, xc, yc, zc, xd, yd, zd;
	double xa3, ya3, za3, xc3, zc3, xd3, yd3, zd3;
	double ya4, za4, xc4, xd4, yd4, zd4;
	double xx, yy, zz, xxx, yyy, zzz;
	twist = twist / 180 * PAI;
	xa = obj1.x - obj2.x;
	ya = obj1.y - obj2.y;
	za = obj1.z - obj2.z;
	xc = obj3.x - obj2.x;
	yc = obj3.y - obj2.y;
	zc = obj3.z - obj2.z;
	xd = (*obj4).x - obj2.x;
	yd = (*obj4).y - obj2.y;
	zd = (*obj4).z - obj2.z;
	if (yc == 0)
		w1 = 0.0;
	else if (xc >= 0 && yc >= 0)
		w1 = acos(xc / sqrt(xc * xc + yc * yc));
	else if (xc <= 0 && yc >= 0)
		w1 = -acos(-xc / sqrt(xc * xc + yc * yc)) + PAI;
	else if (xc <= 0 && yc <= 0)
		w1 = acos(-xc / sqrt(xc * xc + yc * yc)) - PAI;
	else if (xc >= 0 && yc <= 0)
		w1 = -acos(xc / sqrt(xc * xc + yc * yc));

	/* rotating all the point clockwise about Z axis a w1 degree
	   then the atom C is lying in the Z-X plane */

	xa3 = xa * cos(w1) + ya * sin(w1);
	ya3 = -xa * sin(w1) + ya * cos(w1);
	za3 = za;
	xc3 = xc * cos(w1) + yc * sin(w1);

	/*  yc3=-xc*sin(w1)+yc*cos(w1);
	   this value is not use anymore */

	zc3 = zc;
	xd3 = xd * cos(w1) + yd * sin(w1);
	yd3 = -xd * sin(w1) + yd * cos(w1);
	zd3 = zd;

	/* rotating all  point clockwise about Y axis a w2 degree then
	   atom C is lying the X axis */

	if (zc3 >= 0 && xc3 >= 0)
		w2 = -acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 >= 0)
		w2 = acos(xc3 / sqrt(xc3 * xc3 + zc3 * zc3));
	else if (zc3 <= 0 && xc3 <= 0)
		w2 = -acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) + PAI;
	else if (zc3 >= 0 && xc3 <= 0)
		w2 = acos(-xc3 / sqrt(xc3 * xc3 + zc3 * zc3)) - PAI;
	za4 = za3 * cos(w2) + xa3 * sin(w2);
	ya4 = ya3;
	xc4 = -zc3 * sin(w2) + xc3 * cos(w2);
	zd4 = zd3 * cos(w2) + xd3 * sin(w2);
	yd4 = yd3;
	xd4 = -zd3 * sin(w2) + xd3 * cos(w2);

	/* The following program  rotating A,D clockwise about X axis a w3 degree
	   then the atom A is lying in the X-Y plane */

	if (ya4 >= 0 && za4 <= 0)
		w3 = -acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 >= 0 && za4 >= 0)
		w3 = acos(ya4 / sqrt(ya4 * ya4 + za4 * za4));
	else if (ya4 <= 0 && za4 >= 0)
		w3 = -acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) + PAI;
	else if (ya4 <= 0 && za4 <= 0)
		w3 = acos(-ya4 / sqrt(ya4 * ya4 + za4 * za4)) - PAI;

	/*  xa5=xa4; */

	ya5 = ya4 * cos(w3) + za4 * sin(w3);
	xc5 = xc4;
	yd5 = yd4 * cos(w3) + zd4 * sin(w3);
	zd5 = -yd4 * sin(w3) + zd4 * cos(w3);

	if (xc5 > 0) {
		if (ya5 == 0)
			ya5 += 0.0000000001;
		if (yd5 == 0)
			yd5 += 0.0000000001;
		if ((ya5 > 0 && yd5 > 0) || (ya5 < 0 && yd5 < 0)) {
			/* A and D are on the same side of B-C */
			if (zd5 == 0)
				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = PAI - acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
		} else {

			/* A, D are not on the same side of B-C bonds, so the initial 
			   degree is PAI or -3.1415926 */
			if (zd5 == 0)

				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
			if (zd5 > 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = -acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) + PAI;
		}
	} else {
		if (xc5 == 0)
			xc5 += 0.0000000001;
		if (yd5 == 0)
			yd5 += 0.0000000001;
		if ((ya5 > 0 && yd5 > 0) || (ya5 < 0 && yd5 < 0)) {
			/* A and D are on the same side of B-C */
			if (zd5 == 0)
				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = PAI - acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
		} else {

			/* A, D are not on the same side of bond B-C */

			if (zd5 == 0)
				zd5 += 0.0000000001;
			if (yd5 == 0)
				yd5 += 0.0000000001;
			if (zd5 < 0 && yd5 > 0)
				t = -acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 < 0 && yd5 < 0)
				t = -acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) + PAI;
			if (zd5 > 0 && yd5 > 0)
				t = acos(yd5 / sqrt(zd5 * zd5 + yd5 * yd5));
			if (zd5 > 0 && yd5 < 0)
				t = acos(-yd5 / sqrt(zd5 * zd5 + yd5 * yd5)) - PAI;
		}
	}
	t = t * 180 / PAI;
	/* for test purpose */

	t = t * PAI / 180;
	if (!AngleChange)
		twist = -t + twist;

	/* The following program restore the coordinate of atom D to
	   the original coordinate system which mean the first three atoms'
	   coordinates is unchangable */

	xx = xd4;
	yy = yd5 * cos(twist) - zd5 * sin(twist);
	zz = zd5 * cos(twist) + yd5 * sin(twist);
	yyy = yy * cos(w3) - zz * sin(w3);
	zzz = yy * sin(w3) + zz * cos(w3);
	xxx = xx;
	zz = zzz * cos(w2) - xxx * sin(w2);
	yy = yyy;
	xx = zzz * sin(w2) + xxx * cos(w2);
	xxx = xx * cos(w1) + yy * sin(-w1) + obj2.x;
	yyy = xx * sin(w1) + yy * cos(w1) + obj2.y;
	zzz = zz + obj2.z;
	(*obj4).x = xxx;
	(*obj4).y = yyy;
	(*obj4).z = zzz;
	return twist * 180 / PAI;

}



void intercoord(int atomnum, ATOM atom[])
{
	int i, j, k, m;
	int tmpint1, tmpint2, tmpint3;
	int *select;
	int breakindex;

	select = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (select == NULL) {
                fprintf(stdout, "memory allocation error for *index in intercoord()\n");
                exit(1);
        }
	
	for (i = 0; i < atomnum; i++)
		select[i] = -1;

	atom[0].bond = 0.0;
	atom[0].angle = 0.0;
	atom[0].twist = 0.0;
	atom[0].bondatom = -1;
	atom[0].angleatom = -1;
	atom[0].twistatom = -1;
	atom[1].bondatom = 0;
	atom[1].bond = distance(atom[0], atom[1]);
	atom[1].angleatom = -1;
	atom[1].angle = 0.0;
	atom[1].twistatom = -1;
	atom[1].twist = 0.0;
	breakindex = 0;
	for (i = 0; i < 6; i++)
		if (atom[2].con[i] == 1) {
			atom[2].bondatom = 1;
			atom[2].bond = distance(atom[1], atom[2]);
			atom[2].angleatom = 0;
			atom[2].angle = anglecal(atom[2], atom[1], atom[0]);
			breakindex = 1;
			break;
		}
	if (breakindex == 0) {
		atom[2].bondatom = 0;
		atom[2].bond = distance(atom[0], atom[2]);
		atom[2].angleatom = 1;
		atom[2].angle = anglecal(atom[2], atom[0], atom[1]);
	}
	atom[2].twistatom = -1;
	atom[2].twist = 0.0;
	select[0] = 1;
	select[1] = 1;
	select[2] = 1;
	select[3] = 1;
	for (i = 3; i < atomnum; i++) {
		tmpint1 = -1;
		tmpint2 = -1;
		tmpint3 = -1;
		breakindex = 0;
		for (j = 0; j < 6; j++)
			if (select[atom[i].con[j]] == 1) {
				atom[i].bondatom = atom[i].con[j];
				tmpint1 = atom[i].bondatom;
				if (tmpint1 == -1)
					continue;
				for (k = 0; k < 6; k++)
					if (select[atom[tmpint1].con[k]] == 1
						&& atom[tmpint1].con[k] != i) {
						atom[i].angleatom = atom[tmpint1].con[k];
						tmpint2 = atom[i].angleatom;
						if (tmpint2 == -1)
							continue;
						for (m = 0; m < 6; m++)
							if (select[atom[tmpint2].con[m]] == 1
								&& atom[tmpint2].con[m] != tmpint1) {
								atom[i].twistatom = atom[tmpint2].con[m];
								tmpint3 = atom[i].twistatom;
								if (tmpint3 == -1)
									continue;
								if (tmpint3 != -1 && tmpint2 != -1
									&& tmpint3 != -1) {
									breakindex = 1;
									break;
								}
							}
						if (breakindex == 1)
							break;
					}
				if (breakindex == 1)
					break;
			}
		if (breakindex == 0)	/*improper torsional angle */
			for (j = 0; j < 6; j++) {
				if( tmpint1 == -1 ) continue; 
					/* above fix suggested by Charles Karney;
					   probably should emit an error instead?  */
				if (select[atom[tmpint1].con[j]] == 1
					&& atom[tmpint1].con[j] != atom[i].angleatom)
					atom[i].twistatom = atom[tmpint1].con[j];
				break;
			}
		select[i] = 1;
		atom[i].bond = distance(atom[i], atom[atom[i].bondatom]);
		atom[i].angle = anglecal(atom[i], atom[atom[i].bondatom],
								 atom[atom[i].angleatom]);
		atom[i].twist = rotate(atom[i], atom[atom[i].bondatom],
							   atom[atom[i].angleatom],
							   &atom[atom[i].twistatom]);
/*		      printf("\n%5s%5d%5d%5d%8.3lf%8.3lf%8.3lf",
		      atom[i].name, atom[i].bondatom, atom[i].angleatom,
		      atom[i].twistatom, atom[i].bond, atom[i].angle,
		      atom[i].twist); */
	}
}

void cartcoord(int atomnum, ATOM atom[])
{
	int i;
	atom[0].x = 0.0;
	atom[0].y = 0.0;
	atom[0].z = 0.0;
	atom[1].x = atom[1].bond;
	atom[1].y = 0.0;
	atom[1].z = 0.0;
	if (atom[2].bondatom == 0)
		atom[2].x = atom[0].x + atom[2].bond * cos(DEGRAD * atom[2].angle);
	else
		atom[2].x = atom[1].x - atom[2].bond * cos(DEGRAD * atom[2].angle);
	atom[2].y = atom[2].bond * sin(DEGRAD * atom[2].angle);
	atom[2].z = 0.0;

	for (i = 3; i < atomnum; i++) {
/*
 printf("\n%5d%5s%5d %8.3lf%5d %8.3lf%5d %8.3lf", i, atom[i].name, atom[i].bondatom, atom[i].bond,
          atom[i].angleatom, atom[i].angle, atom[i].twistatom, atom[i].twist);
*/
		rotate0(atom[atom[i].twistatom], atom[atom[i].angleatom],
				atom[atom[i].bondatom], &atom[i], atom[i].bond,
				atom[i].angle, atom[i].twist);
	}
}



int connect(char *connect_file, int atomnum, ATOM * atom, int *bondnum,
			BOND * bond, int maxbond)
{
	int i, j;
	int tmpint;
	int numbond = 0;
	int overflow_flag = 0;
	char line[MAXCHAR];
	FILE *fpin;
	double tmpfloat;
	double offset;
	double ref;
	double connectradius[120];
	if ((fpin = fopen(connect_file, "r")) == NULL) {
		fprintf(stdout, "Cannot open connect_file %s in connect(), exit\n", connect_file);
		exit(1);
	}
	for (;;) {
		if (fgets(line, LINELEN_MAX, fpin) == NULL)
			break;
		if (strncmp(".", &line[10], 1) == 0)
			sscanf(&line[3], "%d%lf", &tmpint, &tmpfloat);
		connectradius[tmpint] = tmpfloat;
	}
	for (i = 0; i < atomnum; i++) {
		atom[i].connum = 0;
		atom[i].con[0] = -1;
		atom[i].con[1] = -1;
		atom[i].con[2] = -1;
		atom[i].con[3] = -1;
		atom[i].con[4] = -1;
		atom[i].con[5] = -1;
	}
	for (i = 0; i < atomnum - 1; i++)
		for (j = i + 1; j < atomnum; j++)
			if (((atom[i].resno - atom[j].resno) == 1 && atom[j].ter != 1)
				|| ((atom[j].resno - atom[i].resno) == 1 && atom[i].ter != 1)
				|| (atom[i].resno == atom[j].resno)) {

				ref =
					connectradius[atom[i].atomicnum] +
					connectradius[atom[j].atomicnum];
				if (ref <= 1.5)
					offset = ref * 0.15;
				if (ref > 1.5 && ref <= 1.90)
					offset = ref * 0.11;
				if (ref > 1.90 && ref <= 2.05)
					offset = ref * 0.09;
				if (ref > 2.05)
					offset = ref * 0.08;
				tmpfloat = distance(atom[i], atom[j]);
				if (tmpfloat < (ref + offset) && tmpfloat > ref * 0.5)
					if (atom[i].connum < 6 && atom[j].connum < 6) {
						atom[i].con[atom[i].connum++] = j;
						atom[j].con[atom[j].connum++] = i;
						if (overflow_flag == 0) {
							bond[numbond].bondi = i;
							bond[numbond].bondj = j;
							bond[numbond].type = 0;
						}
						numbond++;
						if (numbond >= maxbond && overflow_flag == 0) {
							printf
								("\nInfo: the bond number exceeds MAXBOND, reallocate memory automatically\n");
							overflow_flag = 1;
						}
					}
			}
	*bondnum = numbond;
	fclose(fpin);
	return *bondnum;
}

void omegarotate(ATOM obj1, ATOM obj2, double *w1, double *w2)
{
        double x, y, z;
        double x2, z2;
/*
double y2;
*/

        x = (obj1).x - obj2.x;
        y = (obj1).y - obj2.y;
        z = (obj1).z - obj2.z;

        if (y == 0)
                *w1 = 0.0;
        else if (x >= 0 && y >= 0)
                *w1 = acos(x / sqrt(x * x + y * y));
        else if (x <= 0 && y >= 0)
                *w1 = -acos(-x / sqrt(x * x + y * y)) + PAI;
        else if (x <= 0 && y <= 0)
                *w1 = acos(-x / sqrt(x * x + y * y)) - PAI;
        else if (x >= 0 && y <= 0)
                *w1 = -acos(x / sqrt(x * x + y * y));

        /* rotating all the point clockwise about Z axis a w1 degree
           then the atom C is lying in the Z-X plane */

        x2 = x * cos(*w1) + y * sin(*w1);
/*
        y2 = -x * sin(*w1) + y * cos(*w1);
*/
        z2 = z;

        /* rotating all  point clockwise about Y axis a w2 degree then
           atom C is lying the X axis */

        if (z2 >= 0 && x2 >= 0)
                *w2 = -acos(x2 / sqrt(x2 * x2 + z2 * z2));
        else if (z2 <= 0 && x2 >= 0)
                *w2 = acos(x2 / sqrt(x2 * x2 + z2 * z2));
        else if (z2 <= 0 && x2 <= 0)
                *w2 = -acos(-x2 / sqrt(x2 * x2 + z2 * z2)) + PAI;
        else if (z2 >= 0 && x2 <= 0)
                *w2 = acos(-x2 / sqrt(x2 * x2 + z2 * z2)) - PAI;

}

void omegarotate2(ATOM atom[], int atomnum, ATOM refatom, double angle, double w1, double w2) {
        double x,y,z;
        double x2, y2, z2;
        double x3, y3, z3;
        double x4, y4, z4;
        double x5, y5, z5;
        double x6, y6, z6;
        int i;
        angle = angle / 180.0 * PAI;
        for(i = 0; i<atomnum; i++) {
                x = atom[i].x - refatom.x;
                y = atom[i].y - refatom.y;
                z = atom[i].z - refatom.z;

                x2 = x * cos(w1) + y * sin(w1);
                y2 = -x * sin(w1) + y * cos(w1);
                z2 = z;

                x3 = -z2 * sin(w2) + x2 * cos(w2);
                y3 = y2;
                z3 = z2 * cos(w2) + x2 * sin(w2);

                x4 = x3;
                y4 = y3 * cos(angle) + z3 * sin(angle);
                z4 = -y3 * sin(angle) + z3 * cos(angle);

                x5 = z4 * sin(w2) + x4 * cos(w2);
                y5 = y4;
                z5 = z4 * cos(w2) - x4 * sin(w2);

                x6 = x5 * cos(w1) - y5 * sin(w1) + refatom.x;
                y6 = x5 * sin(w1) + y5 * cos(w1) + refatom.y;
                z6 = z5 + refatom.z;
                atom[i].x =x6;
                atom[i].y =y6;
                atom[i].z =z6;
        }
}

void alignx(ATOM obj1, ATOM obj2, ATOM atom[], int atomnum)
{
	int i;
        double x, y, z;
        double x2, y2, z2;
	double coordx, coordy, coordz;
	double w1, w2;
/* assume the original is the center of points obj1 and obj2 */
        x = 0.5*((obj1).x - obj2.x);
        y = 0.5*((obj1).y - obj2.y);
        z = 0.5*((obj1).z - obj2.z);

	for(i=0;i<atomnum;i++) {
		atom[i].x -= obj2.x;			
		atom[i].y -= obj2.y;			
		atom[i].z -= obj2.z;			
	}

        /* rotating Obj1 clockwise about Z axis a w1 degree
           then the atom Obj1 is lying in the Z-X plane */

        if (y == 0)
                w1 = 0.0;
        else if (x >= 0 && y >= 0)
                w1 = acos(x / sqrt(x * x + y * y));
        else if (x <= 0 && y >= 0)
                w1 = -acos(-x / sqrt(x * x + y * y)) + PAI;
        else if (x <= 0 && y <= 0)
                w1 = acos(-x / sqrt(x * x + y * y)) - PAI;
        else if (x >= 0 && y <= 0)
                w1 = -acos(x / sqrt(x * x + y * y));

	for(i=0;i<atomnum;i++) {
		coordx = atom[i].x;
		coordy = atom[i].y;
		coordz = atom[i].z;
        	atom[i].x = coordx * cos(w1) + coordy * sin(w1);
        	atom[i].y = -coordx * sin(w1) + coordy * cos(w1);
	}
	x2 = x*cos(w1) + y*sin(w1);
	y2 = -x*sin(w1) + y*cos(w1);
	z2 = z;

        /* rotating Obj1 clockwise about Y axis a w2 degree then
           Obj1 is lying in the X axis */
        if (z2 >= 0 && x2 >= 0)
                w2 = -acos(x2 / sqrt(x2 * x2 + z2 * z2));
        else if (z2 <= 0 && x2 >= 0)
                w2 = acos(x2 / sqrt(x2 * x2 + z2 * z2));
        else if (z2 <= 0 && x2 <= 0)
                w2 = -acos(-x2 / sqrt(x2 * x2 + z2 * z2)) + PAI;
        else if (z2 >= 0 && x2 <= 0)
                w2 = acos(-x2 / sqrt(x2 * x2 + z2 * z2)) - PAI;

	for(i=0;i<atomnum;i++) {
		coordx = atom[i].x;
		coordy = atom[i].y;
		coordz = atom[i].z;
        	atom[i].x = -coordz * sin(w2) + coordx * cos(w2);
        	atom[i].z =  coordz * cos(w2) + coordx * sin(w2);
	}
}


void aligny(ATOM obj1, ATOM obj2, ATOM atom[], int atomnum)
{
	int i;
        double x, y, z;
        double x2, y2, z2;
	double coordx, coordy, coordz;
	double w1, w2;

/* assume the original is the center of points obj1 and obj2 */
        x = 0.5*((obj1).x - obj2.x);
        y = 0.5*((obj1).y - obj2.y);
        z = 0.5*((obj1).z - obj2.z);

	for(i=0;i<atomnum;i++) {
		atom[i].x -= obj2.x;			
		atom[i].y -= obj2.y;			
		atom[i].z -= obj2.z;			
	}
        /* The following program rotates all points along X axis by w1
           degrees, then the atom C is lying in the X-Y plane */

        if (z == 0)
                w1 = 0.0;
        if (y >= 0 && z <= 0)
                w1 = -acos(y / sqrt(y * y + z * z));
        else if (y >= 0 && z >= 0)
                w1 = acos(y / sqrt(y * y + z * z));
        else if (y <= 0 && z >= 0)
                w1 = -acos(-y / sqrt(y * y + z * z)) + PAI;
        else if (y <= 0 && z <= 0)
                w1 = acos(-y / sqrt(y * y + z * z)) - PAI;

	for(i=0;i<atomnum;i++) {
		coordx = atom[i].x;
		coordy = atom[i].y;
		coordz = atom[i].z;
        	atom[i].y = coordy * cos(w1) + coordz * sin(w1);
        	atom[i].z = -coordy * sin(w1) + coordz * cos(w1);
	}
	x2 = x;
	y2 = y*cos(w1) + z*sin(w1);
	z2 = -y*sin(w1) + z*cos(w1);

        /* rotating Point1 clockwise about Z axis a w2 degree then
           it is lying the Y axis */
	if (y2 == 0) w2 = 0.0;
        if (x2 >= 0 && y2 >= 0)
                w2 = acos(x2 / sqrt(x2 * x2 + y2 * y2)) - HALFPAI;
        else if (x2 <= 0 && y2 >= 0)
                w2 = HALFPAI -acos(-x2 / sqrt(x2 * x2 + y2 * y2)) ;
        else if (y2 <= 0 && x2 <= 0)
                w2 = acos(-x2 / sqrt(x2 * x2 + y2 * y2)) + HALFPAI;
        else if (y2 <= 0 && x2 >= 0)
                w2 = -acos(x2 / sqrt(x2 * x2 + y2 * y2)) - HALFPAI;

	for(i=0;i<atomnum;i++) {
		coordx = atom[i].x;
		coordy = atom[i].y;
		coordz = atom[i].z;
        	atom[i].x =  coordx * cos(w2) + coordy * sin(w2);
        	atom[i].y = -coordx * sin(w2) + coordy * cos(w2);
	}
}


void alignz(ATOM obj1, ATOM obj2, ATOM atom[], int atomnum)
{
	int i;
        double x, y, z;
        double x2,y2, z2;
	double coordx, coordy, coordz;
	double w1, w2;

/* assume the original is the center of points obj1 and obj2 */
        x = 0.5*((obj1).x - obj2.x);
        y = 0.5*((obj1).y - obj2.y);
        z = 0.5*((obj1).z - obj2.z);

	for(i=0;i<atomnum;i++) {
		atom[i].x -= obj2.x;			
		atom[i].y -= obj2.y;			
		atom[i].z -= obj2.z;			
	}

        if (y == 0)
                w1 = 0.0;
        else if (x >= 0 && y >= 0)
                w1 = acos(x / sqrt(x * x + y * y));
        else if (x <= 0 && y >= 0)
                w1 = -acos(-x / sqrt(x * x + y * y)) + PAI;
        else if (x <= 0 && y <= 0)
                w1 = acos(-x / sqrt(x * x + y * y)) - PAI;
        else if (x >= 0 && y <= 0)
                w1 = -acos(x / sqrt(x * x + y * y));

        /* rotating Obj1 clockwise about Z axis a w1 degree
           then it is lying in the Z-X plane */

	for(i=0;i<atomnum;i++) {
		coordx = atom[i].x;
		coordy = atom[i].y;
		coordz = atom[i].z;
        	atom[i].x = coordx * cos(w1) + coordy * sin(w1);
        	atom[i].y = -coordx * sin(w1) + coordy * cos(w1);
	}

        x2 = x*cos(w1) + y*sin(w1);
        y2 = -x*sin(w1) + y*cos(w1);
        z2 = z; 

        /* rotating Obj1 clockwise about Y axis a w2 degree then
           it is lying in the Z axis */
	if(x2 == 0 && z2 ==0 ) w2 = 0;
        if (z2 >= 0 && x2 >= 0)
                w2 = acos(z2 / sqrt(x2 * x2 + z2 * z2));
        else if (z2 <= 0 && x2 >= 0)
                w2 = -acos(-z2 / sqrt(x2 * x2 + z2 * z2)) + PAI;
        else if (z2 <= 0 && x2 <= 0)
                w2 = acos(-z2 / sqrt(x2 * x2 + z2 * z2)) - PAI;
        else if (z2 >= 0 && x2 <= 0)
                w2 = -acos(z2 / sqrt(x2 * x2 + z2 * z2)) ;

	for(i=0;i<atomnum;i++) {
		coordx = atom[i].x;
		coordy = atom[i].y;
		coordz = atom[i].z;
        	atom[i].x =  -coordz * sin(w2) + coordx * cos(w2);
        	atom[i].z =  coordz * cos(w2) + coordx * sin(w2);
	}
}

/* find all atoms linking from atid1 to atid2 */
int group(int atomnum, ATOM atom[], int atid1, int atid2) {
int i,j;
int satomnum = 0;
int rflag;
for(i=0;i<atomnum;i++)
        atom[i].select = 0;
atom[atid2].select = 1;
satomnum = 1;
rflag = 1;

while(rflag == 1) {
        rflag = 0;
        for(i=0;i<atomnum;i++) {
                if(atom[i].select == 0) continue;
                for(j=0;j<atom[i].connum;j++)
                        if(atom[atom[i].con[j]].select == 0) {
                                if(atom[i].con[j] != atid1) {
                                        atom[atom[i].con[j]].select = 1;
                                        satomnum ++;
                                        rflag = 1;
                                }
                                else {
                                        if(i == atid2)
                                                continue;
                                        else
                                                return -1;
/* ring close happens*/
                                }
                        }
        }
}
return satomnum;
}

