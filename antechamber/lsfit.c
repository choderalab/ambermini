double lsfit(ATOM atm1[], ATOM atm2[], ATOM atm3[],int n) {
typedef struct {
double x;
double y;
double z;
} CORD;
int i,j,k;
int ict;
int ix,iy,iz;
int iflag;
/* ict: counter for number of iterations. must do at least 3 iterations.
ix, iy,iz: pointer used in iterative least squares. may have the value 1,2,3
the value changes on each iteration. */
double tol=0.00001;
double sig,gam, sg, bb,cc;
CORD *x0;
CORD *x1;
CORD *xyzfit;
double aa[3][3];
double rot[3][3];
double cg0[3];
double cg1[3];
double an;
double rms,ssum;

x0 = (CORD *) malloc(sizeof(CORD) * (n +10));
if (x0 == NULL) {
	fprintf(stdout, "memory allocation error for *name in lsfit()\n");
        exit(1);
}
x1 = (CORD *) malloc(sizeof(CORD) * (n +10));
if (x1 == NULL) {
	fprintf(stdout, "memory allocation error for *name in lsfit()\n");
        exit(1);
}
xyzfit = (CORD *) malloc(sizeof(CORD) * (n +10));
if (xyzfit == NULL) {
	fprintf(stdout, "memory allocation error for *name in lsfit()\n");
        exit(1);
}
ssum=0.;
for(i=0;i<n;i++) {
     	x1[i].x=atm2[i].x;
     	x1[i].y=atm2[i].y;
     	x1[i].z=atm2[i].z;
}
for(i=0;i<n;i++) {
    	x0[i].x=atm1[i].x;
    	x0[i].y=atm1[i].y;
    	x0[i].z=atm1[i].z;
}
an=1.0/(double)n;
for(i=0;i<3;i++) {
	cg0[i]=0.;
  	cg1[i]=0.;
}
  for(i=0;i<n;i++) {
	cg0[0] += x0[i].x;
	cg1[0] += x1[i].x;
	cg0[1] += x0[i].y;
	cg1[1] += x1[i].y;
	cg0[2] += x0[i].z;
	cg1[2] += x1[i].z;
}
for(i=0;i<3;i++) {
	cg0[i]*=an;
     	cg1[i]*=an;
}
for(i=0;i<3;i++)
	for(j=0;j<3;j++)
     		aa[i][j]=0.0;
for(k=0;k<n;k++) {
   aa[0][0]=aa[0][0]+(x1[k].x-cg1[0])*(x0[k].x-cg0[0]);
   aa[1][0]=aa[1][0]+(x1[k].y-cg1[1])*(x0[k].x-cg0[0]);
   aa[2][0]=aa[2][0]+(x1[k].z-cg1[2])*(x0[k].x-cg0[0]);
   aa[0][1]=aa[0][1]+(x1[k].x-cg1[0])*(x0[k].y-cg0[1]);
   aa[1][1]=aa[1][1]+(x1[k].y-cg1[1])*(x0[k].y-cg0[1]);
   aa[2][1]=aa[2][1]+(x1[k].z-cg1[2])*(x0[k].y-cg0[1]);
   aa[0][2]=aa[0][2]+(x1[k].x-cg1[0])*(x0[k].z-cg0[2]);
   aa[1][2]=aa[1][2]+(x1[k].y-cg1[1])*(x0[k].z-cg0[2]);
   aa[2][2]=aa[2][2]+(x1[k].z-cg1[2])*(x0[k].z-cg0[2]);
}
for(i=0;i<3;i++) {
    	for(j=0;j<3;j++)
      		rot[i][j]=0.;
    		rot[i][i]=1.0;
}
ict=0.;
goto a51;
a50:
ix++;
if(ix<4) goto a52;
if(iflag==0) goto a70;
a51:
iflag=0;
ix=1;
a52:
ict=ict+1;
if(ict>1000) goto a70;
iy=ix+1;
if(iy==4)  
	iy=1;
iz=6-ix-iy;
sig=aa[iz-1][iy-1]-aa[iy-1][iz-1];
gam=aa[iy-1][iy-1]+aa[iz-1][iz-1];
sg=sqrt(sig*sig+gam*gam);
if(sg==0) goto a50;
sg=1.0/sg;
if(fabs(sig)<tol*fabs(gam)) goto a50;
for(k=0;k<3;k++) {
	bb=gam*aa[iy-1][k]+sig*aa[iz-1][k];
   	cc=gam*aa[iz-1][k]-sig*aa[iy-1][k];
   	aa[iy-1][k]=bb*sg;
   	aa[iz-1][k]=cc*sg;
   	bb=gam*rot[iy-1][k]+sig*rot[iz-1][k];
   	cc=gam*rot[iz-1][k]-sig*rot[iy-1][k];
   	rot[iy-1][k]=bb*sg;
   	rot[iz-1][k]=cc*sg;
}
iflag=1;
goto a50;
a70:
 /* the following code translates the ligand the center of mass of the
 receptor site and rotates it according to the rotation matrxouneror ix
 calculated by orient */
for(i=0;i<n;i++) {
	xyzfit[i].x=cg0[0]+rot[0][0]*(x1[i].x-cg1[0])+ \
             rot[0][1]*(x1[i].y-cg1[1])+
             rot[0][2]*(x1[i].z-cg1[2]);
	xyzfit[i].y=cg0[1]+rot[1][0]*(x1[i].x-cg1[0])+ \
             rot[1][1]*(x1[i].y-cg1[1])+
             rot[1][2]*(x1[i].z-cg1[2]);
	xyzfit[i].z=cg0[2]+rot[2][0]*(x1[i].x-cg1[0])+ \
             rot[2][1]*(x1[i].y-cg1[1])+
             rot[2][2]*(x1[i].z-cg1[2]);
}
for(i=0;i<n;i++) {
	ssum+=(xyzfit[i].x-x0[i].x)*(xyzfit[i].x-x0[i].x);
	ssum+=(xyzfit[i].y-x0[i].y)*(xyzfit[i].y-x0[i].y);
	ssum+=(xyzfit[i].z-x0[i].z)*(xyzfit[i].z-x0[i].z);
}
rms=ssum/n;
rms=sqrt(rms);
for(i=0;i<n;i++) {
	atm3[i].x=xyzfit[i].x;
  	atm3[i].y=xyzfit[i].y;
  	atm3[i].z=xyzfit[i].z;
}
free(x0);
free(x1);
free(xyzfit);
return rms;
}

