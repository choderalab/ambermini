
!  --- following should not need to be modified unless you are adding
!      more variables to the "locmem" style of memory management

! NOTE: if you change this file, make sure you change the corresponding file
! in both src/sander/memory.h and AmberTools/src/sqm/memory.h

! MFC changing indices into IX and X to more rational names
!       I12 = Iibh
!       I14 = Ijbh
!       I16 = Iicbh
!       I18 = Iiba
!       I20 = Ijba
!       I22 = Iicba

integer       natom,nres,nbonh,nbona,ntheth,ntheta,nphih, &
      nphia,nnb,ntypes,nconp,maxmem,nwdvar,nparm, &
      natc,nattgtfit,nattgtrms,ibelly,natbel,ishake,nmxrs, &
      mxsub,natyp,npdec,i02,i04,i06,i08,i10, &
      iibh,ijbh,iicbh,iiba,ijba,iicba, &
      i24,i26,i28,i30,i32,i34,i36,i38,i40, &
      i42,i44,i46,i48,i50,i52,i54,i56,i58,ibellygp, &
      icnstrgp,itgtfitgp,itgtrmsgp,i64,i65,i66,i68, &
      i70,i72,i74,i76,i78,i79,i80,i82,i84,i86, &
      i100, &
      l15,lwinv,lpol,lcrd,lforce,l36,lvel,lvel2,l45,l50, &
      lcrdr,l60,l65,lmass,l75,l80,l85,l90,l95,l96,l97,l98,l99,lfrctmp, &
      l105,l110,l115,l120,l125,l130,l135,l140,l145,l150, &
      l165,l170,l175,l180,l185,l186,l187,l188,l189,l190, &
      m02,m04,m06,m08,m10,m12,m14,m16,m18,i01, &
      iifstwt,iifstwr,nrealb,nintb,nholb,npairb,lastr,lasti,lasth, &
      lastpr,nbper,ngper,ndper,ifpert,ncopy, &
      imask1,imask2,numadjst,mxadjmsk,noshake, &
! SEBOMD
      divchg  ,hessian  ,gradsebomd, grad1tmp ,grad2tmp, &
      grad3tmp, grad4tmp, &
! Modified by WJM, YD, RL
      l2402, l2403, l2404, ldf, lpol2, lpolbnd  ! Hai Nguyen: add GB index here

!  1        2         3         4         5      6     7     8      9      10
common/memory/ &
 natom   ,nres     ,nbonh    ,nbona   ,ntheth,ntheta,nphih  ,                       & !  7
 nphia   ,nnb      ,ntypes   ,nconp   ,maxmem,nwdvar,nparm  ,                       & ! 14
 natc    ,nattgtfit,nattgtrms,ibelly  ,natbel,ishake,nmxrs  ,                       & ! 21
 mxsub   ,natyp    ,npdec    ,i02     ,i04   ,i06   ,i08    ,i10   ,                & ! 29
 iibh    ,ijbh     ,iicbh    ,iiba    ,ijba  ,iicba ,                               & ! 35
 i24     ,i26      ,i28      ,i30     ,i32   ,i34   ,i36    ,i38   ,i40   ,         & ! 44
 i42     ,i44      ,i46      ,i48     ,i50   ,i52   ,i54    ,i56   ,i58   ,ibellygp,& ! 54
 icnstrgp,itgtfitgp,itgtrmsgp,i64     ,i65   ,i66   ,i68    ,                       & ! 61
 i70     ,i72      ,i74      ,i76     ,i78   ,i79   ,i80    ,i82   ,                & ! 69
 i84     ,i86      ,                                                                & ! 71
 i100    ,                                                                          & ! 72
 l15     ,lwinv    ,lpol     ,lcrd    ,lforce,l36   ,lvel   ,lvel2 ,                & ! 80
 l45     ,l50      ,                                                                & ! 82
 lcrdr   ,l60      ,l65      ,lmass   ,l75   ,l80   ,l85    ,l90   ,l95   ,l96     ,& ! 92
 l97     ,l98      ,l99      ,lfrctmp  ,                                            & ! 96
 l105    ,l110     ,l115     ,l120    ,l125  ,l130  ,l135   ,l140  ,l145  ,l150    ,& !106
 l165    ,l170     ,l175     ,l180    ,l185  ,l186  ,l187   ,l188  ,l189  ,l190    ,& !116
 m02     ,m04      ,m06      ,m08     ,m10   ,m12   ,m14    ,m16   ,m18   ,i01     ,& !126
 iifstwt ,iifstwr  ,nrealb   ,nintb   ,nholb ,npairb,lastr  ,lasti ,lasth ,         & !135
 lastpr  ,nbper    ,ngper    ,ndper   ,ifpert,ncopy ,                               & !141
 imask1  ,imask2   ,numadjst ,mxadjmsk,noshake,                                     & !146
 divchg  ,hessian  ,gradsebomd, grad1tmp ,grad2tmp, grad3tmp, grad4tmp,             & !153 SEBOMD
 l2402   ,l2403    ,l2404    ,ldf     ,lpol2 ,lpolbnd                                        !158
! Modified by WJM, RL
 !Hai Nguyen: add GB arrays here

! MEMORY_H_BC_MEMORY is the size of the MEMORY common block
#define MEMORY_H_BC_MEMORY 156
