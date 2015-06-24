#include "../include/dprec.fh"
!-------------BEGIN    md.h  ------------------------------------------------

! NOTE: if you change this file, make sure you change the corresponding file
! in both src/sander/md.h, AmberTools/src/sqm/md.h, and AmberTools/src/pbsa/md.h

integer BC_MDI  ! size in integers of common block mdi
integer BC_MDR  ! size in Reals of common block mdr

! ... integer variables:

integer nrp,nspm,ig,ntx,ntcx,            &!5
      ntxo,ntt,ntp,ntr,init,             &!10
      ntcm,nscm,nsolut,klambda,          &!14
      ntc,ntcc,ntf,ntid,ntn,             &!19
      ntnb,nsnb,ndfmin,nstlim,nrc,       &!24
      ntrx,npscal,imin,maxcyc,ncyc,      &!29
      ntmin,irest,jfastw,barostat,       &!33
      ibgwat,ienwat,iorwat,mcbarint,     &!37
      iwatpr,nsolw,igb,alpb,iyammp,      &!42
      gbsa,vrand,iwrap,nrespa,irespa,nrespai,icfe,       &!49
      rbornstat,ivcap,iconstreff,                        &!52
      neb,vv,tmode,ipol,iesp,ievb,nodeid,num_noshake,    &!60
      idecomp,icnstph,ntcnstph,maxdup,numexchg,repcrd,numwatkeep,hybridgb, &!68
      ibgion,ienion,profile_mpi,lj1264,                  &!72
      ipb,inp,ntrelax,relaxing,dec_verbose,vdwmodel,     &!78
      csurften,ninterface,no_ntt3_sync,nkija,idistr,     &!83
      nucat !84

common/mdi/nrp,nspm,ig, &                                               !3
      ntx,ntcx,ntxo,ntt,ntp,ntr,init,ntcm,nscm, &                       !12
      nsolut,ntc,ntcc,ntf,ntid,ntn,ntnb,nsnb,ndfmin, &                  !21
      nstlim,nrc,ntrx,npscal,imin,maxcyc,ncyc,ntmin, &                  !29
      irest,jfastw,ibgwat,ienwat,iorwat, &                              !34
      iwatpr,nsolw,igb,alpb,iyammp,gbsa,vrand,numexchg,repcrd, &        !43
      numwatkeep,hybridgb,barostat,mcbarint, &                          !47
      iwrap,nrespa,irespa,nrespai,icfe,rbornstat, &                     !53
      ivcap,iconstreff,idecomp,klambda,icnstph,ntcnstph,maxdup,neb,vv, &!62
      tmode,ipol,iesp,ievb,nodeid,num_noshake,ibgion,ienion, &          !70
      profile_mpi,lj1264,ipb,inp,ntrelax,relaxing,dec_verbose,vdwmodel, &!78
      csurften,ninterface,no_ntt3_sync,nkija,idistr,     &!83
      nucat !84

parameter (BC_MDI=84) ! Number of elements in the common block;
                      ! Be sure to update if you change things

! ... floats:

_REAL_ t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, & !9
      tol,taur,dx0,drms,vlimit,rbtarg(9),tmass,tmassinv,  & !25
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt,  & !32
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon,  & !38
      solvph,rgbmax,fsmax,restraint_wt, &  !42
      skmin,skmax,vfac,gbneckscale,v11,v12,v22,kevb,evbt,Arad, & !52
      gbalphaH,gbbetaH,gbgammaH, & !55 igb 8 pars
      gbalphaC,gbbetaC,gbgammaC, & !58
      gbalphaN,gbbetaN,gbgammaN, & !61
      gbalphaOS,gbbetaOS,gbgammaOS, & !64 
      gbalphaP,gbbetaP,gbgammaP, &    !67
      Sh,Sc,Sn,So,Ss,Sp, &  ! 73 end igb 8 pars
      gamma_ten, & !74
      gb_alpha_hnu, gb_beta_hnu, gb_gamma_hnu, & !77 ! GBneck2nu pars
      gb_alpha_cnu, gb_beta_cnu, gb_gamma_cnu, & !80
      gb_alpha_nnu, gb_beta_nnu, gb_gamma_nnu, & !83
      gb_alpha_osnu, gb_beta_osnu, gb_gamma_osnu, & !86
      gb_alpha_pnu, gb_beta_pnu, gb_gamma_pnu, &    !89
      screen_hnu, screen_cnu, screen_nnu, screen_onu, & !93
      screen_pnu    !94 ! end GBneck2nu pars

common/mdr/t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, &             !9
      tol,taur,dx0,drms,vlimit,rbtarg,tmass,tmassinv, &               !25
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt, &            !32
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon, &             !38
      solvph,rgbmax,fsmax,restraint_wt,skmin,skmax,vfac,gbneckscale, &!46
      v11,v12,v22,kevb,evbt,Arad, & !52
      gbalphaH,gbbetaH,gbgammaH,   & !55    !igb 8 pars
      gbalphaC,gbbetaC,gbgammaC, &   !58
      gbalphaN,gbbetaN,gbgammaN, &   !61
      gbalphaOS,gbbetaOS,gbgammaOS, & !64
      gbalphaP,gbbetaP,gbgammaP, & !67
      Sh,Sc,Sn,So,Ss,Sp, & !73 end igb 8 pars
      gamma_ten, & !74
      gb_alpha_hnu, gb_beta_hnu, gb_gamma_hnu, & !77 ! GBneck2nu pars
      gb_alpha_cnu, gb_beta_cnu, gb_gamma_cnu, & !80
      gb_alpha_nnu, gb_beta_nnu, gb_gamma_nnu, & !83
      gb_alpha_osnu, gb_beta_osnu, gb_gamma_osnu, & !86
      gb_alpha_pnu, gb_beta_pnu, gb_gamma_pnu, &    !89
      screen_hnu, screen_cnu, screen_nnu, screen_onu, & !93
      screen_pnu    !94 ! end GBneck2nu pars

parameter (BC_MDR=94) ! Number of elements in the common block;
                      ! Be sure to update if you change things

! ... strings:

character(len=4) iwtnm,iowtnm,ihwtnm
character(len=256) restraintmask,bellymask,tgtfitmask,&
            tgtrmsmask,noshakemask,crgmask,iwrap_mask
common/mds/ restraintmask,bellymask,tgtfitmask,tgtrmsmask,noshakemask,crgmask,  &
            iwtnm,iowtnm,ihwtnm(2),iwrap_mask

!-------------END    md.h  ------------------------------------------------
