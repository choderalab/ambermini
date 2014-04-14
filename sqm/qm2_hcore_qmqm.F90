! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
subroutine qm2_hcore_qmqm(COORD,H,W,ENUCLR)
!***********************************************************************/R
! Current code, optimisation and inlining by: Ross Walker (TSRI, 2005)
!
! d-orbital implementation by Taisung Lee (Rutgers, 2011)
!
! This routine is responsible for generating the one-electron matrix and
! the two electron integrals via calls to qm2_h1elec and qm2_rotate_qmqm.
! qm2_h1elec has been inlined in this code for speed.
!
! Current Version: Ross Walker (TSRI, 2005)
!
!IN -
! COORD = QM coordinates
!
!OUT-
! H = One electron matix
! W = Two electron integrals
! ENUCLR = Nuclear energy 
!***********************************************************************

 
      use constants, only : zero, one, A2_TO_BOHRS2, A_TO_BOHRS
      use ElementOrbitalIndex, only: MaxValenceOrbitals, MaxValenceDimension
      use qmmm_module, only : qmmm_nml, qmmm_struct, qm2_struct, qm2_params, qm2_rij_eqns, &
                              qmmm_mpi, OVERLAP_CUTOFF      
      use Rotation, only : GetRotationMatrix, Rotate2Center2Electron, RotateCore
!DEBUG      use utilitiesModule, only : Print
      
      implicit none

!Passed in
      _REAL_, intent(in) :: COORD(3,qmmm_struct%nquant_nlink)
      _REAL_, intent(out) :: W(qm2_struct%n2el)
      _REAL_, intent(out) :: ENUCLR
      _REAL_, intent(out) :: H(qm2_struct%matsize)

!Local
      _REAL_ ::SHMAT(MaxValenceOrbitals,MaxValenceOrbitals)
                 
      integer:: i, j, k, ni, i1, i2
      integer:: kr, j1, first_sj, last_pj, ii, j2, jj, ki
      integer:: loop_count, jstart, jend
      integer:: naoi, naoj, qmitype, qmjtype
      _REAL_:: enuc, elec_ke_p
      _REAL_:: half_num, r2, r2InAu, rij, rijInAu, oneOverRij
      _REAL_, allocatable :: WW(:,:)
      
      _REAL_ :: X(3),Y(3),Z(3), RI(22), CORE(10,2)
      _REAL_ :: rotationMatrix(15,45)
      
      integer:: firstIndexAO_i, firstIndexAO_j, lastIndexAO_i, lastIndexAO_j
      integer:: i_dimension, j_dimension
      logical::hasDOrbital

      W=zero
      enuclr=zero
      
! FILL THE DIAGONALS as we don't do them in the loop below.      
      do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
        ! everyone should s orbital
        if (qm2_params%natomic_orbs(i)>=1) then  
              
            firstIndexAO_i=qm2_params%orb_loc(1,i)
            H( qm2_params%pascal_tri2(firstIndexAO_i) )=qm2_params%orb_elec_ke(1,i)       
      
            ! p-orbitals
            if (qm2_params%natomic_orbs(i)>=4) then
       
                do j=firstIndexAO_i+1, firstIndexAO_i+3
                   H( qm2_params%pascal_tri2(j) )=qm2_params%orb_elec_ke(2,i) 
                end do
                
                ! d-orbitals      
                if (qm2_params%natomic_orbs(i)>=9) then
                    do j=firstIndexAO_i+4, firstIndexAO_i+8
                       H( qm2_params%pascal_tri2(j) )=qm2_params%orb_elec_ke(3,i) 
                    end do
                end if  ! d-orbitals
                  
             end if  ! p-orbitals 
        end if ! s-orbital                                                   
      end do


      loop_count=0
#ifdef MPI
      KR = qmmm_mpi%two_e_offset+1
      do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
         jstart =  qmmm_mpi%nquant_nlink_jrange(1,i)
         jend = qmmm_mpi%nquant_nlink_jrange(2,i)
#else
      KR=1                                                                      
      do I=2,qmmm_struct%nquant_nlink
         jstart = 1
         jend = i-1
#endif
 
         firstIndexAO_i=qm2_params%orb_loc(1,I)
         lastIndexAO_i=qm2_params%orb_loc(2,I)                
         naoi = qm2_params%natomic_orbs(i)
         qmitype = qmmm_struct%qm_atom_type(i) 
         NI=qmmm_struct%iqm_atomic_numbers(I)
!   FILL THE ATOM-OTHER ATOM ONE-ELECTRON MATRIX<PSI(LAMBDA)|PSI(SIGMA)>        
         do J=jstart, jend
            loop_count=loop_count+1
            
            firstIndexAO_j=qm2_params%orb_loc(1,J)
            lastIndexAO_j=qm2_params%orb_loc(2,J) 
            naoj = qm2_params%natomic_orbs(j)
            qmjtype = qmmm_struct%qm_atom_type(j)             
            
            hasDOrbital=((naoi.ge.9) .or. (naoj.ge.9)) 

            r2=sum( (coord(1:3,i) - coord(1:3,j)) **2) 
            rij=sqrt(r2)
            rijInAu=rij*A_TO_BOHRS
            oneOverRij=one/rij
            r2InAu=r2*A2_TO_BOHRS2
                          
            if (r2InAu < OVERLAP_CUTOFF) then 

            !Calculate Overlap Integrals using a Gaussian Expansion
            !STO-6G BY R.F. STEWART, J. CHEM. PHYS., 52 431-438, 1970
            !Fill SHMAT with a 4x4 array of overlaps, in order S,PX,PY,PZ
            !     r2InAu   =  INTERATOMIC DISTANCE^2 IN BOHRS2

                
                if (.not. hasDOrbital) then
                  ! the implementation for sp orbitals by Ross Walker 
                  
                  SHMAT=0.0d0
                  call qm2_h1elec(r2InAu,COORD(1:3,I),COORD(1:3,J),  &
                            naoi,naoj,SHMAT, qmitype, qmjtype)
                  I2=0
                  do I1=firstIndexAO_i,lastIndexAO_i
                     II=qm2_params%pascal_tri1(i1)+firstIndexAO_j-1
                     I2=I2+1
                     J2=0
                     JJ=MIN(I1,lastIndexAO_j)
                     do J1=firstIndexAO_j,JJ                                                   
                        II=II+1                                                       
                        J2=J2+1                                                       
                        H(II)=H(II)+SHMAT(I2,J2) 
                     end do
                  end do
                  
                else
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! one-electron part
                  ! the implementation for d-orbital by Taisung Lee 
                  ! basically by coping things from the MNDO program  
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  call qm2_h1elec_d(r2InAu,COORD(1:3,I),COORD(1:3,J),  &
                            naoi,naoj,firstIndexAO_i, firstIndexAO_j,  &
                            qmitype, qmjtype,qm2_struct%norbs,H)                
                
                endif
              
            end if !(R2 < OVERLAP_CUTOFF)
            
            
            ! -----------------------------------------------------
            ! Calculate two-electron integrals W and
            ! electron-nuclear terms 
            ! core-core repulsion energy for atom pair i, j
            ! -----------------------------------------------------             
            RI=0.0D0
            core=0.0D0
            call GetRotationMatrix(coord(1:3,j)-coord(1:3,i), rotationMatrix, hasDOrbital)   
            
            call qm2_rotate_qmqm(loop_count,i,j,NI,qmmm_struct%iqm_atomic_numbers(J),COORD(1,I),COORD(1,J), &
                       W(KR),KI,RI, core)

            if (hasDOrbital) then   ! spd case   
                  
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! 2-center 2-electron part and core-core repulsion
              ! 
              ! the implementaion for d-orbital by Taisung Lee 
              ! basically by coping things from the MNDO program  
              ! Note that the sp part of RI and CORE is alredy done at this point
              ! but we need to re-do them again to get CORE filled.
              ! This should be fixed later by filling CORE in the SP part.
              ! Taisung Lee, Rutgers 2011
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
                
                i_dimension=naoi*(naoi+1)/2
                j_dimension=naoj*(naoj+1)/2
                ki=i_dimension*j_dimension                    
               
                allocate(ww(1:j_dimension, 1:i_dimension)) 
                WW=0.0D0

                call qm2_repp_d(qmitype,qmjtype,rijInAu,RI,CORE,WW,i_dimension,j_dimension,1)
 
                ! WW now holds 2-center 2-electron integrals
                ! the following code copys WW to the global W storage
                ! note that both have the second atom index as the first index
 
                k=0
                do ii=1,i_dimension
                  do jj=1, j_dimension
                    W(kr+k)=WW(jj,ii)
                    k=k+1
                  end do
                end do

                deallocate(ww)
                call Rotate2Center2Electron(W(kr), i_dimension,j_dimension, rotationMatrix)
                
            end if  ! ((naoi.ge.9) .or. (naoj.ge.9))

            ii=qm2_params%pascal_tri2(firstIndexAO_i)
            jj=qm2_params%pascal_tri2(firstIndexAO_j)  
                          
            call RotateCore(firstIndexAO_i,firstIndexAO_j,              &
                naoi,naoj,  &
                ii,jj,core,rotationMatrix,H)

            call qm2_core_core_repulsion(i, j, rij, oneOverRij, RI(1), enuc)    
            enuclr = enuclr + enuc                    
    
            ! shift the pointer to the global W storage by the size of WW
            kr=kr+ki
            
         end do  ! J=1,iminus
      end do !  I=1,qmmm_struct%nquant_nlink

!DEBUG      call print ('One-electron matrix',H,.true.)
      
end subroutine qm2_hcore_qmqm

                                                                       
