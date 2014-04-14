subroutine qm2_transf_roothan_diag(na,nb,n,a,b,eigen_val,error)

   use qmmm_module, only: qm2_struct
   
   implicit none

!Passed in:
   integer, intent(in)    :: na, nb       ! Dimensions of the matrices a, b
   integer, intent(in)    :: n
   _REAL_,  intent(inout) :: a(na,n)      ! IN: H matrix, OUT: eigenvectors
   _REAL_,  intent(in)    :: b(nb,n)      ! The overlap matrix
   _REAL_,  intent(out)   :: eigen_val(n) ! Eigenvalues
   integer, intent(out)   :: error        ! !=0 if error.

!locals:
   integer :: i,j
   _REAL_  :: X(nb,n)                     ! The transformation matrix
   _REAL_  :: XT(n,nb)                    ! X transpose
   _REAL_  :: U(nb,n)                     ! Eigenvectors of S
   _REAL_  :: s(nb,n)                     ! Eigenvalues  of S
   _REAL_  :: ta(nb,n)                    ! Transformed 'a'
   _REAL_  :: tc(nb,n)                    ! Transformed 'c' (eigenvectors)

!
! This routine uses Ross's diagonalizer (qm2_mat_diag, in qm2_scf.f) to solve
! a non-orthogonal eigenvalue problem. For this, we use the tranformed 
! roothan equations. See Szabo & Ostlund, sec. 3.4.5 for details. In short,
! given a problem like HC = SCe (non orthogonal), we do:
!
!    1. Diagonalize S : Eigenvectors --> U; 
!                       Eigenvalues  --> s
!    2. Build the transformation matrix X: X  = U s**(-1/2)
!    3. Build the transformed H matrix:    H' = XT H X
!    4. Diagonalize H': Eigenvectors --> C'
!                       Eigenvalues  --> e
!    5. 'UN'transform the eigenvectors:    C  = X C'
!
! Gustavo Seabra, 2005.
!

