#include "../include/assert.fh"
#include "../include/dprec.fh"
#include "xmin.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver routine for XMIN minimization.
!-----------------------------------------------------------------------

subroutine xmin( natom, x, fg, escf, xmin_iter, maxiter, born_radii, &
     one_born_radii, intdiel, extdiel, Arad, scf_mchg, grms_tol, ntpr )

   use constants, only : zero
   implicit none

   integer, intent(in) :: natom, maxiter, ntpr
   _REAL_ born_radii(*), one_born_radii(*), intdiel, extdiel, Arad, scf_mchg(*)
   _REAL_,  intent(inout) :: x(*)
   _REAL_,  intent(inout) :: fg(*)
   _REAL_,  intent(inout) :: escf
   integer, intent(inout) :: xmin_iter  
   _REAL_,  intent(in)  :: grms_tol
   ! ------ External functions -----------------
   _REAL_   xminC
   external xminC

   ! ------ local variables --------------------
   _REAL_  :: grms            = ZERO
   logical :: is_error
   logical :: is_xmin_done
   _REAL_  :: minimum_energy
   integer :: n_force_calls
   integer :: return_flag
   integer :: status_flag
   _REAL_  :: xmin_time
   integer :: mytaskid =0

   ! The depth of the LBFGS memory for XMIN's LBFGS minimization or TNCG
   ! preconditioning.
   ! The value 0 turns off preconditioning in TNCG minimization.
   integer  :: lbfgs_memory_depth = 3

   ! XMIN's finite difference Hv matrix-vector product method.
   ! Renaming of Kolossvary's numdiff.
   integer  :: mvpm_code = 1

   ! XMIN minimization method.
   integer,          parameter :: XMIN_METHOD_PRCG_CODE  = 1
   integer,          parameter :: XMIN_METHOD_LBFGS_CODE = 2
   integer,          parameter :: XMIN_METHOD_TNCG_CODE  = 3
   integer       :: xmin_method_code = 3

   ! Verbosity of the internal status output from the XMIN package:
   ! 0 = none, 1 = minimization details, 2 = minimization and
   ! line search details plus CG details in TNCG.
   integer,      parameter :: MAXIMUM_XMIN_VERBOSITY = 2
   integer,      parameter :: MINIMUM_XMIN_VERBOSITY = 0
   integer                 :: xmin_verbosity = 0

   integer ::  ls_method, ls_maxiter, ls_iter, xyz_min
   _REAL_  :: ls_maxatmov, beta_armijo, c_armijo, mu_armijo, ftol_wolfe, &
              gtol_wolfe

   logical :: first

   ! ------ External Functions -----------------
   _REAL_ ddot

   first = .true.
   status_flag = 0
   xyz_min = 1
   ls_method = 2
   ls_maxiter = 20
   ls_maxatmov = 0.2
   beta_armijo = 0.5
   c_armijo = 0.4
   mu_armijo = 1.0
   ftol_wolfe = 0.0001
   gtol_wolfe = 0.9

   ! keep xminC() from thinking that atoms are frozen:
   fg(1:3*natom) = 1.d0

   n_force_calls = 0
   is_error = .false.
   is_xmin_done = .false.
   do while ( .not. is_xmin_done .and. .not. is_error )

      minimum_energy = xminC( xyz_min, xmin_method_code, maxiter, grms_tol, &
            natom, lbfgs_memory_depth, mvpm_code, &
            x, escf, fg, grms, xmin_iter, xmin_time, &
            xmin_verbosity, ls_method, ls_maxiter, ls_iter, ls_maxatmov, &
            beta_armijo, c_armijo, mu_armijo, ftol_wolfe, gtol_wolfe,  &
            return_flag, status_flag )

      select case ( return_flag )
      case ( DONE )
         ! Finished minimization.
         is_xmin_done = .true.
         is_error = status_flag < 0
         write(6,'(a)') '  ... geometry converged !'
      case ( CALCENRG, CALCGRAD, CALCBOTH )
         ! Normal Amber control of NB list updates.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            call sqm_energy( natom, x, escf, born_radii, one_born_radii, &
                 intdiel, extdiel, Arad, scf_mchg ) 
            call sqm_forces(natom, fg)
            n_force_calls = n_force_calls + 1
         end if
      case ( CALCENRG_NEWNBL, CALCGRAD_NEWNBL, CALCBOTH_NEWNBL )
         ! Coerce a NB list update.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            if( n_force_calls>0 .and. mod(xmin_iter,ntpr)==0 ) then
               grms = sqrt(ddot(3*natom,fg,1,fg,1)/dble(3*natom))
               if (first) then
                  write(6,'(a)') '      iter         sqm energy              rms gradient'
                  write(6,'(a)') '      ----    -------------------    -----------------------'
                  first = .false.
               end if
               write(6,'(a,i5,f14.4,a,f14.4,a)') 'xmin ', xmin_iter, escf,' kcal/mol', grms, ' kcal/(mol*A)'
            end if
            call sqm_energy( natom, x, escf, born_radii, one_born_radii, &
                 intdiel, extdiel, Arad, scf_mchg ) 
            call sqm_forces(natom, fg)
            n_force_calls = n_force_calls + 1
         end if
      case ( CALCENRG_OLDNBL, CALCGRAD_OLDNBL, CALCBOTH_OLDNBL )
         ! Prevent a NB list update.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            call sqm_energy( natom, x, escf, born_radii, one_born_radii, &
                 intdiel, extdiel, Arad, scf_mchg ) 
            call sqm_forces(natom, fg)
            n_force_calls = n_force_calls + 1
         end if
      case default
         ! error from XMIN or the return_flag is corrupted.
         is_error = status_flag < 0
         ASSERT( is_error )
      end select

   end do

   if ( is_error ) then
      write(6,'(a,i4)') '  XMIN ERROR: Status is ', status_flag
      call mexit(6,1)
   end if

   return

end subroutine xmin
