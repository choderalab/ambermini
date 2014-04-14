! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++
!This subroutine reads the QMMM namelist
!and also calls allocation routines
!for QMMM based on natom.
!
!Author:
!     Ross Walker (SDSC)
!+++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Reads the qmmm namelist and calls the qmmm memory allocation routines
#ifdef SQM
subroutine read_qmmm_nm_and_alloc( natom_inout, igb, atnam, atnum, maxcyc, &
            grms_tol, ntpr, ncharge_in, excharge, chgatnum )
#else
subroutine read_qmmm_nm_and_alloc( igb, ih, ix, x, cut, use_pme, ntb, qmstep, isabfqm, abfqmcharge) ! lam81
#endif

   use findmask
   use constants, only : RETIRED_INPUT_OPTION
   use qmmm_module, only : qmmm_struct, qm2_struct, qmmm_nml, &
                           validate_qm_atoms, qmsort, &
                           allocate_qmmm, get_atomic_number, qmmm_div, &
                           qmmm_opnq, qmmm_vsolv
   use qmmm_vsolv_module, only : read_vsolv_nml
   use qmmm_qmtheorymodule
   use ElementOrbitalIndex, only : numberElements
   use ParameterReader, only : ReadParameterFile
    
   implicit none

!STATIC MEMORY
integer :: max_quantum_atoms  !Needed in read qmmm namelist since namelists cannot contain pointers
parameter ( max_quantum_atoms = 10000 )
!END STATIC MEMORY

!Passed in
   integer :: igb        !Value of igb from cntrl namelist
#ifdef SQM
   integer use_pme, ntb
   integer, intent(inout) :: natom_inout
   character(len=8), intent(in) :: atnam(*)
   integer, intent(in)  :: atnum(*)
   integer, intent(out) :: maxcyc, ntpr
   _REAL_, intent(out)  :: grms_tol
   _REAL_, intent(in) :: excharge(*)
   integer, intent(in) :: chgatnum(*)
   integer, intent(in) :: ncharge_in
#else
   character(len=4) ih(*)
   integer, intent(in) :: ix(*)
   _REAL_, intent(in) :: x(*)
   _REAL_, intent(in) :: cut !MM-MM cutoff in angstroms
   integer, intent(in) :: use_pme, ntb
   integer, intent(in) :: qmstep      ! lam81
   integer, intent(in), optional :: isabfqm(*)   ! lam81
   integer, intent(in), optional :: abfqmcharge  ! lam81
#endif

!local  
   _REAL_ :: qmcut      ! local copied to qmmm_nml%qmcut - specified cutoff to use for QM-MM electrostatics.
                         ! Default = same as regular MM cutoff.
   _REAL_ :: lnk_dis     ! Distance from the QM atom to place link atom.
                         !A value of <0.0 means the link atom gets placed on the MM link pair atom's coordinates
                         !on every MD step.
   _REAL_ :: scfconv     ! local copied to qmmm_nml%scfconv - Convergence criteria for SCF routine. Default = 1.0D-8.
                         ! Minimum (tightest criteria) = 1.0D-16

   !+TJG 01/26/2010
   _REAL_ :: errconv      ! Convergence criteria for maximum value in the error matrix FP-PF
   integer :: ndiis_matrices ! Maximum number of matrices used in a diis extrapolation
   integer :: ndiis_attempts ! The initial number of diis tokens... the number of iterations
                             ! that diis extrapolations will be attempted before giving up on diis
   !-TJG 01/26/2010

   _REAL_ :: pseudo_diag_criteria !Criteria - maximum change in density matrix between successive SCF iterations
                                  !in which to allow pseudo diagonalisation. Default = 0.05.
   integer :: lnk_atomic_no !Atomic number of link atom
   integer :: lnk_method !controls how QM-MM valence terms are dealt with.
   integer :: qmgb       ! local copied to qmmm_nml%qmgb - flag for type of GB do with QM region
   integer :: qmtheory   ! deprecated flag for level of theory to use for QM region
   integer :: qmcharge   ! local copied to qmmm_nml%qmcharge - value of charge on QM system
   integer :: corecharge    ! lam81
   integer :: buffercharge  ! lam81
   integer :: spin       ! local copied to qmmm_nml%spin - spin state of system
   integer :: i,j        ! temporary counter for local loops
   integer :: ifind
   integer :: qmqmdx     ! local copied to qmmm_nml%qmqm_analyt - 1 = analytical, 2 = numerical QM-QM derivatives in qm2
   integer :: verbosity  ! local copied to qmmm_nml%verbosity - Controls amount of info about QM part of calc that is printed (0=Def)
   integer :: tight_p_conv ! local copied to qmmm_nml%tight_p_conv - Controls convergence of density matrix. 0 (Def) = 0.05*sqrt(SCFCRT)
                           ! 1 = Converged both Energy and Density to SCFCONV
   integer :: printcharges !Local copied to qmmm_nml%printcharges as a logical. 1 = true - print mulliken and cm1a and cm2a charges
                           !on every step. 0 = false = don't print charges. Default = 0 (.false.)
   integer :: printdipole  !Local copied to qmmm_nml%printdipole as an integer 1 = QM dipole moment, 2 = QM + MM dipole moment, (0=Def) 
   integer :: print_eigenvalues  !Local copied to qmmm_nml%print_eigenvalues, 0 = no printing, 1 = at end of run (default), 2 = each SCF cycle, 3 = each SCF iteration
   integer :: peptide_corr !Local copied to the logical qmmm_nml%peptide_corr
                           !Add MM correction to peptide linkages 0 = No (Default), 1 = Yes.
   integer :: itrmax       !Local copied to qmmm_nml%itrmax - Maximum number of scf cycles to run
                           !before assuming convergence has failed (default = 1000)
   integer :: printbondorders !Local copied to qmmm_nml%printbondorders as a logical. 
                              ! 1 = true - print bondorders at the end of the
                              ! calculation. 0 = false = dont print bondorders. Default = 0 (.false.)
   integer :: qmshake      !Local copied to qmmm_nml%qmshake - shake QM atoms if ntc>1?
   integer :: qmmmrij_incore !Flag to store rij between qm-mm pairs and related equations in memory.
                             !1 (default) = store in memory. 0 = calc on fly.
   integer :: qmqm_erep_incore !Flag to store QM-QM 1 electron repulsion integrals in memory or to calculate
                               !them on the fly. Only available with QM-QM analytical derivatives.
                               !1 (default) = store in memory. 0 = calc on fly.
   integer :: pseudo_diag      !Whether to allow pseudo diagonalisations to be done when possible in SCF.
                               !0 (default) = Always do full diagonalisations.
                               !1 = do pseudo diagonalisations when possible.
   integer :: qm_ewald          !0 (default) do only regular QM-MM interaction in periodic calculations.
                                !1           do ewald based periodic QM-MM interactions.
                                !2           do ewald based periodic QM-MM but with the QM image charges
                                !            fixed at the previous steps mulliken charges during the SCF.
   integer :: qm_pme            !0 use regular Ewald for doing QM-MM interactions when qm_ewald>0.
                                !1 (default) use PME to do the reciprocal sum.
   integer :: kmaxqx, kmaxqy, kmaxqz !Maximum K space vectors
   integer :: ksqmaxq !Maximum K squared values for spherical cutoff in k space.
   _REAL_ :: kappa  ! the ewald coefficient for QM region ewald calculations
   integer :: writepdb
   integer :: qmmm_int !QM-MM interaction method
   integer :: adjust_q
   integer :: diag_routine !Controls diagonalization routine to use in SCF.
#ifdef OPENMP
   integer :: qmmm_omp_max_threads !Maximum number of openmp threads to use for parallel QMMM routines
#endif
   integer :: density_predict !Controls prediction of density matrix for next SCF step.
   integer :: fock_predict !Controls prediction of Fock matrix for next SCF step.
   integer :: vsolv ! = 0 by default, = 1 for simple vsolv QM/MM, = 2 for adaptive QM/MM with vsolv

   _REAL_ :: fockp_d1 !prefactor for fock matrix prediction.
   _REAL_ :: fockp_d2 !prefactor for fock matrix prediction.
   _REAL_ :: fockp_d3 !prefactor for fock matrix prediction.
   _REAL_ :: fockp_d4 !prefactor for fock matrix prediction.
   _REAL_ :: damp   ! AWG SCF damping
   _REAL_ :: vshift ! AWG SCF level shift parameter
   logical :: mdin_qmmm=.false.

   integer :: idc
   integer :: divpb
!! (GMS)
   _REAL_  :: chg_lambda       ! Charge scaling factor for free energy calculation
!! DFTB options
   integer :: dftb_maxiter     ! Max # of iterations before resetting Broyden (default: 70 ) ==> qmmm_nml%dftb_maxiter
   integer :: dftb_disper      ! Use dispersion?  (default: 0 = false) ==> qmmm_nml%dftb_disper
   integer :: dftb_chg         ! DFTB CM3 charges (default: 0 = Mulliken, 1 = CM3) ==> qmmm_nml%dftb_chg
   _REAL_  :: dftb_telec       ! Electronic temperature, in Kelvins. (Default = 0.0K) ==> qmmm_nml%dftb_telec 
   _REAL_  :: dftb_telec_step  ! Telec step size for convergence accelerator (Default = 0.0K) ==> qmmm_nml%dftb_telec_step
   character(Len=256) :: dftb_3rd_order  ! 3rd order SCC-DFTB (default: 'NONE'== No third order)
                                       !     'PA' == Do 3rd order, Proton Affinities parameterization
                                       !     'PR' ==               Phosphate reactions parameterization
                                       !     'READ' == read the parameters from a user-specified file (TO IMPLEMENT)
   _REAL_ :: r_switch_lo    !Lower bound of the QM/MM switching function
   _REAL_ :: r_switch_hi    !Upper bound of the QM/MM switching function
   integer :: qmmm_switch   !0           Turn off QM/MM switching function 
                            !1           Turn on QM/MM switching function 

   integer :: abfqmmm       ! lam81
   integer :: hot_spot      ! lam8
   _REAL_ :: r_core_in      ! lam81
   _REAL_ :: r_core_out     ! lam81
   _REAL_ :: r_qm_in        ! lam81
   _REAL_ :: r_qm_out       ! lam81
   _REAL_ :: r_buffer_in    ! lam81
   _REAL_ :: r_buffer_out   ! lam81

   character(len=256) :: cut_bond_list_file          ! lam81
   character(len=256) :: oxidation_number_list_file  ! lam81

   integer :: mom_cons_type            ! lam81
   integer :: mom_cons_region          ! lam81

   integer :: fix_atom_list            ! lam81
   integer :: solvent_atom_number      ! lam81

   integer :: selection_type           ! lam81
   integer :: center_type              ! lam81
   integer :: initial_selection_type   ! lam81

   integer :: max_bonds_per_atom       ! lam81
   integer :: n_max_recursive          ! lam81

   _REAL_ :: min_heavy_mass            ! lam81

   _REAL_ :: gamma_ln_qm               ! lam81

   character(len=256) :: read_idrst_file  ! lam81
   character(len=256) :: write_idrst_file ! lam81
   integer :: ntwidrst                    ! lam81

   character(len=256) :: pdb_file      ! lam81
   integer :: ntwpdb                   ! lam81

#include "../include/memory.h"

   !Apparently you can't use a pointer in a namelist :-( Therefore
   !we need a local scratch array that will be big enough that 
   !the iqmatoms list never exceeds it
   integer :: iqmatoms( max_quantum_atoms )
   integer :: core_iqmatoms( max_quantum_atoms )    ! lam81
   integer :: buffer_iqmatoms( max_quantum_atoms )  ! lam81

   character(len=8192) :: qmmask        ! lam81
   character(len=8192) :: coremask      ! lam81
   character(len=8192) :: buffermask    ! lam81

   integer :: qm_subsetatoms( natom )      ! lam81
   integer :: core_subsetatoms( natom )    ! lam81
   integer :: buffer_subsetatoms( natom )  ! lam81

   integer :: center_subsetatoms( natom )  ! lam81

   character(len=8192) :: ext_qmmask_subset      ! lam81
   character(len=8192) :: ext_coremask_subset    ! lam81
   character(len=8192) :: ext_buffermask_subset  ! lam81

   character(len=8192) :: centermask             ! lam81

   character(len=12) :: qm_theory 
        !Options=PM3,AM1,MNDO,PDDG-PM3,PM3PDDG,PDDG-MNDO,PDDGMNDO,
        !        PM3-CARB1,PM3CARB1,DFTB,SCC-DFTB,RM1,PM6,PM3-ZnB,PM3-MAIS
        !        EXTERNAL (for external programs like ADF/GAMESS/TeraChem)
   integer, dimension(:), pointer :: isqm
   integer :: ier=0
   character(len=80) :: parameter_file   
   logical :: qxd
   logical :: test

   namelist /qmmm/ qmcut, iqmatoms,qmmask,qmgb,qm_theory, qmtheory, &
                   qmcharge, qmqmdx, verbosity, tight_p_conv, scfconv, &
                   errconv,ndiis_matrices,ndiis_attempts, &  !+TJG 01/26/2010
                   parameter_file, qxd,&
                   printcharges, printdipole, print_eigenvalues, peptide_corr, itrmax, qmshake, &
                   qmqm_erep_incore, qmmmrij_incore, &
                   lnk_dis, lnk_atomic_no, lnk_method, spin, pseudo_diag,   &
                   pseudo_diag_criteria, &
                   qm_ewald, qm_pme, kmaxqx, kmaxqy, kmaxqz, ksqmaxq, kappa, &
                   writepdb, qmmm_int, adjust_q, diag_routine, &
                   density_predict, fock_predict, &
                   fockp_d1, fockp_d2, fockp_d3, fockp_d4, idc, divpb, &
                   dftb_maxiter, dftb_disper, dftb_3rd_order, dftb_chg, &
                   dftb_telec, dftb_telec_step, printbondorders, &
                   qmmm_switch, r_switch_lo, r_switch_hi, damp, vshift, &
#ifdef OPENMP
                   qmmm_omp_max_threads, &
#endif
#ifdef SQM
                   maxcyc, ntpr, grms_tol, &
#endif
                   chg_lambda, vsolv, &                                              ! lam81
                   abfqmmm, r_core_in, r_core_out, r_buffer_in, r_buffer_out, &      ! lam81
                   r_qm_in, r_qm_out, coremask, buffermask, &                        ! lam81
                   cut_bond_list_file, oxidation_number_list_file, &                 ! lam81
                   mom_cons_type, mom_cons_region, &                                 ! lam81
                   fix_atom_list, solvent_atom_number, &                             ! lam81
                   selection_type, center_type, initial_selection_type, &            ! lam81
                   corecharge, buffercharge, &                                       ! lam81
                   max_bonds_per_atom, n_max_recursive, min_heavy_mass, &            ! lam81
                   gamma_ln_qm, &                                                    ! lam81
                   read_idrst_file, ntwidrst, write_idrst_file, ntwpdb, pdb_file, &  ! lam81
                   ext_qmmask_subset, ext_coremask_subset, ext_buffermask_subset, &  ! lam81
                   centermask, hot_spot                                              ! lam81

!Setup defaults
#ifdef SQM
   qmcut = 9999.d0
   use_pme = 0
   ntb = 0
   maxcyc = 9999
   grms_tol = 0.02
   ntpr=10
#else
   qmcut = cut
#endif
   lnk_dis=1.09d0  !Methyl C-H distance
   lnk_atomic_no=1 !Hydrogen
   lnk_method=1 !treat MMLink as being MM atom.
   qmgb = 2 !Gets set to zero if igb==6 or igb==0.
   qm_theory = ''
   qmtheory = RETIRED_INPUT_OPTION 
   qmcharge = 0
   corecharge = 0    ! lam81
   buffercharge = 0  ! lam81
   spin = 1
   qmqmdx = 1
   verbosity = 0
   parameter_file=''
   qxd=.false.
   !+TJG 01/26/2010
#ifdef SQM
   ! defaults for stand-alone (geom. opt.) code:
   !   dac: for now, just use the same values as in the non-stand-alone
   !        code, but these should be updated once we figure out the best
   !        values for geometry optimization
   tight_p_conv = 0
   scfconv = 1.0D-8
   errconv = 1.0d-1
   ndiis_matrices = 6
   ndiis_attempts = 0
#else
   ! defaults for MD is off:
   tight_p_conv = 0
   scfconv = 1.0D-8
   errconv = 1.0D-1
   ndiis_matrices = 6
   ndiis_attempts = 0
#endif
   !-TJG 01/26/2010
   printcharges = 0
   printbondorders = 0
   printdipole = 0
   print_eigenvalues = 1
   peptide_corr = 0
   itrmax = 1000
   qmshake = 1
   qmmask=''
   coremask=''   ! lam81
   buffermask='' ! lam81
   iqmatoms(1:max_quantum_atoms) = 0
   core_iqmatoms(1:max_quantum_atoms) = 0    ! lam81
   buffer_iqmatoms(1:max_quantum_atoms) = 0  ! lam81
   ext_qmmask_subset=''      ! lam81
   ext_coremask_subset=''    ! lam81
   ext_buffermask_subset=''  ! lam81
   centermask=''             ! lam81
   qm_subsetatoms(1:natom) = 0       ! lam81
   core_subsetatoms(1:natom) = 0     ! lam81
   buffer_subsetatoms(1:natom) = 0   ! lam81
   center_subsetatoms(1:natom) = 0   ! lam81
   qmmmrij_incore = 1
   qmqm_erep_incore = 1
   pseudo_diag = 1
   pseudo_diag_criteria = 0.05d0
   qm_ewald=1 !Default is to do QMEwald, with varying charges, if ntb=0 or use_pme=0 then this will get turned off
   qm_pme = 1 !use pme for QM-MM
   kmaxqx=8; kmaxqy=8; kmaxqz=8    !Maximum K space vectors
   kappa=-1.0
   ksqmaxq=100 !Maximum K squared values for spherical cutoff in k space.
   writepdb = 0 !Set to 1 to write a pdb on the first step with just the QM region in it.
   qmmm_int = 1 !Default, do full interaction without extra Gaussian terms for PM3 / AM1 etc.
   adjust_q = 2 !Default adjust q over all atoms.
   diag_routine = 1 !Use default internal diagonalizer.
#ifdef OPENMP
   qmmm_omp_max_threads = 1 !Use just 1 openmp thread by default.
#endif
   density_predict = 0 !Use density matrix from previous MD step.
   fock_predict = 0 !Do not attempt to predict the Fock matrix.
   fockp_d1 = 2.4d0
   fockp_d2 = -1.2d0
   fockp_d3 = -0.8d0
   fockp_d4 = 0.6d0
   damp = 1.0
   vshift = 0.0
   idc = 0
   divpb = 0
   vsolv = 0 ! by default do not use simple vsolv QM/MM or adaptive QM/MM based on vsolv
   qmmm_switch = 0              !Use QM/MM switching function
   r_switch_hi = qmcut          !Set the default value to be equal to qmcut
   r_switch_lo = r_switch_hi - 2.0D0  !Set the default value to be 2 Angstrom shorter than r_switch_hi

   !DFTB
   dftb_maxiter     = 70   
   dftb_disper      = 0
   dftb_chg         = 0
   dftb_telec       = 0.0d0
   dftb_telec_step  = 0.0d0
   chg_lambda  = 1.0d0
   dftb_3rd_order   = 'NONE'

   !ABFQMMM          ! lam81
   abfqmmm      = 0  ! lam81
   hot_spot     = 0  ! lam81
   r_core_in    = 0  ! lam81
   r_core_out   = 0  ! lam81
   r_qm_in      = 0  ! lam81
   r_qm_out     = 0  ! lam81
   r_buffer_in  = 0  ! lam81
   r_buffer_out = 0  ! lam81

   cut_bond_list_file = ''          ! lam81
   oxidation_number_list_file = ''  ! lam81

   mom_cons_type = 1            ! lam81
   mom_cons_region = 1          ! lam81

   fix_atom_list = 0            ! lam81
   solvent_atom_number = 3      ! lam81

   selection_type = 1           ! lam81
   center_type = 1              ! lam81
   initial_selection_type = 0   ! lam81

   max_bonds_per_atom = 4       ! lam81
   n_max_recursive = 10000      ! lam81

   min_heavy_mass = 4.0         ! lam81

   gamma_ln_qm = 0.0d0          ! lam81

   read_idrst_file = ''                 ! lam81
   write_idrst_file = 'abfqmmm.idrst'   ! lam81
   ntwidrst = 0                         ! lam81 

   pdb_file = 'abfqmmm.pdb'     ! lam81
   ntwpdb = 0                   ! lam81

   !Read qmmm namelist
   rewind 5

   call nmlsrc('qmmm',5,ifind)
   if (ifind /= 0) mdin_qmmm=.true.

   !Read qmmm namelist
   rewind 5
   if ( mdin_qmmm ) then
     read(5,nml=qmmm)
   else
     write(6, '(1x,a,/)') 'Could not find qmmm namelist'
     call mexit(6,1)
   endif

   !AWG NEW
   call CheckRetiredQmTheoryInputOption(qmtheory)
   call set(qmmm_nml%qmtheory, qm_theory)
   !AWG END NEW
   
   !  Read-in the user-defined parameter file 
    !  TL (Rutgers, 2011)
    call ReadParameterFile(parameter_file)
    
    ! turn on OPNQ if necessary 
    qmmm_opnq%useOPNQ=qxd

#ifdef SQM
   ! Disable EXTERN in SQM since
   ! it does not make sense for SQM to be calling the external ADF interface.
   if (qmmm_nml%qmtheory%EXTERN) then                                                  
      call sander_bomb('read_qmmm_namelist','External interface is not supported in SQM.', &       
           '(qm_theory = ''EXTERN'')')                                                 
   end if

   if (qmmm_nml%qmtheory%SEBOMD) then
      call sander_bomb('read_qmmm_namelist','SEBOMD interface is not supported in SQM.', &
           '(qm_theory = ''SEBOMD'')')
   end if

   if (ncharge_in > 0) then ! we have external charge
      if (maxcyc > 0) then
         ! external charge calculation is not supported in gradient minimization
         call sander_bomb('read_qmmm_namelist','maxcyc > 0 but external charge found.', &
                       'external charge calculation is for single point energy calculations only.')
      end if
      natom = natom_inout + ncharge_in
      qmmm_struct%nquant = natom_inout
      qmmm_struct%nlink  = 0
      qmmm_struct%nquant_nlink = natom_inout
      qmmm_struct%qm_mm_pairs = ncharge_in
      do i=1,natom_inout
         iqmatoms(i) = i
      end do
      !update natom_inout with the number of external charges
      natom_inout = natom
   else
      natom = natom_inout
      qmmm_struct%nquant = natom
      qmmm_struct%nlink  = 0
      qmmm_struct%nquant_nlink = natom
      do i=1,natom
         iqmatoms(i) = i
      end do
   end if
   !natom = natom_in
   !qmmm_struct%nquant = natom
   !qmmm_struct%nlink  = 0
   !qmmm_struct%nquant_nlink = natom
   !do i=1,natom
   !   iqmatoms(i) = i
   !end do
#else

   if (qmmm_nml%qmtheory%SEBOMD) then
     qmmm_struct%nquant = 0
   else

    if( (qmmask /= '') .or. (coremask /= '') ) then  !  get the quantum atoms from the mask
       if(abfqmmm == 1 .and. qmmm_struct%abfqmmm /= 1) then                 ! lam81
        write(6,'(a)') ''                                                   ! lam81
        write(6,'(/80("-")/"   ADAPTIVE BUFFERED FORCE QM/MM",/80("-")/)')  ! lam81
        if(hot_spot == 1) then                                              ! lam81
         write(6,'(a)') ''                                                  ! lam81
         write(6,'(/80("-")/"   HOT SPOT IS ACTIVE           ",/80("-")/)') ! lam81
        end if                                                              ! lam81
       end if

       if(qmmm_struct%abfqmmm /= 1) then                                                  ! lam81
        if(abfqmmm == 1) then                                                             ! lam81
         write(6,'(a)') ''                                                                ! lam81
         write(6,'(a)') '------------------------'                                        ! lam81
         write(6,'(a)') 'Specification of regions'                                        ! lam81
         write(6,'(a)') '------------------------'                                        ! lam81
         write(6,'(a)') ''                                                                ! lam81
         write(6,'(a)') 'QM atoms:'                                                       ! lam81
         write(6,'(a)') '---------'                                                       ! lam81
         write(6,'(a)') ''                                                                ! lam81
         if(qmmask /= '') then                                                            ! lam81
          write(6,'(a)') 'INFO: loading the quantum atoms as groups'                      ! lam81
         else                                                                             ! lam81
          write(6,'(a)') 'INFO: quantum atoms for adaptive QM/MM are not defined'         ! lam81
          write(6,'(a)') 'INFO: core atoms will be used in reduced calculation'           ! lam81
         end if                                                                           ! lam81
        else                                                                              ! lam81
          write(6,'(a)') 'LOADING THE QUANTUM ATOMS AS GROUPS'                            ! lam81
        end if                                                                            ! lam81
       end if                                                                             ! lam81

       allocate(isqm( natom ), stat=ier)
       REQUIRE(ier==0)

       call atommask( natom, nres, 0, ih(m04), ih(m06), &
          ix(i02), ih(m02), x(lcrd), qmmask, isqm )

       if(abfqmmm == 1 .and. qmstep == 0) then      ! lam81
        if(r_core_in == 0) then                     ! lam81
         qmmm_struct%r_core_in = r_core_out         ! lam81
        else                                        ! lam81
         qmmm_struct%r_core_in = r_core_in          ! lam81
        end if                                      ! lam81
        if(r_core_out <= r_core_in) then            ! lam81
         qmmm_struct%r_core_out = r_core_in         ! lam81
        else                                        ! lam81
         qmmm_struct%r_core_out = r_core_out        ! lam81
        end if                                      ! lam81
        if(r_qm_in == 0) then                       ! lam81
         qmmm_struct%r_qm_in = r_qm_out             ! lam81
        else                                        ! lam81
         qmmm_struct%r_qm_in = r_qm_in              ! lam81
        end if                                      ! lam81
        if(r_qm_out <= r_qm_in) then                ! lam81
         qmmm_struct%r_qm_out = r_qm_in             ! lam81
        else                                        ! lam81
         qmmm_struct%r_qm_out = r_qm_out            ! lam81
        end if                                      ! lam81
        if(r_buffer_in == 0) then                   ! lam81
         qmmm_struct%r_buffer_in = r_buffer_out     ! lam81
        else                                        ! lam81
         qmmm_struct%r_buffer_in = r_buffer_in      ! lam81
        end if                                      ! lam81
        if(r_buffer_out <= r_buffer_in) then        ! lam81
         qmmm_struct%r_buffer_out = r_buffer_in     ! lam81
        else                                        ! lam81
         qmmm_struct%r_buffer_out = r_buffer_out    ! lam81
        end if                                      ! lam81

        qmmm_struct%cut_bond_list_file = cut_bond_list_file                  ! lam81
        qmmm_struct%oxidation_number_list_file = oxidation_number_list_file  ! lam81

        qmmm_struct%mom_cons_type = mom_cons_type                      ! lam81
        qmmm_struct%mom_cons_region = mom_cons_region                  ! lam81

        qmmm_struct%fix_atom_list = fix_atom_list                      ! lam81
        qmmm_struct%solvent_atom_number = solvent_atom_number          ! lam81

        qmmm_struct%selection_type = selection_type                    ! lam81
        qmmm_struct%center_type = center_type                          ! lam81
        qmmm_struct%initial_selection_type = initial_selection_type    ! lam81

        qmmm_struct%max_bonds_per_atom = max_bonds_per_atom            ! lam81
        qmmm_struct%n_max_recursive = n_max_recursive                  ! lam81

        qmmm_struct%min_heavy_mass = min_heavy_mass                    ! lam81

        qmmm_struct%gamma_ln_qm = gamma_ln_qm                          ! lam81

        qmmm_struct%read_idrst_file = read_idrst_file                  ! lam81
        qmmm_struct%write_idrst_file = write_idrst_file                ! lam81
        qmmm_struct%ntwidrst = ntwidrst                                ! lam81

        qmmm_struct%pdb_file = pdb_file                                ! lam81
        qmmm_struct%ntwpdb = ntwpdb                                    ! lam81

       end if                                       ! lam81

       if (qmstep /= 0) isqm(1:natom) = isabfqm(1:natom)  ! lam81

       qmmm_struct%nquant = sum(isqm(1:natom))
       if( (qmmm_struct%abfqmmm /= 1) .and. (qmmask /= '') ) then                                        ! lam81
        write(6,'(a,a,a,i5,a)') '     Mask ', qmmask(1:len_trim(qmmask)), &                              ! lam81
          ' matches ',qmmm_struct%nquant,' atoms'                                                        ! lam81
       end if                                                                                            ! lam81
       if(abfqmmm == 1 .and. qmmm_struct%abfqmmm /= 1) then                                              ! lam81
        write(6,'(a,i2)') '     qm-charge: ', qmcharge                                                   ! lam81
        qmmm_nml%qmcharge = qmcharge                                                                     ! lam81
       end if

       j = 0
       do i=1,natom
          if( isqm(i)>0 ) then
             j = j+1
             iqmatoms(j) = i
          end if
       end do

       if(abfqmmm == 1 .and. qmmm_struct%abfqmmm /= 1) then                   ! lam81

        if( ext_qmmask_subset /= '' ) then                                                             ! lam81
         write(6,'(a)') 'INFO: qm subset was specified for the extended qm region'                     ! lam81
         write(6,'(a)') 'INFO: loading qm subset atoms'                                                ! lam81
         
         call atommask( natom, nres, 0, ih(m04), ih(m06), &                                            ! lam81
          ix(i02), ih(m02), x(lcrd), ext_qmmask_subset, isqm )                                         ! lam81

         qmmm_struct%qm_nsubset = sum(isqm(1:natom))                                                   ! lam81
         write(6,'(a,a,a,i5,a)') '     Mask ', ext_qmmask_subset(1:len_trim(ext_qmmask_subset)), &     ! lam81
         ' matches ',qmmm_struct%qm_nsubset,' atoms'                                                   ! lam81
         if(qmmm_struct%qm_nsubset == 0) then                                                          ! lam81
          write(6,'(a)') 'INFO: qm subset is an empty set => all atoms can be in extended qm region'   ! lam81
         end if                                                                                        ! lam81

         j = 0                           ! lam81
         do i=1,natom                    ! lam81
            if( isqm(i)>0 ) then         ! lam81
               j = j+1                   ! lam81
               qm_subsetatoms(j) = i     ! lam81
            end if                       ! lam81
         end do                          ! lam81

        else                                                                                ! lam81
         write(6,'(a)') 'INFO: qm subset was not specified for the extended qm region'      ! lam81

         qmmm_struct%qm_nsubset = 0                                                         ! lam81
        end if                                                                              ! lam81

        write(6,'(a)') ''                                                    ! lam81
        write(6,'(a)') 'CORE atoms:'                                         ! lam81
        write(6,'(a)') '-----------'                                         ! lam81
        write(6,'(a)') ''                                                    ! lam81
        if( coremask /= '' ) then                                             ! lam81
          write(6,'(a)') 'INFO: loading the core atoms for adaptive QM/MM'    ! lam81

          call atommask( natom, nres, 0, ih(m04), ih(m06), &               ! lam81
           ix(i02), ih(m02), x(lcrd), coremask, isqm )                     ! lam81

          qmmm_struct%core_nquant = sum(isqm(1:natom))                     ! lam81
          write(6,'(a,a,a,i5,a)') '     Mask ', coremask(1:len_trim(coremask)), &     ! lam81
          ' matches ',qmmm_struct%core_nquant,' atoms'                                ! lam81

          j = 0                           ! lam81
          do i=1,natom                    ! lam81
             if( isqm(i)>0 ) then         ! lam81
                j = j+1                   ! lam81
                core_iqmatoms(j) = i      ! lam81
             end if                       ! lam81
          end do                          ! lam81

        else                                                                            ! lam81
         write(6,'(a)') 'INFO: core atoms for adaptive QM/MM are not defined'           ! lam81
         write(6,'(a)') 'INFO: reduced calculation will be full MM calculation'         ! lam81
         write(6,'(a)') 'INFO: using FF parameters from topology file'                  ! lam81
         
         qmmm_struct%core_nquant = 0 ! lam81

        end if                                        ! lam81
        write(6,'(a,i2)') '     core-charge: ', corecharge              ! lam81
        qmmm_nml%corecharge = corecharge                                ! lam81

        if( ext_coremask_subset /= '' ) then                                                             ! lam81
         write(6,'(a)') 'INFO: core subset was specified for the extended core region'                   ! lam81
         write(6,'(a)') 'INFO: loading core subset atoms'                                                ! lam81

         call atommask( natom, nres, 0, ih(m04), ih(m06), &                                              ! lam81
          ix(i02), ih(m02), x(lcrd), ext_coremask_subset, isqm )                                         ! lam81

         qmmm_struct%core_nsubset = sum(isqm(1:natom))                                                   ! lam81
         write(6,'(a,a,a,i5,a)') '     Mask ', ext_coremask_subset(1:len_trim(ext_coremask_subset)), &   ! lam81
         ' matches ',qmmm_struct%core_nsubset,' atoms'                                                   ! lam81
         if(qmmm_struct%core_nsubset == 0) then                                                          ! lam81
          write(6,'(a)') 'INFO: core subset is an empty set => all atoms can be in extended core region' ! lam81
         end if                                                                                          ! lam81

         j = 0                               ! lam81
         do i=1,natom                        ! lam81
            if( isqm(i)>0 ) then             ! lam81
               j = j+1                       ! lam81
               core_subsetatoms(j) = i       ! lam81
            end if                           ! lam81
         end do                              ! lam81

        else                                                                                  ! lam81
         write(6,'(a)') 'INFO: core subset was not specified for the extended core region'    ! lam81

         qmmm_struct%core_nsubset = 0                                                         ! lam81
        end if                                                                                ! lam81

        write(6,'(a)') ''                                               ! lam81
        write(6,'(a)') 'BUFFER atoms:'                                  ! lam81
        write(6,'(a)') '-------------'                                  ! lam81
        write(6,'(a)') ''                                               ! lam81
        if( buffermask /= '' ) then                                     ! lam81
          write(6,'(a)') 'INFO: loading the buffer atoms for adaptive QM/MM'  ! lam81

         call atommask( natom, nres, 0, ih(m04), ih(m06), &               ! lam81
          ix(i02), ih(m02), x(lcrd), buffermask, isqm )                   ! lam81

         qmmm_struct%buffer_nquant = sum(isqm(1:natom))                     ! lam81
         write(6,'(a,a,a,i5,a)') '     Mask ', buffermask(1:len_trim(buffermask)), &   ! lam81
         ' matches ',qmmm_struct%buffer_nquant,' atoms'                                ! lam81

         j = 0                           ! lam81
         do i=1,natom                    ! lam81
            if( isqm(i)>0 ) then         ! lam81
               j = j+1                   ! lam81
               buffer_iqmatoms(j) = i    ! lam81
            end if                       ! lam81
         end do                          ! lam81

        else                                                                 ! lam81
         write(6,'(a)') 'INFO: buffer atoms for adaptive QM/MM are not defined'    ! lam81
         if(qmmm_struct%r_buffer_out == 0) then                              ! lam81
          write(6,'(a)') 'WARNING: buffer radius was set to zero and no buffer atoms were given'  ! lam81
          write(6,'(a)') '         this may lead to unconverged QM forces'                        ! lam81
         end if                                                                                   ! lam81

         qmmm_struct%buffer_nquant = 0                             ! lam81

        end if
        write(6,'(a,i2)') '     buffer-charge: ', buffercharge     ! lam81
        qmmm_nml%buffercharge = buffercharge                       ! lam81

        if( ext_buffermask_subset /= '' ) then                                                             ! lam81
         write(6,'(a)') 'INFO: buffer subset was specified for the extended buffer region'                 ! lam81
         write(6,'(a)') 'INFO: loading buffer subset atoms'                                                ! lam81

         call atommask( natom, nres, 0, ih(m04), ih(m06), &                                                ! lam81
          ix(i02), ih(m02), x(lcrd), ext_buffermask_subset, isqm )                                         ! lam81

         qmmm_struct%buffer_nsubset = sum(isqm(1:natom))                                                   ! lam81
         write(6,'(a,a,a,i5,a)') '     Mask ', ext_buffermask_subset(1:len_trim(ext_coremask_subset)), &   ! lam81
         ' matches ',qmmm_struct%buffer_nsubset,' atoms'                                                   ! lam81
         if(qmmm_struct%buffer_nsubset == 0) then                                                            ! lam81
          write(6,'(a)') 'INFO: buffer subset is an empty set => all atoms can be in extended buffer region' ! lam81
         end if                                                                                              ! lam81

         j = 0                               ! lam81
         do i=1,natom                        ! lam81
            if( isqm(i)>0 ) then             ! lam81
               j = j+1                       ! lam81
               buffer_subsetatoms(j) = i     ! lam81
            end if                           ! lam81
         end do                              ! lam81

        else                                                                                        ! lam81
         write(6,'(a)') 'INFO: buffer subset was not specified for the extended buffer region'      ! lam81

         qmmm_struct%buffer_nsubset = 0                                                             ! lam81
        end if                                                                                      ! lam81

        write(6,'(a)') ''                                               ! lam81
        write(6,'(a)') 'CENTER atoms:'                                  ! lam81
        write(6,'(a)') '-------------'                                  ! lam81
        write(6,'(a)') ''                                               ! lam81
        if(centermask /= '') then                                       ! lam81
          write(6,'(a)') 'INFO: center atom list was specified'         ! lam81
          write(6,'(a)') 'INFO: loading the center atoms for adaptive QM/MM'  ! lam81

         call atommask( natom, nres, 0, ih(m04), ih(m06), &               ! lam81
         ix(i02), ih(m02), x(lcrd), centermask, isqm )                    ! lam81

         qmmm_struct%center_nsubset = sum(isqm(1:natom))                   ! lam81
         write(6,'(a,a,a,i5,a)') '     Mask ', centermask(1:len_trim(centermask)), &   ! lam81
         ' matches ',qmmm_struct%center_nsubset,' atoms'                               ! lam81
         if(qmmm_struct%center_nsubset == 0) then                                      ! lam81
          if(qmmm_struct%core_nquant > 0) then                                         ! lam81
           write(6,'(a)') 'INFO: center atom list is an empty set => center atom list = user defined core list' ! lam81
          else                                                                         ! lam81
           write(6,'(a)') 'INFO: center atom list is an empty set => center atom list = user defined qm list' ! lam81
          end if                                                                       ! lam81
         end if                                                                        ! lam81

         j = 0                               ! lam81
         do i=1,natom                        ! lam81
            if( isqm(i)>0 ) then             ! lam81
               j = j+1                       ! lam81
               center_subsetatoms(j) = i     ! lam81
            end if                           ! lam81
         end do

        else
         if(qmmm_struct%core_nquant > 0) then                                         ! lam81
          write(6,'(a)') 'INFO: center atom list is not defined => center atom list = user defined core list' ! lam81
         else                                                                         ! lam81
          write(6,'(a)') 'INFO: center atom list is not defined => center atom list = user defined qm list' ! lam81
         end if                                                                       ! lam81
         qmmm_struct%center_nsubset = 0                                               ! lam81
        end if                                                                        ! lam81

       end if                                                      ! lam81

       if(abfqmmm == 1) qmmm_struct%abfqmmm = 1       ! lam81
       if(abfqmmm == 1 .and. hot_spot == 1) then      ! lam81
         qmmm_struct%hot_spot = 1                     ! lam81
       end if                                         ! lam81

       deallocate(isqm, stat=ier)
       REQUIRE(ier==0)

    else  !  get the count from the input iqmatoms array

       do i = 1,max_quantum_atoms
          if( iqmatoms(i) == 0 ) exit
       end do
       qmmm_struct%nquant = i-1

    end if

   end if

   if(abfqmmm == 1) then                                                          ! lam81
    if(qmstep == 0) then                                                          ! lam81
     if(qmmm_struct%nquant > 0) then                                              ! lam81
      call validate_qm_atoms(iqmatoms,qmmm_struct%nquant,natom)                   ! lam81
      call int_legal_range('QMMM: (number of qm atoms) ', &                       ! lam81
        qmmm_struct%nquant, 1, max_quantum_atoms )                                ! lam81
      call qmsort(iqmatoms)                                                       ! lam81
      allocate(qmmm_struct%iqmatoms(qmmm_struct%nquant))                          ! lam81
      do i = 1, qmmm_struct%nquant                                                ! lam81
       qmmm_struct%iqmatoms(i)=iqmatoms(i)                                        ! lam81
      end do                                                                      ! lam81
     end if                                                                       ! lam81
     if(qmmm_struct%qm_nsubset > 0) then                                          ! lam81
      call validate_qm_atoms(qm_subsetatoms,qmmm_struct%qm_nsubset,natom)         ! lam81
      call int_legal_range('QMMM: (number of qm subset atoms) ', &                ! lam81
        qmmm_struct%qm_nsubset, 1, natom )                                        ! lam81
      call qmsort(qm_subsetatoms)                                                 ! lam81
      allocate(qmmm_struct%qm_subsetatoms(qmmm_struct%qm_nsubset))                ! lam81
      do i = 1, qmmm_struct%qm_nsubset                                            ! lam81
       qmmm_struct%qm_subsetatoms(i)=qm_subsetatoms(i)                            ! lam81
      end do                                                                      ! lam81
     end if                                                                       ! lam81
     if(qmmm_struct%core_nquant > 0) then                                         ! lam81
      call validate_qm_atoms(core_iqmatoms,qmmm_struct%core_nquant,natom)         ! lam81
      call int_legal_range('QMMM: (number of core atoms) ', &                     ! lam81
        qmmm_struct%core_nquant, 1, max_quantum_atoms )                           ! lam81
      allocate(qmmm_struct%core_iqmatoms(qmmm_struct%core_nquant))                ! lam81
      do i = 1, qmmm_struct%core_nquant                                           ! lam81
       qmmm_struct%core_iqmatoms(i)=core_iqmatoms(i)                              ! lam81
      end do                                                                      ! lam81
     end if                                                                       ! lam81
     if(qmmm_struct%core_nsubset > 0) then                                        ! lam81
      call validate_qm_atoms(core_subsetatoms,qmmm_struct%core_nsubset,natom)     ! lam81
      call int_legal_range('QMMM: (number of core subset atoms) ', &              ! lam81
        qmmm_struct%core_nsubset, 1, natom )                                      ! lam81
      call qmsort(core_subsetatoms)                                               ! lam81
      allocate(qmmm_struct%core_subsetatoms(qmmm_struct%core_nsubset))            ! lam81
      do i = 1, qmmm_struct%core_nsubset                                          ! lam81
       qmmm_struct%core_subsetatoms(i)=core_subsetatoms(i)                        ! lam81
      end do                                                                      ! lam81
     end if                                                                       ! lam81
     if( (qmmm_struct%nquant == 0) .and. (qmmm_struct%core_nquant == 0) ) then    ! lam81
      write(6,*)                                                                  ! lam81
      write(6,*) 'ERROR: qmmask and coremask are empty sets!'                     ! lam81
      stop                                                                        ! lam81
     end if                                                                       ! lam81
     if(qmmm_struct%buffer_nquant > 0) then                                       ! lam81
      call validate_qm_atoms(buffer_iqmatoms,qmmm_struct%buffer_nquant,natom)     ! lam81
      call int_legal_range('QMMM: (number of buffer atoms) ', &                   ! lam81
        qmmm_struct%buffer_nquant, 1, max_quantum_atoms )                         ! lam81
      allocate(qmmm_struct%buffer_iqmatoms(qmmm_struct%buffer_nquant))            ! lam81
      do i = 1, qmmm_struct%buffer_nquant                                         ! lam81
       qmmm_struct%buffer_iqmatoms(i)=buffer_iqmatoms(i)                          ! lam81
      end do                                                                      ! lam81
     end if                                                                       ! lam81
     if(qmmm_struct%buffer_nsubset > 0) then                                      ! lam81
      call validate_qm_atoms(buffer_subsetatoms,qmmm_struct%buffer_nsubset,natom) ! lam81
      call int_legal_range('QMMM: (number of buffer subset atoms) ', &            ! lam81
        qmmm_struct%buffer_nsubset, 1, natom )                                    ! lam81
      call qmsort(buffer_subsetatoms)                                             ! lam81
      allocate(qmmm_struct%buffer_subsetatoms(qmmm_struct%buffer_nsubset))        ! lam81
      do i = 1, qmmm_struct%buffer_nsubset                                        ! lam81
       qmmm_struct%buffer_subsetatoms(i)=buffer_subsetatoms(i)                    ! lam81
      end do                                                                      ! lam81
     end if                                                                       ! lam81
     if(qmmm_struct%center_nsubset > 0) then                                      ! lam81
      call validate_qm_atoms(center_subsetatoms,qmmm_struct%center_nsubset,natom) ! lam81
      call int_legal_range('QMMM: (number of center subset atoms) ', &            ! lam81
        qmmm_struct%center_nsubset, 1, natom )                                    ! lam81
      call qmsort(center_subsetatoms)                                             ! lam81
      allocate(qmmm_struct%center_subsetatoms(qmmm_struct%center_nsubset))        ! lam81
      do i = 1, qmmm_struct%center_nsubset                                        ! lam81
       qmmm_struct%center_subsetatoms(i)=center_subsetatoms(i)                    ! lam81
      end do                                                                      ! lam81
     end if                                                                       ! lam81
        ! the regions must be disjoint sets                                       ! lam81
     do i = 1, qmmm_struct%nquant                                                 ! lam81
      do j = 1, qmmm_struct%core_nquant                                           ! lam81
       if(qmmm_struct%iqmatoms(i) == qmmm_struct%core_iqmatoms(j)) then           ! lam81
        write(6,*)                                                                ! lam81
        write(6,*) 'ERROR: qmmask and coremask must be disjoint sets!'            ! lam81
        write(6,*) 'ERROR: atom number', qmmm_struct%iqmatoms(i),' is common!'    ! lam81
        stop                                                                      ! lam81
       end if                                                                     ! lam81
      end do                                                                      ! lam81
      do j = 1, qmmm_struct%buffer_nquant                                         ! lam81
       if(qmmm_struct%iqmatoms(i) == qmmm_struct%buffer_iqmatoms(j)) then         ! lam81
        write(6,*)                                                                ! lam81
        write(6,*) 'ERROR: qmmask and buffermask must be disjoint sets!'          ! lam81
        write(6,*) 'ERROR: atom number', qmmm_struct%iqmatoms(i),' is common!'    ! lam81
        stop                                                                      ! lam81
       end if                                                                     ! lam81
      end do                                                                      ! lam81
     end do                                                                       ! lam81
     do i = 1, qmmm_struct%core_nquant                                            ! lam81
      do j = 1, qmmm_struct%buffer_nquant                                         ! lam81
       if(qmmm_struct%core_iqmatoms(i) == qmmm_struct%buffer_iqmatoms(j)) then    ! lam81
        write(6,*)                                                                ! lam81
        write(6,*) 'ERROR: coremask and buffermask must be disjoint sets!'        ! lam81
        write(6,*) 'ERROR: atom number', qmmm_struct%core_iqmatoms(i),' is common!' ! lam81
        stop                                                                      ! lam81
       end if                                                                     ! lam81
      end do                                                                      ! lam81
     end do                                                                       ! lam81
     return                                                                       ! lam81
    else                                                                          ! lam81
     if(associated(qmmm_struct%iqmatoms)) deallocate(qmmm_struct%iqmatoms)        ! lam81
     if(associated(qmmm_struct%core_iqmatoms)) deallocate(qmmm_struct%core_iqmatoms)     !lam81
     if(associated(qmmm_struct%buffer_iqmatoms)) deallocate(qmmm_struct%buffer_iqmatoms) !lam81
    end if                                                                        ! lam81
   end if                                                                         ! lam81

!Initialize nlink to 0
   qmmm_struct%nlink = 0
   qmmm_struct%nquant_nlink = qmmm_struct%nquant

#endif
   if (.not. qmmm_nml%qmtheory%SEBOMD) then ! no SEBOMD
! Test to see if QM atom selection is legal.
   call validate_qm_atoms(iqmatoms,qmmm_struct%nquant,natom)

!  check we don't bust our statically allocated max_quantum_atoms
   call int_legal_range('QMMM: (number of quantum atoms) ', &
      qmmm_struct%nquant, 1, max_quantum_atoms )

   call qmsort(iqmatoms) !ensure the list of qm atoms is sorted numerically
   endif ! no SEBOMD

! --- Variable QM solvent region - has to be very early here because we
!     will be changing nquant and iqmatoms.  ---


   qmmm_nml%vsolv = vsolv
   if (qmmm_nml%vsolv > 0) then
#ifdef SQM
      write(6,*) 'SQM does not support the use of nearest_qm_solvent.'
      call mexit(6,1)
#else
      call read_vsolv_nml(qmmm_vsolv, nres)
      !We need to work out how many atoms nearest_qm_solvent * natoms_per_solvent_residue equals.
      !We then need to find the nearest atoms and update nquant and iqmatoms respectively.
      call qmmm_vsolv_setup(qmmm_struct%nquant, max_quantum_atoms, iqmatoms, &
                            nres, ih(m02), ix(i02),ix(i70), natom)
      ! check again that we don't bust our statically allocated max_quantum_atoms
      call int_legal_range('QMMM: (number of quantum atoms) ', &
           qmmm_struct%nquant, 1, max_quantum_atoms )
#endif
   end if
! --- End Variable QM water region ---

   call float_legal_range('QMMM: (QM-MM Cutoff) ', qmcut,0.0D0,1.0D30)
   call int_legal_range('QMMM: (variable solvent VSOLV) ', vsolv,0,3)
   call int_legal_range('QMMM: (QM GB Method) ', qmgb,0,3)
   call int_legal_range('QMMM: (QM-QM Derivatives) ', qmqmdx,1,2)
   call int_legal_range('QMMM: (Verbosity) ', verbosity,0,5)
   call int_legal_range('QMMM: (Max SCF Iterations) ', itrmax,1,10000000)
   call int_legal_range('QMMM: (Shake on QM atoms) ', qmshake,0,1)
   call int_legal_range('QMMM: (Density Matrix Convergence) ', tight_p_conv,0,1)
   call float_legal_range('QMMM: (SCF Convergence) ', scfconv,1.0D-16,1.0D0)
   !+TJG 01/26/2010
   call float_legal_range('QMMM: (Error Matrix Convergence) ', errconv,1.0D-16,1.0D0)
   call int_legal_range('QMMM: (Max num matrices in DIIS extrapolation) ', ndiis_matrices,1,20)
   call int_legal_range('QMMM: (Max num of DIIS attempts) ', ndiis_attempts,0,1000)
   !-TJG 01/26/2010
   call int_legal_range('QMMM: (PRINT CHARGES) ', printcharges,0,1)
   call int_legal_range('QMMM: (PRINT BONDORDERS) ',printbondorders,0,1)
   call int_legal_range('QMMM: (PRINT QM/Dipole) ', printdipole,0,2)
   if (qmmm_nml%qmtheory%EXTERN) then                                                  
     !AWG: Allow any spin multiplicity for external QM programs
     call int_legal_range('QMMM: (Spin multiplicity) ', spin, 1,100)
   else
     call int_legal_range('QMMM: (Spin multiplicity) ', spin,1,1)
     !RCW: Currently limit spin state to singlets only since the code for spin>1 does not exist / work at present.
     !     WARNING - IF WE LATER ALLOW SPIN>1 qm2_densit will need updating.
   end if
   call int_legal_range('QMMM: (Peptide Correction) ',peptide_corr,0,1)
   call int_legal_range('QMMM: (QM-MM RIJ in Core) ',qmmmrij_incore,0,1)
   call int_legal_range('QMMM: (QM-QM E-Rep in Core) ',qmqm_erep_incore,0,1)
   call int_legal_range('QMMM: (Link Atomic Number) ',lnk_atomic_no,1,numberElements)
   call int_legal_range('QMMM: (QM-MM Link Method) ',lnk_method,1,2)
   call int_legal_range('QMMM: (Pseudo Diag) ',pseudo_diag,0,1)
   call int_legal_range('QMMM: (QM Ewald) ',qm_ewald,0,2)
   call int_legal_range('QMMM: (QM PME) ',qm_pme,0,1)
   call int_legal_range('QMMM: (QM Ewald kmaxqx) ',kmaxqx,1,99999999)
   call int_legal_range('QMMM: (QM Ewald kmaxqy) ',kmaxqy,1,99999999)
   call int_legal_range('QMMM: (QM Ewald kmaxqz) ',kmaxqz,1,99999999)
   call int_legal_range('QMMM: (QM Ewald ksqmaxq) ',ksqmaxq,1,kmaxqx*kmaxqy*kmaxqz)
   call int_legal_range('QMMM: (QM-MM qmmm_int) ',qmmm_int,0,5)
   call int_legal_range('QMMM: (QM-MM adjust_q) ',adjust_q,0,2)
   call int_legal_range('QMMM: (QM-MM diag_routine) ',diag_routine,0,7)
#ifdef OPENMP
   call int_legal_range('QMMM: (QM-MM qmmm_omp_max_threads) ',qmmm_omp_max_threads,1,32)
#endif
   call int_legal_range('QMMM: (QM-MM density_predict) ',density_predict,0,1)
   call int_legal_range('QMMM: (QM-MM fock_predict) ',fock_predict,0,1)
   call float_legal_range('QMMM: (Pseudo Diag Criteria) ',pseudo_diag_criteria,1.0D-12,1.0D0)
   call int_legal_range('QMMM: (QM-MM qmmm_switch) ',qmmm_switch,0,1)
   call float_legal_range('QMMM: (QM-MM r_switch_lo) ',r_switch_lo,0.0D0,1.0D30)
   call float_legal_range('QMMM: (QM-MM r_switch_hi) ',r_switch_hi,0.0D0,1.0D30)
   if (lnk_dis>0.0d0) then
     !if lnk_dis is less than 0.0d0 then the link atom is just placed on top of
     !the MM link pair atom.
     call float_legal_range('QMMM: (Link Atom Distance) ',lnk_dis,0.7D0,4.0D0)
   endif

!! GMS
   call float_legal_range('QMMM: (QM-MM chg_lambda)'     , chg_lambda     , 0.0D0 , 1.0D0  )
   qmmm_nml%chg_lambda  = chg_lambda
!! DFTB
   call int_legal_range(  'QMMM: (QM-MM dftb_maxiter ) ' , dftb_maxiter   , 1     , 10000  )
   call int_legal_range(  'QMMM: (QM-MM dftb_disper) '   , dftb_disper    , 0     , 1      )
   call int_legal_range(  'QMMM: (QM-MM dftb_chg   ) '   , dftb_chg       , 0     , 1      )
   call float_legal_range('QMMM: (QM-MM dftb_telec)'     , dftb_telec     , 0.0D0 , 1.0D4  )

   if (dftb_3rd_order /= 'NONE') then
      call check_dftb_3rd_order(dftb_3rd_order)
      qmmm_nml%dftb_3rd_order = dftb_3rd_order
   endif

   qmmm_nml%dftb_maxiter   = dftb_maxiter
   qmmm_nml%dftb_disper      = dftb_disper
   qmmm_nml%dftb_chg         = dftb_chg
   qmmm_nml%dftb_telec       = dftb_telec
   qmmm_nml%dftb_telec_step  = dftb_telec_step

   if (dftb_chg > 0) printcharges=1

   qmmm_nml%qmcut = qmcut
   qmmm_nml%qmcut2 = qmcut*qmcut
   qmmm_nml%lnk_dis = lnk_dis
   qmmm_nml%lnk_atomic_no = lnk_atomic_no
   qmmm_nml%lnk_method = lnk_method

!DFTB Limitations - current things not supported in DFTB
!These are silent limitations that are non fatal - just to
!avoid problems with the default. Fatal errors are handled
!later in this routine.
   if (qmmm_nml%qmtheory%DFTB) then
     qmmmrij_incore=0
     qmqm_erep_incore=0
     pseudo_diag=0
     qmqmdx=1
   end if

!Divcon limitations - current things not supported in sander.DIVCON
!These are silent changes to remove defaults. Fatal errors are handled
!later in this routine.
   if (idc /= 0) then
     qmmmrij_incore=0
     qmqm_erep_incore=0
     pseudo_diag=0
   end if

   !qmgb values:
   ! 0 - do GB but leave QM charges as zero and add nothing to Fock matrix. This is like a vacuum
   !     QM molecule in a solvated MM system.
   ! 1 - do GB using the prmtop fixed resp charges for the GB calculation.
   ! 2 - do GB using Mulliken charges that are consistent with the GB field by modifying the fock
   !     matrix at every SCF step. (default)
   ! 3 - do GB using QM gas phase Mulliken charges - This is really a debugging option since the charges
   !     will not be consistent with the GB field since the fock matrix is not modified. This similarly
   !     means that the gradients will not be accurate. A warning will be printed at every QM call if this
   !     option is selected.

   !Make sure igb in &cntrl namelist is compatible with qmgb setting.
   if (igb==0 .or. igb==6) then
      !no qmgb available
      qmgb = 0
   end if
   !Print warning about qmgb being for debugging only.
   if (qmgb==3) then
     write(6,*) "QMMM: ------------------------------ WARNING --------------------------------"
     write(6,*) "QMMM: qmgb = 3 is designed for debugging purposes only. It gives GB"
     write(6,*) "QMMM:          energies based on gas phase QM Mulliken charges. These charges"
     write(6,*) "QMMM:          are NOT consistent with the GB field felt by the QM region and"
     write(6,*) "QMMM:          so any gradients calculated using this approach will"
     write(6,*) "QMMM:          NOT BE ACCURATE."
     write(6,*) "QMMM:          This option is really designed for:"
     write(6,*) "QMMM:                SINGLE POINT ENERGY EVALUATIONS ONLY"
     write(6,*) "QMMM: ------------------------------ WARNING --------------------------------"
   end if
   qmmm_nml%qmgb = qmgb

   qmmm_struct%AM1_OR_PM3 = (qmmm_nml%qmtheory%AM1 .or. qmmm_nml%qmtheory%AM1D .or. qmmm_nml%qmtheory%PM3 &
                             .OR. qmmm_nml%qmtheory%PDDGPM3 .OR. qmmm_nml%qmtheory%PM3CARB1 .OR. &
                             qmmm_nml%qmtheory%RM1 .OR. qmmm_nml%qmtheory%PDDGPM3_08 .OR. qmmm_nml%qmtheory%PM3ZNB)
   qmmm_struct%PDDG_IN_USE = (qmmm_nml%qmtheory%PDDGPM3 .OR. qmmm_nml%qmtheory%PDDGMNDO  &
                             .OR. qmmm_nml%qmtheory%PDDGPM3_08 )
   if (qmmm_int == 3 .OR. qmmm_int == 4) then
      if (qmmm_nml%qmtheory%PM3) then
         qmmm_struct%PM3MMX_INTERFACE = .true.
      else
         qmmm_struct%PM3MMX_INTERFACE = .false.
         call sander_bomb('read_qmmm_namelist','qmmm_int == 3/4 (Modified QM-MM interface) but qm_theory /= PM3.', &
                       'qmmm_int == 3/4 is only supported along with PM3 Hamiltonian.')
      end if
   else
      qmmm_struct%PM3MMX_INTERFACE = .false.
   end if
   qmmm_nml%qmcharge = qmcharge
#ifndef SQM
   if(abfqmmm == 1) then                   ! lam81
    qmmm_nml%qmcharge = abfqmcharge        ! lam81
   end if                                  ! lam81
#endif
   qmmm_nml%spin = spin
   qmmm_nml%verbosity = verbosity
   qmmm_nml%itrmax = itrmax
   qmmm_nml%qmshake = qmshake
   qmmm_nml%pseudo_diag_criteria = pseudo_diag_criteria

   ! Analytical or numerical integral derivatives for semiempirical methods:
   ! Note: for d orbitals only numerical derivatives are available
   !       Thus switch to numerical derivatives for MNDO/d and AM1/d
   !       For PM6 qmmm_nml%qmqm_analyt will be adjusted after parameters are
   !       read in qm2_load_params_and_allocate since we can do analytical
   !       derivatives if we don't have an element with d orbitals
   if (qmqmdx /= 1 .or. qmmm_nml%qmtheory%MNDOD .or. qmmm_nml%qmtheory%AM1D) then
      ! Do numerical QM-QM derivatives in qm2
      qmmm_nml%qmqm_analyt = .false. 
   else
      !Do analytical QM-QM dericatives in qm2
      qmmm_nml%qmqm_analyt = .true.  
   end if

   if (tight_p_conv /= 1) then
      ! Loose density matrix convergence (0.05*sqrt(SCFCRT))
      qmmm_nml%tight_p_conv = .false. 
   else
      ! Tight density matrix convergence (SCFCRT)
      qmmm_nml%tight_p_conv = .true.  
   end if

   !Write a warning about excessively tight convergence requests.
   if ( scfconv < 1.0D-12 ) then
     write(6,'(" QMMM: WARNING - SCF Conv = ",G8.2)') scfconv
     write(6,*) "QMMM:           There is a risk of convergence problems when the"
     write(6,*) "QMMM:           requested convergence is less that 1.0D-12 kcal/mol."
   end if
   qmmm_nml%scfconv = scfconv

   !How tight do we want the density convergence?
   if (qmmm_nml%tight_p_conv) then
      qmmm_nml%density_conv = qmmm_nml%scfconv
   else
      qmmm_nml%density_conv = 0.05D0 * sqrt(qmmm_nml%scfconv)
   end if

   !+TJG 01/26/2010
   qmmm_nml%errconv = errconv
   qmmm_nml%ndiis_matrices = ndiis_matrices
   qmmm_nml%ndiis_attempts = ndiis_attempts
   !-TJG 01/26/2010

   if ( printcharges /= 1) then
      qmmm_nml%printcharges=.false.
   else
      qmmm_nml%printcharges=.true.
   end if

   qmmm_nml%printdipole=printdipole

   qmmm_nml%print_eigenvalues = print_eigenvalues

   if ( printbondorders /= 1) then
      qmmm_nml%printbondorders=.false.
   else
      qmmm_nml%printbondorders=.true.
   end if
   if ( peptide_corr == 0) then
      qmmm_nml%peptide_corr = .false.
   else
      qmmm_nml%peptide_corr =  .true.
   end if 
   if ( qmmmrij_incore == 0 .or. qmmm_int==0 .or. qmmm_int == 5 ) then
      qmmm_nml%qmmmrij_incore = .false.
   else
      qmmm_nml%qmmmrij_incore = .true. !Only available with qmmm_int>1 and qmmm_int /= 5
   end if
   if ( qmqm_erep_incore == 0 .or. qmqmdx == 2 ) then
      qmmm_nml%qmqm_erep_incore = .false.
   else
      !Only available with analytical derivatives.
      qmmm_nml%qmqm_erep_incore = .true.
   end if
   if ( pseudo_diag == 1 ) then
     qmmm_nml%allow_pseudo_diag = .true.
   else
     qmmm_nml%allow_pseudo_diag = .false.
   end if
 
   qmmm_nml%qm_ewald = qm_ewald
   qmmm_nml%ksqmaxq = ksqmaxq
   qmmm_nml%kmaxqx = kmaxqx
   qmmm_nml%kmaxqy = kmaxqy
   qmmm_nml%kmaxqz = kmaxqz
   qmmm_nml%kappa = kappa
   !If ntb=0 or use_pme =0 then we can't do qm_ewald so overide what the user may
   !have put in the namelist and set the value to false.
   if (ntb==0 .or. use_pme==0) then
     qmmm_nml%qm_ewald = 0
     qmmm_nml%qm_pme = .false.
   end if
   if (qmmm_nml%qm_ewald>0 .and. qm_pme>0) then
     qmmm_nml%qm_pme=.true.
   else
     qmmm_nml%qm_pme=.false.
   end if

   if ( writepdb == 0 ) then
     qmmm_nml%writepdb=.false.
   else
     qmmm_nml%writepdb=.true.
   end if

   qmmm_nml%adjust_q = adjust_q

   qmmm_nml%qmmm_int = qmmm_int

   if (qmmm_nml%qmmm_int == 5) then
      ! Mechanical embedding
      ! Do not use QM and QM/MM Ewald or PME
      write(6,'(a)') 'QMMM: Mechanical embedding in use'
      if ( qmmm_nml%qm_pme ) then
         write(6,'(a)') 'QMMM: WARNING'
         write(6,'(a)') 'QMMM: Switching off QM PME'
         qmmm_nml%qm_pme = .false.
      end if
      if ( qmmm_nml%qm_ewald > 0) then
         write(6,'(a)') 'QMMM: WARNING'
         write(6,'(a)') 'QMMM: Switching off QM Ewald'
         qmmm_nml%qm_ewald = 0
      end if
      ! Prevent adjust_q, there is nothing to adjust
      qmmm_nml%adjust_q = 0
      ! Set qmcut as zero
      qmmm_nml%qmcut = 0.1d0
      qmmm_nml%qmcut2 = qmmm_nml%qmcut * qmmm_nml%qmcut
   end if

   qmmm_nml%idc = idc
   qmmm_nml%divpb = divpb

   if (qmmm_switch == 1) then
     qmmm_nml%qmmm_switch = .true.
   else
     qmmm_nml%qmmm_switch = .false.
   end if
   qmmm_nml%r_switch_lo = r_switch_lo
   qmmm_nml%r_switch_hi = r_switch_hi
   if ( (qmmm_nml%r_switch_hi - qmmm_nml%r_switch_lo) < 0.0D0 ) then
     call sander_bomb('read_qmmm_namelist', &
       & 'r_switch_hi is smaller than r_switch_lo!', &
       & 'please try a different set.') 
   end if

!Setup some specific calculation flags that depend on namelist variables.
!Need to make sure these get copied to other threads in an MPI run.

   !Will we be calculating the Mulliken charges on every SCF iteration?
   !Default is no. Will be set to true in a bit if certain options, such as qm_ewald
   !require it.
#ifdef SQM
   qm2_struct%calc_mchg_scf = .true.
#else
   if (qmmm_nml%qm_ewald==1 .or. qmmm_nml%qmgb>1) then
     !We will be needing the mulliken charges on every SCF iteration
     qm2_struct%calc_mchg_scf = .true.
   else
     qm2_struct%calc_mchg_scf = .false.
   end if
#endif

   !DFTB Calculates Mulliken charges anyway so we might as well store them in the correct place.
   if (qmmm_nml%qmtheory%DFTB) qm2_struct%calc_mchg_scf = .true.

!At this point we know nquant and natom so we can allocate our arrays that depend on nquant or natom
!Note if this is a LES run qmmm_struct%nquant is nqaunt
!Note non master mpi threads need to call this allocation routine manually themselves.
   qmmm_nml%nquant = qmmm_struct%nquant
   qmmm_struct%natom = natom
   call allocate_qmmm( qmmm_nml, qmmm_struct, natom )

#ifdef SQM
   !qmmm_struct%iqm_atomic_numbers(1:natom) = atnum(1:natom)
   !qmmm_nml%iqmatoms(1:natom) = iqmatoms(1:natom)
   qmmm_struct%iqm_atomic_numbers(1:qmmm_struct%nquant_nlink) = atnum(1:qmmm_struct%nquant_nlink)
   qmmm_nml%iqmatoms(1:qmmm_struct%nquant_nlink) = iqmatoms(1:qmmm_struct%nquant_nlink)
   if (ncharge_in > 0) then
      qmmm_struct%qm_xcrd = 0.0D0
      j=0
      do i=1,qmmm_struct%qm_mm_pairs
         !write(0,*) "MM charge ", i
         qmmm_struct%qm_xcrd(1,i) = excharge(j+1)
         qmmm_struct%qm_xcrd(2,i) = excharge(j+2)
         qmmm_struct%qm_xcrd(3,i) = excharge(j+3)
         qmmm_struct%qm_xcrd(4,i) = excharge(j+4)
         j=j+4
      end do

      if (qmmm_struct%PM3MMX_INTERFACE) then
         allocate(qmmm_struct%qm_mm_pair_atom_numbers(qmmm_struct%qm_mm_pairs), stat=ier)
         REQUIRE(ier == 0)
         do i = 1, qmmm_struct%qm_mm_pairs
            qmmm_struct%qm_mm_pair_atom_numbers(i) = chgatnum(i)
         end do
      end if
   end if
#else
#   ifndef NO_SANDER_DIVCON
   !DIVCON SPECIFIC STUFF
   if(qmmm_nml%divpb == 1)then
      !the +100 is for link atoms, needs some way to figure out actual number of link atoms
      allocate(qmmm_div%all_atom_numbers(natom+100), stat=ier)
      REQUIRE(ier == 0)

      do i=1,natom
         call get_atomic_number(ih(m04+i-1),x(lmass+i-1),qmmm_div%all_atom_numbers(i))
      enddo
   endif
   !END DIVCON SPECIFIC STUFF
#   endif

   do i = 1, qmmm_struct%nquant
      qmmm_nml%iqmatoms(i) = iqmatoms(i)
      !Get the atomic numbers (used to be done in rdparm2...)
      j = iqmatoms(i)
      call get_atomic_number( ih(m04+j-1),x(lmass+j-1),qmmm_struct%iqm_atomic_numbers(i) )
   end do

#endif

   ! From now on the code uses qmmm_struct%iqmatoms
   ! which will be extended to contain link atom info
   qmmm_struct%iqmatoms(:) = qmmm_nml%iqmatoms(:)

   ! Now we have a list of atom numbers for QM atoms we can build a true false (natom long) list
   ! specifying what the quantum atoms are. Useful for doing quick .OR. operations against other
   ! lists.
   qmmm_struct%atom_mask = .false. !Note, sets entire natom long array to false

   do i = 1, qmmm_struct%nquant
     qmmm_struct%atom_mask(qmmm_nml%iqmatoms(i)) = .true.
   end do

   qmmm_nml%diag_routine = diag_routine
#ifdef OPENMP
   qmmm_nml%qmmm_omp_max_threads = qmmm_omp_max_threads

   !For the time being the number of threads to use for diag and pdiag
   !routines is set to max_threads - later this will be optimized if
   !diag_routine=0.
   qmmm_omp%diag_threads = qmmm_omp_max_threads
   qmmm_omp%pdiag_threads = qmmm_omp_max_threads
#endif
   qmmm_nml%density_predict = density_predict
   qmmm_nml%fock_predict = fock_predict
   qmmm_nml%fockp_d1 = fockp_d1
   qmmm_nml%fockp_d2 = fockp_d2
   qmmm_nml%fockp_d3 = fockp_d3
   qmmm_nml%fockp_d4 = fockp_d4
   qmmm_nml%damp = damp
   qmmm_nml%vshift = vshift

! --- CHECK FOR LIMITATIONS ---

  !--- Mechanical embedding is not supported for GB at the moment. ---
  if ( (igb > 0) .and. (qmmm_nml%qmmm_int == 5) ) then
     call sander_bomb('read_qmmm_nm_and_alloc','Mechanical embedding currently not supported with GB models.', &
                      'Cannot have igb > 0 and qmmm_int = 5')
  end if

  !--- You cannot mix Fock prediction with density prediction. ---
  if (qmmm_nml%fock_predict > 0 .and. qmmm_nml%density_predict > 0) then
      call sander_bomb('read_qmmm_nm_and_alloc','Fock matrix and Density matrix prediction are mutually exclusive.', &
                       'Cannot have fock_predict > 0 and density_predict > 0')
  end if

  !--- For Fock prediction the 4 pre-factors must sum to 1.0d0 ---
  if (qmmm_nml%fock_predict == 1) then
    if (abs(1.0d0-(qmmm_nml%fockp_d1 + qmmm_nml%fockp_d2 + qmmm_nml%fockp_d3 + qmmm_nml%fockp_d4)) > 1.0d-6) then
       write(6,*) 'QMMM: Failure, fockp_d1 to d4 must sum to 1.0d0 - current sum is', &
                  (qmmm_nml%fockp_d1 + qmmm_nml%fockp_d2 + qmmm_nml%fockp_d3 + qmmm_nml%fockp_d4)
       call sander_bomb('read_qmmm_nm_and_alloc','Fock matrix prediction coefficients do not sum to 1.0.', &
                         'adjust fockp_d1 to fockp_d4 so that they sum to 1.0.')
    end if
  end if

  !--- You cannot use variable solvent with GB calculations.
  if (qmmm_nml%qmgb > 0 .and. qmmm_nml%vsolv > 0) then
      call sander_bomb('read_qmmm_nm_and_alloc','Nearest QM solvent and qmgb are mutually exclusive.', &
                       'Cannot have nearest_qm_solvent > 0 and qmgb > 0')
  end if

  !--- DIVCON LIMITATIONS ---
  !Divcon currently only works with gas phase simulations.
  if (qmmm_nml%idc>0) then
    if (ntb /= 0) then
      !This covers qmewald, qm_pme as well.
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but periodic boundaries are in use.', &
                       'Periodic boundaries are currently only supported when idc == 0')
    end if
    if (qmmm_nml%qmgb/=0) then
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but qmgb/=0. QMMM GB is currently', &
                       'only supported when idc == 0')
    end if
    if (qmmm_nml%peptide_corr) then
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but peptide_corr /= 0.', &
                       'Peptide correction for divcon is handled by the divcon.in file. Set peptide_corr = 0 to proceed.')
    end if
    if (qmmm_nml%printcharges) then
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but printcharges /= 0.', &
                       'Printcharges is not available with idc > 0.')
    end if
    if (qmmm_nml%vsolv > 0) then
        call sander_bomb('read_qmmm_nm_and_alloc','Nearest QM solvent is not available with idc>0.', &
                         'Cannot have nearest_qm_solvent > 0 and idc > 0')
    end if
#ifdef MPI
    !No support for parallel Divcon
    write(6,*) 'Divcon capability (idc>0) can only run in serial mode for now'
    call mexit(6,1)
#endif
    !Write a warning about qm_theory being ignored with Divcon.
    write(6,'("|QMMM: WARNING DIVCON IN USE")')
    write(6,'("|QMMM: qm_theory IS IGNORED WHEN USING DIVCON - QM HAMILTONIAN MUST BE SELECTED")')
    write(6,'("|QMMM: IN DIVCON.IN FILE.")')
  end if
  !--- END DIVCON LIMITATIONS ---

  !--- DFTB LIMITATIONS ---
  if (qmmm_nml%qmtheory%DFTB ) then
    if (qmmm_nml%peptide_corr) then
      call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=DFTB but peptide_corr /= 0.', &
                       'Peptide correction is not available, or required for  DFTB. Set peptide_corr = 0 to proceed.')
    end if
  end if
  !--- END DFTB LIMITATIONS ---

  !--- PM6 LIMITATIONS ---
  if (qmmm_nml%qmtheory%PM6 ) then
    if (qmmm_nml%peptide_corr) then
       ! AWG: Do not use peptide correction with PM6 since it has not been parametrized for PM6
       call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=PM6 but peptide_corr /= 0.', &
                       'Peptide correction is not available for PM6. Set peptide_corr = 0 to proceed.')
    end if
  end if
  !--- END DFTB LIMITATIONS ---

  !--- EXTERNAL INTERFACE LIMITATIONS ---
  ! ADF/GAMESS support works through an external interface. Currently there
  ! are a number of limitations.
  if (qmmm_nml%qmtheory%EXTERN) then
    ! 1) PME and EWALD are not supported with EXTERN.
    if (qmmm_nml%qm_ewald /= 0) then
      call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=EXTERN but qm_ewald /= 0.', &
                       'The external interface does not currently support EWALD or PME.')
    end if
    ! 2) GB is not currently supported with EXTERN.
    if (qmmm_nml%qmgb /= 0) then
      call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=EXTERN but qmgb /= 0.', &
                       'The external interface does not currently support Generalized Born.')
    end if
  end if
  !--- END EXTERNAL INTERFACE LIMITATIONS ---
 

! --- END CHECK FOR LIMITATIONS ---

  return

end subroutine read_qmmm_nm_and_alloc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

