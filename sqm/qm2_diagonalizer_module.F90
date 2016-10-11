#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"
module qm2_diagonalizer_module
! ----------------------------------------------------------------------
! PURPOSE: Module collecting subroutines for testing diagonalizer
!          performance and choosing appropriate diagonalizer
!          and performing required memory allocation
!
!          Contains code originally written by Ross Walker
! 
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Bugfixes and Optimization for AmberTools 15 by Ross Walker
! Date  : October 2011
! ----------------------------------------------------------------------

  implicit none

  private
  public :: qm2_diagonalizer_setup

contains

  subroutine qm2_diagonalizer_setup(diag_routine, allow_pseudo_diag, &
       verbosity, master, &
#ifdef MPI
       communicator, &
#endif
       qm2_struct, qmmm_scratch, &
       silence)

    use qmmm_module, only : qm2_structure, qmmm_scratch_structure, qmmm_struct  ! lam81

    implicit none

    integer, intent(inout) :: diag_routine
    logical, intent(in)    :: allow_pseudo_diag
    integer, intent(in)    :: verbosity
    logical, intent(in)    :: master
#ifdef MPI
    integer, intent(in)    :: communicator
#endif
    type(qm2_structure), intent(inout) :: qm2_struct
    type(qmmm_scratch_structure), intent(inout) :: qmmm_scratch
    logical, intent(in)    :: silence
    
#ifdef MPI
    include 'mpif.h'
#endif

    _REAL_ :: aloc_real_scr
    integer :: ier
    integer :: ilaenv !is a function
    integer :: ilaenv_blocksize

    ! *** EIGEN VECTORS ***
    allocate (qm2_struct%eigen_vectors(qm2_struct%norbs,qm2_struct%norbs),stat=ier)
    REQUIRE(ier==0)
    allocate (qm2_struct%eigen_values(qm2_struct%norbs),stat=ier)
    REQUIRE(ier==0)

    !*** FULL DIAGONALIZATIONS and PSEUDO DIAGONALIZATIONS ***
    if (master .and. qmmm_struct%abfqmmm /= 1 .and. .not. silence) then
       write(6,*)
       write(6,'(''| QMMM: *** Diagonalization Routine Information ***'')')
       if (allow_pseudo_diag) &
            write(6,'(''| QMMM: Pseudo diagonalizations are allowed.'')')
    end if

    !*** PSEUDO DIAGONALIZATIONS ***
    if (master .and. allow_pseudo_diag) then
       !only master needs to do the matrix diagonalizations
       allocate( qmmm_scratch%pdiag_scr_norbs_norbs(qm2_struct%norbs,qm2_struct%norbs), &
            qmmm_scratch%pdiag_scr_noccupied_norbs(qm2_struct%nopenclosed,qm2_struct%norbs), &
            qmmm_scratch%pdiag_vectmp1(qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
            qmmm_scratch%pdiag_vectmp2(qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
            qmmm_scratch%pdiag_vectmp3(qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
            qmmm_scratch%pdiag_vecjs(2,qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
            stat=ier )
       REQUIRE(ier == 0)
    end if

    !*** FULL DIAGONALIZATIONS ***
    
    !  0 = Automatically pick fastest routine.
    !  1 = Use internal diagonalization routine. (default)
    !  2 = Use lapack dspev.
    !  3 = Use lapack dspevd.
    !  4 = Use lapack dspevx.
    !  5 = Use lapack dsyev.
    !  6 = Use lapack dsyevd.
    !  7 = Use lapack dsyevr.
    !  8 = Use lapack dsyevx. (not currently implemented)

    if (diag_routine == 0) then

       ! automatic selection of diagonalization routine based on timings.
#ifdef NO_DETAILED_TIMINGS
       ! not available with NO_DETAILED_TIMINGS - since we need a functioning wallclock call.
       call sander_bomb('qm2_load_params_and_allocate',&
            'diag_routine=0 is not available when compiled with NO_DETAILED_TIMINGS', &
            'Please manually select a diagonalization routine.')
#endif

       if (master) then
          if(qmmm_struct%abfqmmm /= 1 .and. .not. silence) &
             write(6,'(''| QMMM: Auto diagonalization routine selection is in use.'')')
          if(qmmm_struct%abfqmmm /= 1 .and. .not. silence) &
             write(6,'(''| QMMM:'')')
          ! only master does diagonalization at present.
          call qm2_time_diag_routines(diag_routine, verbosity, silence)
       end if

#ifdef MPI
       ! have a barrier here to make all threads wait for timing to have been done. We will also
       ! broadcast the final chosen routine. Not strictly necessary at present but may help with
       ! debugging later on.
       call mpi_bcast(diag_routine, 1, MPI_INTEGER, 0, communicator, ier)
       ! Note we also need to broadcast the value of allow_pseudo_diag in case it changed.
       call mpi_bcast(allow_pseudo_diag, 1, mpi_logical, 0, communicator, ier) 
#endif

    else

       if (.not. silence .and. master .and. qmmm_struct%abfqmmm /= 1) &
            write(6,'(''| QMMM: Auto diagonalization routine selection is disabled.'')')

    end if


    ! Now we either have diag_routine set in qmmm_nml or was set by the auto selection routine
    ! above we are now ready to do all our allocation etc.

    if (diag_routine == 1) then

       ! Internal routine.
       qmmm_scratch%lapack_dc_int_scr_aloc=0
       qmmm_scratch%lapack_dc_real_scr_aloc=0
       if (master) then 
          ! Only master does diagonalisation
          if(qmmm_struct%abfqmmm /= 1 .and. .not. silence) &
              write(6,'(''| QMMM: Using internal diagonalization routine (diag_routine=1).'')')
          allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ier)
          REQUIRE(ier==0)
       end if

    else if (diag_routine == 2) then

       ! DSPEV
       ! dspev needs a real workspace array of size (1,3n) where n = norbs
       ! We simply use the dimensions 2 to 4 of mat_diag_workspace for this
       qmmm_scratch%lapack_dc_int_scr_aloc  = 0
       qmmm_scratch%lapack_dc_real_scr_aloc = 0
       if (master) then
          ! Only master does diagonalisation
          if (.not. silence .and. qmmm_struct%abfqmmm /= 1) &
              write(6,'(''| QMMM: Using dspev routine (diag_routine=2).'')')
          allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,4),stat=ier)
          REQUIRE(ier==0)
         end if

      else if (diag_routine == 3) then

         ! DSPEVD
         ! We will use a divide and conquer algorithm for the matrix diagonalisation.
         ! now allocate the scratch arrays.
         ! Divide and conquer algorithm needs a minimum of:
         !   3+5*norbs ints and 1+6*norbs+norbs**2 reals of scratch space
         ! But it also supports the option of calling it to request how much space it needs.
         ! We shall do this by default.
         if (master) then
            ! Only master does diagonalisation
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            if(.not. silence .and. qmmm_struct%abfqmmm /= 1) &
                write(6,'(''| QMMM: Using dspevd routine (diag_routine=3).'')')
            if (verbosity >= 2 .and. .not. silence .and. qmmm_struct%abfqmmm /= 1) then
                write(6,'(''| QMMM: Calling dspevd to query required scratch array sizes.'')')
            end if
            call dspevd('V','U',qm2_struct%norbs,qm2_struct%fock_matrix, &
            qmmm_scratch%mat_diag_workspace(1,1), &
            qm2_struct%eigen_vectors, qm2_struct%norbs, &
            aloc_real_scr, -1, &
            qmmm_scratch%lapack_dc_int_scr_aloc , -1, &
            ier)
            qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
            if (verbosity >= 2) then
               if(.not. silence .and. qmmm_struct%abfqmmm /= 1) &
                   write(6,'(''| QMMM: dspevd required REAL scratch = '',i12,'' REALS.'')') &
                    qmmm_scratch%lapack_dc_real_scr_aloc
               if(.not. silence .and. qmmm_struct%abfqmmm /= 1) &
                   write(6,'(''| QMMM: dspevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                    qmmm_scratch%lapack_dc_int_scr_aloc
            end if
!             qmmm_scratch%lapack_dc_real_scr_aloc=1+6*qm2_struct%norbs+(qm2_struct%norbs*qm2_struct%norbs)
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
!            qmmm_scratch%lapack_dc_int_scr_aloc=3+5*qm2_struct%norbs
            allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
            REQUIRE(ier==0)
         else
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            qmmm_scratch%lapack_dc_real_scr_aloc=0
         end if

      else if (diag_routine == 4) then

         ! DSPEVX
         if (master) then 
            ! Only master does diagonalisation
            if(.not. silence .and. qmmm_struct%abfqmmm /= 1) &
                write(6,'(''| QMMM: Using dspevx routine (diag_routine=4).'')')
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            qmmm_scratch%lapack_dc_real_scr_aloc = 8*qm2_struct%norbs
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
            ! Note for dspevx we pad lapack_dc_int_scr by an extra norbs so we can use the 
            ! first norbs for IFAIL.
            qmmm_scratch%lapack_dc_int_scr_aloc=6*qm2_struct%norbs
            allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
            REQUIRE(ier==0)
         else
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            qmmm_scratch%lapack_dc_real_scr_aloc=0
         end if

      else if (diag_routine == 5) then

         ! DSYEV  - unpacked diagonalizer
         ! DSYEV has complex requirements of scratch arrays. 
         ! But it also supports the option of calling it to request how much space it needs.
         ! We shall do this by default.
         if (master) then !Only master does diagonalisation
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            if(.not. silence .and. qmmm_struct%abfqmmm /= 1) &
                write(6,'(''| QMMM: Using dsyev routine (diag_routine=5).'')')
            if (verbosity >= 2 .and. .not. silence .and. qmmm_struct%abfqmmm /= 1) then
               write(6,'(''| QMMM: Calling dsyev to query required scratch array sizes.'')')
            end if
            call dsyev('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, &
                       qm2_struct%norbs, qmmm_scratch%mat_diag_workspace(1,1), &
                       aloc_real_scr, -1, ier)
            qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            if (verbosity >= 2 .and. .not. silence .and. qmmm_struct%abfqmmm /= 1) then
               write(6,'(''| QMMM: dsyev required REAL scratch = '',i12,'' REALS.'')') &
                        qmmm_scratch%lapack_dc_real_scr_aloc
            end if
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
          else
             qmmm_scratch%lapack_dc_int_scr_aloc=0
             qmmm_scratch%lapack_dc_real_scr_aloc=0
          end if

       else if (diag_routine == 6) then

          ! DSYEVD - unpacked divide and conquor diagonalizer
          ! DSYEVD has complex requirements of scratch arrays. 
          ! But it also supports the option of calling it to request how much space it needs.
          ! We shall do this by default.
          if (master) then
             ! Only master does diagonalisation
             allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
             REQUIRE(ier==0)
             if(.not. silence .and. qmmm_struct%abfqmmm /= 1) &
                 write(6,'(''| QMMM: Using dsyevd routine (diag_routine=6).'')')
             if (verbosity >= 2 .and. .not. silence .and. qmmm_struct%abfqmmm /= 1) then
               write(6,'(''| QMMM: Calling dsyevd to query required scratch array sizes.'')')
             end if
             call dsyevd('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, qm2_struct%norbs, &
                         qmmm_scratch%mat_diag_workspace(1,1), aloc_real_scr, -1, &
                         qmmm_scratch%lapack_dc_int_scr_aloc , -1, ier)
             qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
             if (verbosity >= 2 .and. .not. silence .and. qmmm_struct%abfqmmm /= 1) then
                write(6,'(''| QMMM: dsyevd required REAL scratch = '',i12,'' REALS.'')') &
                     qmmm_scratch%lapack_dc_real_scr_aloc
                write(6,'(''| QMMM: dsyevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                     qmmm_scratch%lapack_dc_int_scr_aloc
             end if
             allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
             REQUIRE(ier==0)
             allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
             REQUIRE(ier==0)
          else
             qmmm_scratch%lapack_dc_int_scr_aloc=0
             qmmm_scratch%lapack_dc_real_scr_aloc=0
          end if

       else if (diag_routine == 7) then

          ! DSYEVR - Reduction to tridiagonal form before diagonalization.
          ! DSYEVR has complex requirements of scratch arrays. 
          ! NOTE: for dsyevr we can't just put the matrix into the eigenvectors array as
          !       it uses the matrix for scratch space. We will assume we can use the
          !       pdiag_scr_norbs_norbs scratch - array which will get allocated if we
          !       are doing pseudo diags and/or we are choosing diag_routine == 7.
          ! NOTE2: The lapack_dc_int array is actually made 2xnorbs bigger than what is
          !        requested since the first 2xnorbs are used for ISUPPZ.
          ! But it also supports the option of calling it to request how much space it needs.
          ! We shall do this by default.
          if (master) then !Only master does diagonalisation
             allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
             REQUIRE(ier==0)
             if(.not. silence .and. qmmm_struct%abfqmmm /= 1) &
                 write(6,'(''| QMMM: Using dsyevr routine (diag_routine=7).'')')
!            if (verbosity >= 2) then
!              write(6,'(''| QMMM: Calling dsyevr to query required scratch array sizes.'')')
!            end if

! RCW 2008/01/22 - weird segfaults in qm2_hcore_qmqm when we call dsyevr to calculate the memory requirements
!                  so for the moment we will just use the default amounts.
            !The dimension of the array WORK.  LWORK >= max(1,26*N).
            !For optimal efficiency, LWORK >= (NB+6)*N,
            !where NB is the max of the blocksize for DSYTRD and DORMTR
            !returned by ILAENV.
            ilaenv_blocksize = ILAENV( 1, 'DSYTRD', 'U', qm2_struct%norbs, -1, -1, -1 )
            ilaenv_blocksize = max(ilaenv_blocksize, ILAENV( 1, 'DORMTR', 'U', qm2_struct%norbs, -1, -1, -1 ))
            !RCW: 2015/04/07 weird behavior with MKL - it seems that sometimes
            !(ilaenv_blocksize+6)*qm2_struct%norbs can come back as smaller than
            !26*norbs - this causes the routine to quit with an error. So for now
            !use the maximum of the blocksize calc and 26*norbs
            qmmm_scratch%lapack_dc_real_scr_aloc=max((ilaenv_blocksize+6)*qm2_struct%norbs,26*qm2_struct%norbs)

            !The dimension of the array IWORK.  LIWORK >= max(1,10*N).
            qmmm_scratch%lapack_dc_int_scr_aloc = 10*qm2_struct%norbs

!            abstol = 2.0d0 * dlamch('S')
!            call dsyevr('V','A','U',qm2_struct%norbs,qmmm_scratch%pdiag_scr_norbs_norbs,qm2_struct%norbs, &
!                 0.0d0, 0.0d0, 0, 0, abstol, int_dummy, qmmm_scratch%mat_diag_workspace(1,1), &
!                 qm2_struct%eigen_vectors, qm2_struct%norbs, qmmm_scratch%lapack_dc_int_scr_aloc, &
!                 aloc_real_scr, -1, qmmm_scratch%lapack_dc_int_scr_aloc, -1, ier)
!            qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)

            !for dsyevr we pad lapack_dc_int_scr with 2*norbs to allow use in ISUPPZ.
            qmmm_scratch%lapack_dc_int_scr_aloc = qmmm_scratch%lapack_dc_int_scr_aloc+2*qm2_struct%norbs
            if (verbosity >= 2 .and. .not. silence .and. qmmm_struct%abfqmmm /= 1) then
              write(6,'(''| QMMM: dsyevr required REAL scratch = '',i12,'' REALS.'')') &
                        qmmm_scratch%lapack_dc_real_scr_aloc
              write(6,'(''| QMMM: dsyevr required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                        qmmm_scratch%lapack_dc_int_scr_aloc
            end if     
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
            allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
            REQUIRE(ier==0)
            if (.not. allow_pseudo_diag) then
              !qmmm_scratch%pdiag_scr_norbs_norbs will not have been allocated so make sure
              !we allocate it.
              allocate(qmmm_scratch%pdiag_scr_norbs_norbs(qm2_struct%norbs,qm2_struct%norbs),stat=ier)
              REQUIRE(ier==0)
           end if
        else
           qmmm_scratch%lapack_dc_int_scr_aloc=0
           qmmm_scratch%lapack_dc_real_scr_aloc=0
        end if

     else

        call sander_bomb('qm2_load_params_and_allocate', &
             'METHOD NOT CURRENTLY IMPLEMENTED', &
             'Selected Diagonalization routine is not currently implemented.')

     end if !diag_routine

   end subroutine qm2_diagonalizer_setup


  subroutine qm2_time_diag_routines(diag_routine,verbosity,silence)

    ! This routine tries a series of different diagonalizers and returns the value
    ! of diag_method that is fastest as the first argument.
#ifdef OPENMP
    use qmmm_module, only : qmmm_scratch, qm2_struct, qmmm_nml, qmmm_omp
#else
    use qmmm_module, only : qmmm_scratch, qm2_struct, qmmm_nml
#endif
    use constants, only : zero

    implicit none

    integer, intent(out) :: diag_routine
    integer, intent(in) :: verbosity
    logical, intent(in) :: silence

    ! local variables
    _REAL_ :: abstol, smallsum, small
    _REAL_ :: dlamch  !is a function
    _REAL_ :: start_time, end_time, current_time, fastest_time
    _REAL_ :: aloc_real_scr
    integer :: ilaenv !is a function
    integer :: ilaenv_blocksize
    integer :: diag_iterations
    integer :: i, ier
    integer :: diag_routine_to_test

#ifdef OPENMP
    ! For the time being we just use the number of threads specified by
    ! the qmmm_max_omp_threads for openmp - later we will update this
    ! to actually test performance for each thread count.
    call omp_set_num_threads(qmmm_omp%diag_threads)
#endif

    abstol = 2.0d0 * dlamch('S') !tolerance for dspevr

    ! We need to do enough calls to generate decent statistics. The Fortran timers do not appear to
    ! be incredibly good so we will call the diagonalization a number of times based on the size of
    ! the matrix.

    if (qm2_struct%norbs <= 50) then
       diag_iterations = 20
    elseif (qm2_struct%norbs <= 100 ) then
       diag_iterations = 15
    elseif (qm2_struct%norbs <= 150 ) then
       diag_iterations = 10
    elseif (qm2_struct%norbs <= 200 ) then
       diag_iterations = 5
    elseif (qm2_struct%norbs <= 250 ) then
       diag_iterations = 4
    elseif (qm2_struct%norbs <= 300 ) then
       diag_iterations = 3
    elseif (qm2_struct%norbs <= 350 ) then
       diag_iterations = 2
    else
       diag_iterations = 1
    endif

    if (.not. silence) then
       write(6,'(''| QMMM: Timing diagonalization routines:'')')
       write(6,'(''| QMMM:                              norbs = '',i8)') qm2_struct%norbs
       write(6,'(''| QMMM:    diag iterations used for timing = '',i8)') diag_iterations
       write(6,'(''| QMMM:'')')
    end if

    ! --------------------
    ! 1) INTERNAL ROUTINE
    ! --------------------

    allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ier)
    REQUIRE(ier==0)

    diag_routine_to_test = 1
    current_time=0.0d0
    do i = 1,diag_iterations
       ! Need to rebuild the fock matrix each time since some diag routines can destroy it.
       call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
       call wallclock(start_time)
       call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs, &
            qm2_struct%eigen_vectors,qm2_struct%eigen_values,abstol)
       call wallclock(end_time)
       current_time = current_time + (end_time-start_time)
    end do

    if (.not. silence) &
       write(6,'(''| QMMM:              Internal diag routine = '',F8.6,'' seconds'')') current_time

    ! initially just assume that the internal routine is the fastest.
    fastest_time = current_time
    diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to
                                        !converge then diag_routine_to_test will have been set back to
                                        !the internal diagonalizer by the diag routine.

    deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
    REQUIRE(ier==0)

  ! --------------------
  ! 2) DSPEV ROUTINE
  ! --------------------
  ! DSPEV
  ! dspev needs a real workspace array of size (1,3n) where n = norbs
  ! We simply use the dimensions 2 to 4 of mat_diag_workspace for this
    allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,4),stat=ier)
    REQUIRE(ier==0)
  
    diag_routine_to_test = 2
    current_time=0.0d0
    do i = 1,diag_iterations
       !Need to rebuild the fock matrix each time since some diag routines can destroy it.
       call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
       call wallclock(start_time)
       call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs, &
            qm2_struct%eigen_vectors,qm2_struct%eigen_values,abstol)
       call wallclock(end_time)
       current_time = current_time + (end_time-start_time)
    end do

    if (.not. silence) &
       write(6,'(''| QMMM:                 Dspev diag routine = '',F8.6,'' seconds'')') current_time

    !check if this routine is faster than the ones tried so far.
    if (fastest_time >= current_time) then
       diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to
                                           !converge then diag_routine_to_test will have been set back to
                                           !the internal diagonalizer by the diag routine.
       fastest_time = current_time
    end if
  
    deallocate (qmmm_scratch%mat_diag_workspace, stat=ier)
    REQUIRE(ier==0)

    ! --------------------
    ! 3) DSPEVD ROUTINE
    ! --------------------
    if (verbosity >= 2 .and. .not. silence) then
       write(6,'(''| QMMM: Calling dspevd to query required scratch array sizes.'')')
    end if
    allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
    REQUIRE(ier==0)
    call dspevd('V','U',qm2_struct%norbs,qm2_struct%fock_matrix, &
         qmmm_scratch%mat_diag_workspace(1,1), &
         qm2_struct%eigen_vectors, qm2_struct%norbs, aloc_real_scr, -1, &
         qmmm_scratch%lapack_dc_int_scr_aloc , -1, ier)
    
    qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
    
    if (verbosity >= 2 .and. .not. silence) then
       write(6,'(''| QMMM: dspevd required REAL scratch = '',i12,'' REALS.'')') &
            qmmm_scratch%lapack_dc_real_scr_aloc
       write(6,'(''| QMMM: dspevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
            qmmm_scratch%lapack_dc_int_scr_aloc
    end if
    allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
    REQUIRE(ier==0)
    allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
    REQUIRE(ier==0)
    
    diag_routine_to_test = 3
    current_time=0.0d0
    do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
       call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
       call wallclock(start_time)
       call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs, &
            qm2_struct%eigen_vectors,qm2_struct%eigen_values,abstol)
       call wallclock(end_time)
       current_time = current_time + (end_time-start_time)
    end do

    if (.not. silence) &
       write(6,'(''| QMMM:                Dspevd diag routine = '',F8.6,'' seconds'')') current_time
    
    !check if this routine is faster than the ones tried so far.
    if (fastest_time >= current_time) then
       diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to
       !converge then diag_routine_to_test will have been set back to
       !the internal diagonalizer by the diag routine.
       fastest_time = current_time
    end if
    deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
    REQUIRE(ier==0)
    qmmm_scratch%lapack_dc_real_scr_aloc = 0
    qmmm_scratch%lapack_dc_int_scr_aloc = 0
    ! --------------------
    ! 4) DSPEVX ROUTINE
    ! --------------------
    ! Allocate mat diag workspace as 2 here since (x,2) is used for error reporting by dspevx
    allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,2),stat=ier)
    REQUIRE(ier==0)
    qmmm_scratch%lapack_dc_real_scr_aloc = 8*qm2_struct%norbs
    allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
    REQUIRE(ier==0)
    qmmm_scratch%lapack_dc_int_scr_aloc=6*qm2_struct%norbs
    allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
    REQUIRE(ier==0)
    diag_routine_to_test = 4
    current_time=0.0d0
    do i = 1,diag_iterations
       !Need to rebuild the fock matrix each time since some diag routines can destroy it.
       call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
       call wallclock(start_time)
       call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs, &
            qm2_struct%eigen_vectors,qm2_struct%eigen_values,abstol)
       call wallclock(end_time)
       current_time = current_time + (end_time-start_time)
    end do

    if (.not. silence) &
       write(6,'(''| QMMM:                Dspevx diag routine = '',F8.6,'' seconds'')') current_time

    !check if this routine is faster than the ones tried so far.
    if (fastest_time >= current_time) then
       diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to
                                           !converge then diag_routine_to_test will have been set back to
                                           !the internal diagonalizer by the diag routine.
       fastest_time = current_time
    end if
    deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
    REQUIRE(ier==0)
    qmmm_scratch%lapack_dc_real_scr_aloc = 0
    qmmm_scratch%lapack_dc_int_scr_aloc = 0
    
  ! --------------------
  ! 5) DSYEV ROUTINE
  ! --------------------
    if (verbosity >= 2 .and. .not. silence) then
       write(6,'(''| QMMM: Calling dsyev to query required scratch array sizes.'')')
    end if
    allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
    REQUIRE(ier==0)
    call dsyev('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, &
         qm2_struct%norbs, qmmm_scratch%mat_diag_workspace(1,1), &
         aloc_real_scr, -1, ier)
    qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
    if (verbosity >= 2 .and. .not. silence) then
       write(6,'(''| QMMM: dsyev required REAL scratch = '',i12,'' REALS.'')') &
            qmmm_scratch%lapack_dc_real_scr_aloc
    end if
    allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
    REQUIRE(ier==0)

    diag_routine_to_test = 5
    current_time=0.0d0
    do i = 1,diag_iterations
       !Need to rebuild the fock matrix each time since some diag routines can destroy it.
       call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
       call wallclock(start_time)
       call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs, &
            qm2_struct%eigen_vectors,qm2_struct%eigen_values,abstol)
       call wallclock(end_time)
       current_time = current_time + (end_time-start_time)
    end do

    if (.not. silence) &
       write(6,'(''| QMMM:                 Dsyev diag routine = '',F8.6,'' seconds'')') current_time

    !check if this routine is faster than the ones tried so far.
    if (fastest_time >= current_time) then
       diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to
                                           !converge then diag_routine_to_test will have been set back to
                                           !the internal diagonalizer by the diag routine.
       fastest_time = current_time
    end if
    deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
    REQUIRE(ier==0)
    qmmm_scratch%lapack_dc_real_scr_aloc = 0
    qmmm_scratch%lapack_dc_int_scr_aloc = 0
    
  ! --------------------
  ! 6) DSYEVD ROUTINE
  ! --------------------
    if (verbosity >= 2 .and. .not. silence) then
       write(6,'(''| QMMM: Calling dsyevd to query required scratch array sizes.'')')
    end if
    allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
    REQUIRE(ier==0)
    call dsyevd('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, qm2_struct%norbs, &
         qmmm_scratch%mat_diag_workspace(1,1), aloc_real_scr, -1, &
         qmmm_scratch%lapack_dc_int_scr_aloc , -1, ier)
    qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
    if (verbosity >= 2 .and. .not. silence) then
       write(6,'(''| QMMM: dsyevd required REAL scratch = '',i12,'' REALS.'')') &
            qmmm_scratch%lapack_dc_real_scr_aloc
       write(6,'(''| QMMM: dsyevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
            qmmm_scratch%lapack_dc_int_scr_aloc
    end if
    allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
    REQUIRE(ier==0)
    allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
    REQUIRE(ier==0)
    
    diag_routine_to_test = 6
    current_time=0.0d0
    do i = 1,diag_iterations
       !Need to rebuild the fock matrix each time since some diag routines can destroy it.
       call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
       call wallclock(start_time)
       call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs, &
            qm2_struct%eigen_vectors,qm2_struct%eigen_values,abstol)
       call wallclock(end_time)
       current_time = current_time + (end_time-start_time)
    end do

    if (.not. silence) &
       write(6,'(''| QMMM:                Dsyevd diag routine = '',F8.6,'' seconds'')') current_time
    
    !check if this routine is faster than the ones tried so far.
    if (fastest_time >= current_time) then
       diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to
                                           !converge then diag_routine_to_test will have been set back to
                                           !the internal diagonalizer by the diag routine.
       fastest_time = current_time
    end if
    deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
    REQUIRE(ier==0)
    qmmm_scratch%lapack_dc_real_scr_aloc = 0
    qmmm_scratch%lapack_dc_int_scr_aloc = 0
    
  ! --------------------
  ! 7) DSYEVR ROUTINE
  ! --------------------
    allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
    REQUIRE(ier==0)
    ! RCW 2008/01/22 - weird segfaults in qm2_hcore_qmqm when we call dsyevr to calculate the memory requirements
    !                  so for the moment we will just use the default amounts.
    !The dimension of the array WORK.  LWORK >= max(1,26*N).
    !For optimal efficiency, LWORK >= (NB+6)*N,
    !where NB is the max of the blocksize for DSYTRD and DORMTR
    !returned by ILAENV.
    ilaenv_blocksize = ILAENV( 1, 'DSYTRD', 'U', qm2_struct%norbs, -1, -1, -1 )
    ilaenv_blocksize = max(ilaenv_blocksize, ILAENV( 1, 'DORMTR', 'U', qm2_struct%norbs, -1, -1, -1 ))
    !RCW: 2015/04/07 weird behavior with MKL - it seems that sometimes
    !(ilaenv_blocksize+6)*qm2_struct%norbs can come back as smaller than
    !26*norbs - this causes the routine to quit with an error. So for now
    !use the maximum of the blocksize calc and 26*norbs
    qmmm_scratch%lapack_dc_real_scr_aloc=max((ilaenv_blocksize+6)*qm2_struct%norbs,26*qm2_struct%norbs)
    !The dimension of the array IWORK.  LIWORK >= max(1,10*N).
    qmmm_scratch%lapack_dc_int_scr_aloc = 10*qm2_struct%norbs
    !for dsyevr we pad lapack_dc_int_scr with 2*norbs to allow use in ISUPPZ.
    qmmm_scratch%lapack_dc_int_scr_aloc = qmmm_scratch%lapack_dc_int_scr_aloc+2*qm2_struct%norbs
    if (verbosity >= 2 .and. .not. silence) then
       write(6,'(''| QMMM: dsyevr required REAL scratch = '',i12,'' REALS.'')') &
            qmmm_scratch%lapack_dc_real_scr_aloc
       write(6,'(''| QMMM: dsyevr required INTEGER scratch = '',i12,'' INTEGERS.'')') &
            qmmm_scratch%lapack_dc_int_scr_aloc
    end if
    allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
    REQUIRE(ier==0)
    allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
    REQUIRE(ier==0)
    if (.not. qmmm_nml%allow_pseudo_diag) then
       !qmmm_scratch%pdiag_scr_norbs_norbs will not have been allocated so make sure
       !we allocate it.
       allocate(qmmm_scratch%pdiag_scr_norbs_norbs(qm2_struct%norbs,qm2_struct%norbs),stat=ier)
       REQUIRE(ier==0)
    end if
    
    diag_routine_to_test = 7
    current_time=0.0d0
    do i = 1,diag_iterations
       !Need to rebuild the fock matrix each time since some diag routines can destroy it.
       call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
       call wallclock(start_time)
       call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs, &
            qm2_struct%eigen_vectors,qm2_struct%eigen_values,abstol)
       call wallclock(end_time)
       current_time = current_time + (end_time-start_time)
    end do

    if (.not. silence) &
       write(6,'(''| QMMM:                Dsyevr diag routine = '',F8.6,'' seconds'')') current_time

    !check if this routine is faster than the ones tried so far.
    if (fastest_time >= current_time) then
       diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to
                                           !converge then diag_routine_to_test will have been set back to
                                           !the internal diagonalizer by the diag routine.
       fastest_time = current_time
    end if
    deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
    REQUIRE(ier==0)
    deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
    REQUIRE(ier==0)
    if (.not. qmmm_nml%allow_pseudo_diag) then
       deallocate(qmmm_scratch%pdiag_scr_norbs_norbs,stat=ier)
       REQUIRE(ier==0)
    end if
    
    qmmm_scratch%lapack_dc_real_scr_aloc = 0
    qmmm_scratch%lapack_dc_int_scr_aloc = 0
 
    if (.not. silence) write(6,'(''| QMMM:'')')

#ifdef OPENMP
    !For the time being we just use the number of threads specified by
    !the qmmm_max_omp_threads for openmp - later we will update this
    !to actually test performance for each thread count.
    call omp_set_num_threads(qmmm_omp%pdiag_threads)
#endif

    !-----------------------
    !PSEUDO DIAGONALIZATIONS
    !-----------------------
    !Finally we will test how long pseudo diagonalizations
    !take - this is not an exact science at all since the ratio
    !of pseudo diags to real diags does not remain constant. 
    !We will simply check how it compares like for like with
    !full diagonalizations. If the ful diagonalization is faster
    !than the pseudo diagonalization then we will turn off the
    !pseudo diagonalization. This will give faster execution but
    !may miss cases where it is still beneficial to turn off
    !the pseudo diag.
    if (qmmm_nml%allow_pseudo_diag) then
       allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
       REQUIRE(ier==0)
       qmmm_scratch%mat_diag_workspace(:,:) = 0.d0

       call qm2_smallest_number(small,smallsum)
      smallsum = max(10.0D0 * sqrt(smallsum),1.4000D-7)
      
      current_time=0.0d0
      do i = 1,diag_iterations
         !Need to rebuild the fock matrix each time since some diag routines can destroy it.
         call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
         call wallclock(start_time)
         call qm2_pseudo_diag(qm2_struct%fock_matrix,qm2_struct%eigen_vectors,qm2_struct%nopenclosed, &
              qmmm_scratch%mat_diag_workspace(1,1),qm2_struct%norbs,smallsum, &
              qmmm_scratch%pdiag_scr_norbs_norbs,qmmm_scratch%pdiag_scr_noccupied_norbs, &
              qmmm_scratch%pdiag_vectmp1,qmmm_scratch%pdiag_vectmp2,qmmm_scratch%pdiag_vectmp3, &
              qmmm_scratch%pdiag_vecjs)
         call wallclock(end_time)
         current_time = current_time + (end_time-start_time)
      end do
      
      if (.not. silence) &
         write(6,'(''| QMMM:                Pseudo diag routine = '',F8.6,'' seconds'')') current_time
      
      !check if this routine is slower than the ones tried so far.
      if (fastest_time < current_time) then
         !The pseudo diag routine is slower than regular diags.
         if (.not. silence) then
            write(6,'(''| QMMM: Pseudo diagonalization appears to be slower than regular'')')
            write(6,'(''| QMMM: diagonalization. Setting pseudo_diag=0 for optimum performance.'')')
         end if
         qmmm_nml%allow_pseudo_diag = .false.
         
         deallocate( qmmm_scratch%pdiag_scr_norbs_norbs, &
              qmmm_scratch%pdiag_scr_noccupied_norbs, &
              qmmm_scratch%pdiag_vectmp1,qmmm_scratch%pdiag_vectmp2, &
              qmmm_scratch%pdiag_vectmp3, qmmm_scratch%pdiag_vecjs, &
              stat=ier )
         REQUIRE(ier == 0)
      end if
      deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
      REQUIRE(ier==0)
   end if


   if (.not. silence) &
      write(6,'(''| QMMM:'')')
#ifdef OPENMP
   !For the time being we just use the number of threads specified by
   !the qmmm_max_omp_threads for openmp - later we will update this
   !to actually test performance for each thread count.
   call omp_set_num_threads(1)
#endif

  !Flush mdout after these timings so people can see something while SQM is
  !running
  if (.not. silence) &
     call flush(6)
   
  end subroutine qm2_time_diag_routines
  

  subroutine qm2_time_diag_routines_random(matsize, fock_matrix)

    implicit none

    integer, intent(in) :: matsize
    _REAL_, intent(out) :: fock_matrix(matsize)

    integer :: i
    _REAL_ :: random_nmbr

    ! Build a fock matrix to diagonalize
    ! Since we don't have enough to build a real fock matrix here what we will
    ! do is just fill the matrix elements with a random distribution.
    ! the matrix is normally dense with elements ranging from -100 to +100
    ! based on the 1NLN test case.
    do i = 1, matsize
!      call amrand(random_nmbr)
! random generator is not necessarily available, and since amrset wasn't called
! beforehand this would only (maybe) work in the libsqm form. Since this is just
! for timing, "fake" the random numbers. JMS
       random_nmbr = i * 0.89101 / matsize
       random_nmbr = 100.0d0 - ( 200.0d0 * random_nmbr )
       fock_matrix(i) = random_nmbr
    end do

  end subroutine qm2_time_diag_routines_random
      

end module qm2_diagonalizer_module
