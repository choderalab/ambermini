      subroutine dsarpack(n_dim,n_eig_in,n_eig_out,ncv_in,itr_in,
     &                    eigval_tol,eigvals,eigvecs,spectrum,
     &                    need_eigvecs,ierr,debug_arpack,
     &                    v,workl,workd,d,resid,ax,select,
     &                    xyz,grad,return_flag,label)
c
      implicit none
c
c     %-----------------%
c     | Dummy Arguments |
c     %-----------------%
c
      integer n_dim,n_eig_in,n_eig_out,ncv_in,itr_in,spectrum,
     &        need_eigvecs,ierr,debug_arpack,return_flag,label
      Double precision eigval_tol
      Double precision eigvals(n_eig_in),eigvecs(n_dim * n_eig_in)
      Double precision v(n_dim,ncv_in),
     &                 workl(ncv_in*(ncv_in+8)),workd(3*n_dim),
     &                 d(ncv_in,2),resid(n_dim),ax(n_dim),
     &                 xyz(n_dim),grad(n_dim)
      logical select(ncv_in)
c
      save
c
c     %---------------%
c     | Include Files |
c     %---------------%
c
      include 'debug.h'
c
c     This code shows how to use ARPACK to find a few eigenvalues 
c     (lambda) and corresponding eigenvectors (x) for the standard 
c     eigenvalue problem:
c          
c                        A*x = lambda*x
c 
c     where A is an n by n real symmetric matrix.
c
c     The main points illustrated here are 
c
c        1) How to declare sufficient memory to find NEV 
c           eigenvalues of largest magnitude.  Other options
c           are available.
c
c        2) Illustration of the reverse communication interface 
c           needed to utilize the top level ARPACK routine DSAUPD 
c           that computes the quantities needed to construct
c           the desired eigenvalues and eigenvectors(if requested).
c
c        3) How to extract the desired eigenvalues and eigenvectors
c           using the ARPACK routine DSEUPD.
c
c     The only thing that must be supplied in order to use this
c     routine on your problem is to change the array dimensions 
c     appropriately, to specify WHICH eigenvalues you want to compute 
c     and to supply a matrix-vector product
c
c                         w <-  Av
c
c     in place of the call to AV( ) below.
c
c     Once usage of this routine is understood, you may wish to explore
c     the other available options to improve convergence, to solve generalized
c     problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory.
c     This codes implements  
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is derived from the central difference discretization
c         of the 2-dimensional Laplacian on the unit square with
c         zero Dirichlet boundary condition.
c     ... OP = A  and  B = I.
c     ... Assume "call av (n,x,y)" computes y = A*x
c     ... Use mode 1 of DSAUPD.
c
c\BeginLib
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: %Z%
c FILE: %M%   SID: %I%   DATE OF SID: %G%   RELEASE: %R%
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
c     %-------------------------------------------------------%
c     | Storage Declarations:                                 |
c     |                                                       |
c     | The maximum dimensions for all arrays are             |
c     | set here to accommodate a problem size of             |
c     | N .le. MAXN                                           |
c     |                                                       |
c     | NEV is the number of eigenvalues requested.           |
c     |     See specifications for ARPACK usage below.        |
c     |                                                       |
c     | NCV is the largest number of basis vectors that will  |
c     |     be used in the Implicitly Restarted Arnoldi       |
c     |     Process.  Work per major iteration is             |
c     |     proportional to N*NCV*NCV.                        |
c     |                                                       |
c     | You must set:                                         |
c     |                                                       |
c     | MAXN:   Maximum dimension of the A allowed. (dynamic) |
c     | MAXNEV: Maximum NEV allowed. (dynamic)                |
c     | MAXNCV: Maximum NCV allowed. (dynamic)                |
c     %-------------------------------------------------------%
c
C     %--------------------------------------%
C     | F90 Allocatable Arrays (on the heap) |
C     %--------------------------------------%
c
C     Double precision,allocatable,save :: v(:,:)
C     integer,save :: v_row_allocated = 0, v_col_allocated = 0
c
c     %----------------------------------------------%
c     | Originally, as F77 parameters, the following |
c     | integers were used to dimension work arrays. |
c     | They are replaced by dummy arguments used to |
c     | dimension the work arrays as F90 automatic   |
c     | arrays, but the integers are still used for  |
c     | passing the dimensions to lower level ARPACK |
c     | routines dsaupd, dseupd and dmout.           |
c     %----------------------------------------------%
c
      integer          maxn, maxnev, maxncv, ldv
c
c     %-------------------------------------------%
c     | Local F90 Automatic Arrays (on the stack) |
c     %-------------------------------------------%
c
      Double precision
C    &                 workl(ncv_in*(ncv_in+8)),
C    &                 workd(3*n_dim), d(ncv_in,2), resid(n_dim),
C    &                 ax(n_dim),
     &                 cg_dstat(4)
C     logical          select(ncv_in)
      integer          iparam(11), ipntr(11),
     &                 cg_istat(4)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info,
     &                 i, j, nx, ishfts, maxitr, mode1, nconv
      integer          L12, L18, ARPACK_ERROR, status_flag
      data             L12, L18, ARPACK_ERROR /1, 2, -2/
C     integer          v_row_needed, v_col_needed
      logical          rvec
      Double precision      
     &                 tol, sigma
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                 zero
      parameter        (zero = 0.0D+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision           
     &                 dnrm2
      external         dnrm2, daxpy, hessvec
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic        abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if ( label.eq.0 ) go to 1
      go to (12,18) label
  1   continue
c
c     %------------------------------------------------%
c     | Values used to calculate work array dimensions |
c     %------------------------------------------------%
c
      maxn = n_dim
      maxnev = n_eig_in
      maxncv = ncv_in
      ldv = maxn
c
c     %---------------------------------------------------%
c     | The include debug.h statement above and           |
c     | assignments here initiate trace output from the   |
c     | internal actions of ARPACK.  See debug.doc in the |
c     | DOCUMENTS directory for usage.  Initially, the    |
c     | most useful information will be a breakdown of    |
c     | time spent in the various stages of computation   |
c     | given by setting msaupd = 1.                      |
c     %---------------------------------------------------%
c
      ndigit = -5
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      if ( debug_arpack.eq.1 ) then
        msaupd = 1
      else
        msaupd = 0
      endif
      msaup2 = 0
      mseigt = 0
      mseupd = 0
c     
c   *** Allocatable array v will be allowed to grow to its largest size;
c   ***  it is never deallocated:
C     v_row_needed = n_dim        !!! ldv
C     v_col_needed = ncv_in       !!! maxncv
C     if( allocated(v) )then
C       if( (v_row_needed .gt. v_row_allocated)
C    & .or. (v_col_needed .gt. v_col_allocated) )then
C         deallocate(v,stat=ierr)
C         if( ierr .ne. 0 )then
C           write( logfil, '(a,i16,1x,i8)' )
C    &       'ARPACK: could not deallocate v'
C           go to 9000
C         endif
C       endif
C     endif
C     if( .not. allocated(v) )then
C       allocate( v(v_row_needed,v_col_needed), stat=ierr )
C       if( ierr .ne. 0 )then
C         write( logfil, '(a,2i10)' )
C    &     'ARPACK: could not allocate v'
C         go to 9000
C       endif
C       v_row_allocated = v_row_needed
C       v_col_allocated = v_col_needed
C     endif
C     v = zero !!! zero out entire v array
c     
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
      n = n_dim
c
c     %----------------------------------------------%
c     |                                              | 
c     | Specifications for ARPACK usage are set      | 
c     | below:                                       |
c     |                                              |
c     |    1) NEV = N_EIG_IN  asks for N_EIG_IN      |  
c     |       eigenvalues to be computed.            | 
c     |                                              |
c     |    2) NCV = NCV_IN sets the length of the    |
c     |       Arnoldi factorization                  |
c     |                                              |
c     |    3) This is a standard problem             |
c     |         (indicated by bmat  = 'I')           |
c     |                                              |
c     |    4) Ask for the NEV eigenvalues of         |
c     |       smallest magnitude                     |
c     |         (indicated by which = 'SM')          |
c     |       See documentation in SSAUPD for the    |
c     |       other options SA, LA, LM, BE.          | 
c     |                                              |
c     | Note: NEV and NCV must satisfy the following |
c     | conditions:                                  |
c     |              NEV <= MAXNEV                   |
c     |          NEV + 1 <= NCV <= MAXNCV            |
c     %----------------------------------------------%
c
      nev   = n_eig_in
      ncv   = ncv_in 
      bmat  = 'I'
      if ( spectrum .eq. 1 ) then
         which = 'SM'
      else if ( spectrum .eq. 2 ) then
         which = 'SA'
      else if ( spectrum .eq. 3 ) then
         which = 'LM'
      else if ( spectrum .eq. 4 ) then
         which = 'LA'
      else if ( spectrum .eq. 5 ) then
         which = 'BE'
      else
          print *, ' ERROR with _SSIMP: Spectrum .NE. (SM|SA|LA|LM|BE)'
         go to 9000
      end if
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling DSAUPD                    |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DSAUPD.                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  | 
c     |                                                     |
c     | The work array WORKL is used in DSAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol = eigval_tol 
      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting PARAM(1) = 1).              |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = itr_in 
      mode1 = 1
c
      iparam(1) = ishfts
c                
      iparam(3) = maxitr
c                  
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c
            status_flag = 0
 11         continue
               call hessvec ( n, workd(ipntr(1)), workd(ipntr(2)),
     &                        xyz, grad, return_flag, status_flag )
               if ( status_flag.eq.0 ) go to 13
               if ( status_flag.lt.0 ) go to 9000
               label = L12
               return
 12         go to 11
 13         continue
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
         go to 9000
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine DSEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c           
         if ( need_eigvecs .eq. 1 ) then
            rvec = .true.
         else
            rvec = .false.
         end if
c
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
c
c        %----------------------------------------------%
c        | Eigenvalues are returned in the first column |
c        | of the two dimensional array D and the       |
c        | corresponding eigenvectors are returned in   |
c        | the first NCONV (=IPARAM(5)) columns of the  |
c        | two dimensional array V if requested.        |
c        | Otherwise, an orthogonal basis for the       |
c        | invariant subspace corresponding to the      |
c        | eigenvalues in D is returned in V.           |
c        %----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of DSEUPD. |
c           %------------------------------------%
c
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd. '
            print *, ' '
            go to 9000
c
         else if ( debug_arpack.eq.1 ) then
c
            nconv =  iparam(5)
            n_eig_out = nconv
            if ( nconv .le. 0 ) then
               print *, ' '
               print *, ' ARPACK: Not a single mode converged.'
               print *, ' '
               go to 9000
            endif
c
C           %--------------------------------------------%
C           | "UnDO" DO 20 j=1,nconv loop, because it is |
C           | illegal to jump in and out from a DO loop. |
C           %--------------------------------------------%
c
            j = 1
 16         continue
c
c              %---------------------------%
c              | Compute the residual norm |
c              |                           |
c              |   ||  A*x - lambda*x ||   |
c              |                           |
c              | for the NCONV accurately  |
c              | computed eigenvalues and  |
c              | eigenvectors.  (iparam(5) |
c              | indicates how many are    |
c              | accurate to the requested |
c              | tolerance)                |
c              %---------------------------%
c
               status_flag = 0
 17            continue
                  call hessvec ( n, v(1,j), ax, xyz, grad,
     &                           return_flag, status_flag )
                  if ( status_flag.eq.0 ) go to 19
                  if ( status_flag.lt.0 ) go to 9000
                  label = L18
                  return
 18            go to 17
 19            continue
c
               call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
               d(j,2) = dnrm2(n, ax, 1)
               d(j,2) = d(j,2) / abs(d(j,1))
c
               j = j + 1
               if ( j .gt. nconv ) go to 20
c
               go to 16
c
 20         continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call dmout(6, nconv, 2, d, maxncv, -6,
     &           'Ritz values and relative residuals')
c
c           %-------------------------------------------%
c           | Print additional convergence information. |
c           %-------------------------------------------%
c
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit',
     &                  ' Arnoldi update, try increasing NCV.'
               print *, ' '
            end if
c
            print *, ' '
            print *, ' _SSIMP '
            print *, ' ====== '
            print *, ' '
            print *, ' Size of the matrix is ', n
            print *, ' The number of Ritz values requested is ', nev
            print *, ' The number of Arnoldi vectors generated',
     &               ' (NCV) is ', ncv
            print *, ' What portion of the spectrum: ', which
            print *, ' The number of converged Ritz values is ',
     &                 nconv
            print *, ' The number of Implicit Arnoldi update',
     &               ' iterations taken is ', iparam(3)
            print *, ' The number of OP*x is ', iparam(9)
            print *, ' The convergence criterion is ', tol
            print *, ' '
         end if
c
c        %----------------------------%
c        | Return eigvals and eigvecs |
c        %----------------------------%
c
         nconv =  iparam(5)
         n_eig_out = nconv
         if ( nconv .le. 0 ) then
            print *, ' '
            print *, ' ARPACK: Not a single mode converged.'
            print *, ' '
            go to 9000
         endif
c
         do 40 j=1, nconv
             eigvals(j) = d(j,1)
c
             do 30 i=1, n
                eigvecs((j-1)*n+i) = v(i,j)
 30          continue
 40      continue
c
      end if
c
c     %--------------------------------%
c     | Done with subroutine dsarpack. |
c     %--------------------------------%
c
      label = 0
      return
c
 9000 continue !!! Error
c
      if( status_flag.eq.0 ) status_flag = ARPACK_ERROR
c
      label = status_flag
      return
c
      end
c 
c ------------------------------------------------------------------
