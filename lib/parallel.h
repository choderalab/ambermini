#ifndef MPI_MAX_PROCESSORS
#define MPI_MAX_PROCESSORS 512
#endif
c
c     Common block necessary for the MPI implementation
c
      integer numtasks,mytaskid,iparpt,iparpt3,rcvcnt,rcvcnt3
      logical mpi_orig
      integer ierr,notdone
      dimension iparpt(0:MPI_MAX_PROCESSORS)
      dimension iparpt3(0:MPI_MAX_PROCESSORS)
      dimension rcvcnt(0:MPI_MAX_PROCESSORS)
      dimension rcvcnt3(0:MPI_MAX_PROCESSORS)
c
      common/parallel/numtasks,mytaskid,notdone,
     .     iparpt,iparpt3,rcvcnt,rcvcnt3,mpi_orig

