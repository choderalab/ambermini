! ============= NETCDF ROUTINES FOR SANDER/PMEMD =====================
!> @file AmberNetcdf.F90 
!> @author Daniel R. Roe, 2013
!> @define BINTRAJ Use NetCDF.
!> @define _NETCDF_DEBUG Activate extra debugging info.
module AmberNetcdf_mod
  implicit none
#ifdef BINTRAJ
  private
  integer, save       :: mdout
  character(10), save :: programName
  character(10), save :: programVersion
! DEFINES
  character(5), parameter  :: NCFRAME = "frame"
  character(7), parameter  :: NCSPATIAL = "spatial"
  character(4), parameter  :: NCATOM = "atom"
  character(12), parameter :: NCCELL_SPATIAL = "cell_spatial"
  character(12), parameter :: NCCELL_LENGTHS = "cell_lengths"
  character(12), parameter :: NCCELL_ANGULAR = "cell_angular"
  character(11), parameter :: NCCELL_ANGLES = "cell_angles"
  character(11), parameter :: NCCOORDS = "coordinates"
  character(10), parameter :: NCVELO = "velocities"
  character(6), parameter  :: NCFRC = "forces"
  character(5), parameter  :: NCTEMPERATURE = "temp0"
  character(4), parameter  :: NCTIME = "time"
  character(5), parameter  :: NCLABEL = "label"
  integer, parameter       :: NCLABELLEN = 5
  character(14), parameter :: NCREMD_DIMENSION = "remd_dimension"
  character(12), parameter :: NCREMD_DIMTYPE = "remd_dimtype"
  character(12), parameter :: NCREMD_INDICES = "remd_indices"
  character(11), parameter :: NCREMD_GROUPS = "remd_groups"
  character(5), parameter  :: NCEPTOT = "eptot"
  character(4), parameter  :: NCBINS = "bins"
  double precision, parameter :: velocityScale = 20.455d0

  logical, public, save :: verbose_netcdf = .true.

  public NC_create, NC_setupAmberNetcdf, NetcdfFileExists, NC_error, &
         NC_openRead, NC_openWrite, &
         NC_setupMdcrd, NC_defineRemdIndices, NC_close, &
         checkNCerror, NC_readRestartBox, NC_setupRestart, &
         NC_readRestartIndices, NC_checkRestart, NC_setupCoordsVelo, &
         NC_checkTraj, NC_setupReservoir, NC_readReservoir

contains
!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_ERROR()
!> @brief check for netcdf error.
!> @param err netcdf error code.
!> @param location purpose/location of the call.
!> @return true if error, false if not.
logical function NC_error(err, location)
  use netcdf
  implicit none
  integer, intent(in)                :: err 
  character(*), optional, intent(in) :: location

  NC_error=.false.
  if (err .ne. nf90_noerr) then
    if (verbose_netcdf) then
      write(mdout, '(a,a)') 'NetCDF error: ', trim(nf90_strerror(err))
      if (present(location)) then
        write(mdout, '(a,a)') '  at ', location
      end if
    end if
    NC_error=.true.
  endif
end function NC_error

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION CHECKNCERROR()
!> @brief Passive check for netcdf error.
subroutine checkNCerror(err, location)
  use netcdf
  implicit none
  integer, intent(in)                :: err
  character(*), optional, intent(in) :: location
  if (err .ne. nf90_noerr .and. verbose_netcdf) then
    write(mdout, '(a,a)') 'NetCDF error: ', trim(nf90_strerror(err))
    if (present(location)) then
      write(mdout, '(a,a)') '  at ', location
    end if
  end if
end subroutine checkNCerror

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NETCDFFILEEXISTS()
!> @brief check if file exists and is valid netcdf.
!> @param fname File name to check.
!> @return true if file exists and is netcdf, false if does not exist
!!         or is not netcdf.
logical function NetcdfFileExists(fname)
  use netcdf
  implicit none
  character(*), intent(in) :: fname
  integer :: ncid, ierr
  NetcdfFileExists=.false.
  if (nf90_open( fname, NF90_NOWRITE, ncid ) .ne. NF90_NOERR) return
  ierr = nf90_close( ncid )
  NetcdfFileExists=.true.
end function NetcdfFileExists

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_OPENREAD()
!> @brief Open netcdf file for reading.
!> @param fname Netcdf file to open
!> @param ncid Netcdf ID of opened file.
!> @return false if successful, true if error occurs.
logical function NC_openRead(fname, ncid)
  use netcdf
  implicit none
  character(*), intent(in) :: fname
  integer, intent(out)     :: ncid
  NC_openRead=.true.
  if (NC_error( nf90_open( fname, NF90_NOWRITE, ncid ))) return
  NC_openRead=.false.
end function NC_openRead

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_OPENWRITE()
!> @brief Open netcdf file for writing.
!> @param fname Netcdf file to open
!> @param ncid Netcdf ID of opened file.
!> @return false if successful, true if error occurs.
logical function NC_openWrite(fname, ncid)
  use netcdf
  implicit none
  character(*), intent(in) :: fname
  integer, intent(out)     :: ncid
  NC_openWrite=.true.
  if (NC_error( nf90_open( fname, NF90_WRITE, ncid ))) return
  NC_openWrite=.false.
end function NC_openWrite

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_CLOSE()
!> @brief Close netcdf file.
subroutine NC_close(ncid)
  use netcdf
  implicit none
  integer, intent(in) :: ncid
  call checkNCerror(nf90_close( ncid ), 'Closing NetCDF file')
end subroutine NC_close

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_CHECKCONVENTIONS
!> @param ncid Netcdf ID of previously opened file.
!> @param checkForRestart If true check restart, otherwise check traj.
!> @return false if ncid is a netcdf traj/restart, true if not.
logical function NC_checkConventions(ncid, checkForRestart)
  use netcdf
  implicit none
  integer, intent(in) :: ncid
  logical, intent(in) :: checkForRestart
  character(80) :: attribute
  NC_checkConventions=.true.
  if (NC_error( nf90_get_att(ncid, NF90_GLOBAL, 'Conventions', attribute),&
                  'Getting Conventions attribute')) return
  if (checkForRestart) then
    if (attribute.ne."AMBERRESTART") then
      write(mdout,'(a)') 'ERROR: NetCDF restart has Conventions that are not AMBERRESTART.'
      return
    endif 
  else ! checkForTraj
    if (attribute.ne."AMBER") then
      write(mdout,'(a)') 'ERROR: NetCDF traj has Conventions that are not AMBER.'
      return
    endif
  endif
  ! Check conventions version.
  ! NOTE: This may need to be separated if traj/restart have diff versions
  if (NC_error( nf90_get_att(ncid, NF90_GLOBAL, 'ConventionVersion', attribute),&
                  'Getting Convention version')) return
  if (attribute .ne. "1.0") &
    write(mdout, '(a)') 'Warning: Netcdf file  has ConventionVersion that is not 1.0:'
  NC_checkConventions=.false.
end function NC_checkConventions

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_GETTITLE
!> @brief Get netcdf file title.
logical function NC_GetTitle(ncid, title)
  use netcdf
  implicit none
  integer, intent(in)       :: ncid
  character(*), intent(out) :: title
  NC_GetTitle=.true.
  if (NC_error( nf90_get_att(ncid, NF90_GLOBAL, 'title', title),&
                  'Getting title')) return
  NC_GetTitle=.false.
end function NC_GetTitle

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION GETDIMINFO
!> @return dimension ID of specified dimension.
integer function GetDimInfo(ncid, attribute, length)
  use netcdf
  implicit none
  integer, intent(in)      :: ncid
  character(*), intent(in) :: attribute
  integer, intent(out)     :: length
  integer :: dimID
  GetDimInfo=-1
  length=0
  ! Get dimID
  if ( NC_error(nf90_inq_dimid(ncid, attribute, dimID),&
                  'Getting dimID for attribute') ) return
  ! Get dimension length
  if ( NC_error(nf90_inquire_dimension(ncid, dimID, len=length),&
                  'Getting dimension length') ) return
  GetDimInfo = dimID
end function GetDimInfo

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION GETVARINFO
!> @return variable ID of specified variable if it exists.
!> @return -1 if the variable does not exist.
integer function GetVarInfo(ncid, attribute)
  use netcdf
  implicit none
  integer, intent(in)      :: ncid
  character(*), intent(in) :: attribute
  integer VID
  GetVarInfo=-1
  if ( nf90_inq_varid(ncid, attribute, VID) .eq. NF90_NOERR ) &
    GetVarInfo=VID
end function GetVarInfo

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION CHECKATTRTEXT
!> @brief Check that text associated with attribute of vid matches textIn
logical function CheckAttrText(ncid, vid, attribute, textIn)
  use netcdf
  implicit none
  integer, intent(in)      :: ncid, vid
  character(*), intent(in) :: attribute, textIn
  ! Local vars
  character(80) :: text
  CheckAttrText=.true.
  if (NC_error( nf90_get_att(ncid, vid, attribute, text),&
                  attribute)) return
  if ( text .ne. textIn ) &
    write(mdout,'(6a)') 'Warning: In netcdf file, expected ', textIn, ' for attribute ', attribute, ', got ', trim(text)
  CheckAttrText=.false.
end function CheckAttrText

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPFRAME
logical function NC_setupFrame(ncid, nframes)
  implicit none
  integer, intent(in)  :: ncid
  integer, intent(out) :: nframes
  integer frameDID
  NC_setupFrame=.true.
  frameDID = GetDimInfo( ncid, NCFRAME, nframes )
  if (frameDID.eq.-1) return
  NC_setupFrame=.false.
end function NC_setupFrame

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPCOORDSVELO
!> @brief Get # of atoms, and coordinate and velocity VIDs.
!> @detail Setup natoms, coordVID, and velocityVID. Check units 
!!         and spatial dimensions.
logical function NC_setupCoordsVelo(ncid, natoms, coordVID, velocityVID)
  use netcdf
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: natoms, coordVID, velocityVID
  ! Local vars
  integer :: spatialDID, atomDID, spatial, spatialVID

  NC_setupCoordsVelo=.true.
  atomDID = GetDimInfo( ncid, NCATOM, natoms )
  if (atomDID .eq. -1) return
  ! Get coord info
  coordVID = GetVarInfo( ncid, NCCOORDS )
  if ( coordVID .ne. -1 ) then
#   ifdef _NETCDF_DEBUG
    write(mdout,'(a)') 'Netcdf file has coordinates.'
#   endif
    if (CheckAttrText(ncid, coordVID, "units", "angstrom")) return
  endif 
  ! Get spatial info
  spatialDID = GetDimInfo( ncid, NCSPATIAL, spatial )
  if (spatialDID .eq. -1) return
  if (spatial .ne. 3) then
    write(mdout,'(a,i6)') 'Error: Netcdf: Expected 3 spatial dimensions, got ', spatial
    return
  endif
  if ( NC_error(nf90_inq_varid(ncid, NCSPATIAL, spatialVID),&
                  'Getting spatial VID') ) return
  ! Get velocity info
  ! NOTE: Check scale factor?
  velocityVID = GetVarInfo( ncid, NCVELO )
  if ( velocityVID .ne. -1 ) then
#   ifdef _NETCDF_DEBUG
    write(mdout,'(a)') 'Netcdf file has velocities.'
#   endif
    if (CheckAttrText(ncid, velocityVID, "units", "angstrom/picosecond")) return
  endif
  NC_setupCoordsVelo=.false.
end function NC_setupCoordsVelo

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPTIME
!> @brief Setup timeVID and check units.
logical function NC_setupTime(ncid, timeVID)
  use netcdf
  implicit none
  integer, intent(in)  :: ncid
  integer, intent(out) :: timeVID
  NC_setupTime=.true.
  if (nf90_inq_varid(ncid, NCTIME, timeVID) .ne. NF90_NOERR) return
  if (CheckAttrText(ncid, timeVID, "units", "picosecond")) return
  NC_setupTime=.false.
end function NC_setupTime

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF SUBROUTINE NC_SETUPBOX
!> @brief Check box info, setup cellLengthVID and cellAngleVID.
!> @param cellLengthVID Cell lengths variable, set to -1 if no box or error.
!> @param cellAngleVID Cell angles variable, set to -1 if no box or error.
subroutine NC_setupBox(ncid, cellLengthVID, cellAngleVID)
  use netcdf
  implicit none
  integer, intent(in)  :: ncid
  integer, intent(out) :: cellLengthVID, cellAngleVID

  cellAngleVID = -1
  cellLengthVID = GetVarInfo( ncid, NCCELL_LENGTHS )
  if ( cellLengthVID.ne.-1 ) then
    if (NC_error(nf90_inq_varid(ncid,NCCELL_ANGLES,cellAngleVID),&
                 'Getting cell angles')) then
      cellLengthVID = -1
      return
    endif
#   ifdef _NETCDF_DEBUG
    write(mdout,'(a)') 'Netcdf Box information found.'
#   endif
  endif 
end subroutine NC_setupBox

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPTEMPERATURE
!> @brief Look for replica temperature info.
integer function NC_setupTemperature(ncid)
  use netcdf
  implicit none
  integer, intent(in) :: ncid
  integer TempVID

  NC_setupTemperature=-1
  if ( nf90_inq_varid(ncid,NCTEMPERATURE,TempVID) .eq. NF90_NOERR ) then
#   ifdef _NETCDF_DEBUG
    write(mdout,'(a)') 'Netcdf file has replica temperatures.'
#   endif
    NC_setupTemperature = TempVID
  endif
end function NC_setupTemperature

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF SUBROUTINE NC_SETUPAMBERNETCDF()
!> @brief Set output unit, program name, and program version.
!> @detail This should be called once by the calling program prior to
!!         any other calls to this module.
!> @param mdout_unit File unit number for mdout output.
!> @param progname Calling program name.
!> @param progver Calling program version.
subroutine NC_setupAmberNetcdf(mdout_unit, progname, progver)
  implicit none
  ! Passed vars
  integer, intent(in)      :: mdout_unit
  character(*), intent(in) :: progname, progver
  mdout = mdout_unit
  programName = progname
  programVersion = progver
end subroutine NC_setupAmberNetcdf

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_DEFINE_VAR
!> @brief Wrapper around nf90_def_var for multiple dimensions.
!> @param vcid NetCDF ID of file to define variable for.
!> @param varname Name of variable to define.
!> @param dtype Data type of variable.
!> @param ndim Number of dimensions in variable.
!> @param dimID Array holding size of each dimension.
integer function NC_define_var(ncid, varname, dtype, ndim, dimID)
  use netcdf
  implicit none
  ! Passed vars
  integer, intent(in)      :: ncid
  character(*), intent(in) :: varname
  integer, intent(in)      :: dtype
  integer, intent(in)      :: ndim
  integer, dimension(NF90_MAX_DIMS), intent(in) :: dimID

  if (ndim .eq. 0) then
    if (NC_error( nf90_def_var(ncid, varname, dtype, NC_define_var) )) &
      NC_define_var=-1
  else if (ndim .eq. 1) then
    if (NC_error( nf90_def_var(ncid, varname, dtype, (/ dimID(1) /), &
                                 NC_define_var) )) NC_define_var=-1
  else if (ndim .eq. 2) then
    if (NC_error( nf90_def_var(ncid, varname, dtype, (/ dimID(1), dimID(2) /), &
                                 NC_define_var) )) NC_define_var=-1
  else ! ndim .eq. 3
    if (NC_error( nf90_def_var(ncid, varname, dtype, (/ dimID(1), dimID(2), dimID(3) /),&
                                 NC_define_var) )) NC_define_var=-1
  endif
end function NC_define_var

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_CREATE()
!> @brief Create a NetCDF traj/restart.
!> @param filename Name of file to create.
!> @param owrite Overwrite character.
!> @param isRestart If true, netcdf restart.
!> @param natomIn Number of atoms.
!> @param hasCoords If true, file will have coordinates.
!> @param hasVelocity If true, file will have velocity info.
!> @param hasBox If true, file will have box information.
!> @param hasTemperature If true, file will have temperature information
!> @param hasTime If true, file will have time information.
!> @param title Traj/restart title.
!> @param ncid Resulting netcdf ID of file.
!> @param timeVID Resulting VID of time.
!> @param coordVID Resulting VID of coordinates.
!> @param velocityVID Resulting VID of velocities.
!> @param cellLengthVID Resulting VID of box lengths.
!> @param cellAngleVID Resulting VID of box angles.
!> @param TempVID Resulting VID of temperature(s).
!> @return true on error, false on success.
logical function NC_create(filename, owrite, isRestart, natomIn, hasCoords, &
                           hasVelocity, hasBox, hasTemperature, hasTime, &
                           hasFrc, title, &
                           ncid, timeVID, coordVID, velocityVID, frcVID, &
                           cellLengthVID, cellAngleVID, TempVID,&
                           frameDID)
  use netcdf
  implicit none

  ! Passed variables
  character(*), intent(in)  :: filename
  character,    intent(in)  :: owrite
!  character,    intent(in)  :: facc
  logical,      intent(in)  :: isRestart
  integer,      intent(in)  :: natomIn
  logical,      intent(in)  :: hasCoords, hasVelocity, hasBox, hasTemperature
  logical,      intent(in)  :: hasFrc
  logical,      intent(in)  :: hasTime
  character(*), intent(in)  :: title
  integer,      intent(out) :: ncid, timeVID, coordVID, velocityVID
  integer,      intent(out) :: frcVID
  integer,      intent(out) :: cellLengthVID, cellAngleVID, TempVID
  integer, optional, intent(out) :: frameDID

  ! Local variables
  integer, dimension(NF90_MAX_DIMS) :: dimensionID
  integer  :: NDIM, dataType, cmode, ierr
  integer  :: atomDID, spatialDID, spatialVID
  integer  :: cell_spatialDID, labelDID, cell_angularDID
  integer  :: cell_spatialVID, cell_angularVID

  NC_create=.true. ! Assume error until success
  if ( len(filename) .eq. 0 ) return
# ifdef _NETCDF_DEBUG
  write(mdout,'(a,a,a,i8,4(a,l1))') "| DEBUG: NC_create: ", trim(filename), &
        "  natom=", natomIn, "  V=", hasVelocity, "  box=", hasBox, &
        "  temp=", hasTemperature, "  time=", hasTime
# endif
  cmode = nf90_64bit_offset
  if (owrite == 'N') cmode = ior(cmode, nf90_noclobber)
  !if (owrite == 'U' .and. facc == 'A') cmode = ior(cmode, nf90_noclobber)
  ierr = nf90_create(path=filename, cmode=cmode, ncid=ncid)
  if (ierr == nf90_eexist) then
     write(mdout,'(a,a)') 'Error: NetCDF file exists and -O not specified: ',&
           filename
     return
  endif
  if (NC_error( ierr )) return
 
  ! Set number of dimensions based on file type
  if (isRestart) then ! Restart
    NDIM = 2
    dataType = NF90_DOUBLE
  else                ! Trajectory
    NDIM = 3
    dataType = NF90_FLOAT
  endif

  !  Frame dimension for traj
  if (.not.isRestart) then
    if (.not.present(frameDID)) then
      write(mdout, '(a)') 'Internal Error: NC_create called without frameDID.'
      return
    endif
    if (NC_error( nf90_def_dim( ncid, NCFRAME, NF90_UNLIMITED, frameDID ),&
                    'Defining frame dimension.' )) return
    dimensionID(1) = frameDID
  endif
  ! Time variable and units
  if (hasTime) then
    timeVID = NC_define_var(ncid, NCTIME, dataType, NDIM-2, dimensionID)
    if ( timeVID .eq. -1) return
    if (NC_error( nf90_put_att( ncid, timeVID, "units", "picosecond" ),&
                    'Writing time VID units.' )) return
  endif
  ! Spatial dimension and variable
  if ( NC_error( nf90_def_dim( ncid, NCSPATIAL, 3, spatialDID),&
                   'Defining spatial dimension')) return
  dimensionID(1) = spatialDID;
  spatialVID = NC_define_var(ncid, NCSPATIAL, NF90_CHAR, 1, dimensionID)
  if (spatialVID .eq. -1) return ! 'Defining spatial variable'
  ! Atom dimension
  if (NC_error( nf90_def_dim( ncid, NCATOM, natomIn, atomDID ),&
                  'Defining atom dimension.' )) return
  ! Setup dimensions for Coords/Velocity
  ! NOTE: MUST BE MODIFIED IF NEW TYPES ADDED
  if (isRestart) then ! Restart
    dimensionID(1) = spatialDID
    dimensionID(2) = atomDID
  else                ! Trajectory
    dimensionID(1) = spatialDID;
    dimensionID(2) = atomDID;
    dimensionID(3) = frameDID;
  endif
  ! Coord variable
  if (hasCoords) then
    coordVID = NC_define_var(ncid, NCCOORDS, dataType, NDIM, dimensionID)
    if (coordVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, coordVID, "units", "angstrom"),&
                    'Writing coordinates variable units.')) return
  endif
  ! Velocity variable
  if (hasVelocity) then
    velocityVID = NC_define_var(ncid, NCVELO, dataType, NDIM, dimensionID)
    if (velocityVID .eq. -1) return ! 'Defining velocities variable.')) return
    if (NC_error( nf90_put_att( ncid, velocityVID, "units", "angstrom/picosecond"),&
                    'Writing velocities variable units.' )) return
    if (NC_error( nf90_put_att( ncid, velocityVID, "scale_factor", velocityScale),&
                    'Writing velocities scale factor.')) return
  endif
  if (hasFrc) then
    frcVID = NC_define_var(ncid, NCFRC, dataType, NDIM, dimensionID)
    if (frcVID .eq. -1) return ! Failed defining force variable
    if (NC_error( nf90_put_att( ncid, frcVID, "units", "kilocalorie/mole/angstrom"), &
                     'Writing force variable units.')) return
  end if
  ! Box Info
  if (hasBox) then
    ! Cell Spatial
    if (NC_error( nf90_def_dim( ncid, NCCELL_SPATIAL, 3, cell_spatialDID),&
                    'Defining cell spatial dimension.' )) return
    dimensionID(1) = cell_spatialDID
    cell_spatialVID = NC_define_var(ncid, NCCELL_SPATIAL, NF90_CHAR, 1, dimensionID)
    if (cell_spatialVID .eq. -1) return ! 'Defining cell spatial variable.'
    ! Cell angular
    if (NC_error( nf90_def_dim( ncid, NCLABEL, NCLABELLEN, labelDID),&
                    'Defining label dimension.' )) return
    if (NC_error( nf90_def_dim( ncid, NCCELL_ANGULAR, 3, cell_angularDID),&
                    'Defining cell angular dimension' )) return
    dimensionID(1) = labelDID;
    dimensionID(2) = cell_angularDID;
    cell_angularVID = NC_define_var(ncid, NCCELL_ANGULAR, NF90_CHAR, 2, dimensionID)
    if (cell_angularVID .eq. -1) return ! 'Defining cell angular variable.'
    ! Setup dimensions for Box
    ! NOTE: This must be modified if more types added
    dimensionID(1) = cell_spatialDID
    if (.not.isRestart) dimensionID(2) = frameDID
    cellLengthVID = NC_define_var(ncid, NCCELL_LENGTHS, NF90_DOUBLE, NDIM-1,&
                                  dimensionID)
    if (cellLengthVID .eq. -1) return !  'Defining cell length variable.'
    if (NC_error( nf90_put_att(ncid, cellLengthVID, "units", "angstrom"),&
                    'Writing cell length variable units' )) return
    dimensionID(1) = cell_angularDID;
    cellAngleVID = NC_define_var(ncid, NCCELL_ANGLES, NF90_DOUBLE, NDIM-1,&
                                 dimensionID)
    if (cellAngleVID .eq. -1) return ! 'Defining cell angle variable.'
    if (NC_error( nf90_put_att(ncid, cellAngleVID, "units", "degree"),&
                    'Writing cell angle variable units.' )) return
  endif
  ! Attributes
  if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"title",title),&
                 'Writing title')) return
  if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"application","AMBER"),&
                 'Writing application')) return
  if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"program",trim(programName)),&
                 'Writing program')) return
  if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"programVersion",&
                              trim(programVersion)),&
                 'Writing program version')) return
  ! Conventions
  if (isRestart) then
    if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","AMBERRESTART"),&
                   'Writing restart conventions')) return
  else
    if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","AMBER"),&
                   'Writing trajectory conventions')) return
  endif
  if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"ConventionVersion","1.0"),&
                 'Writing conventions version')) return
  ! Replica Temperature
  if (hasTemperature) then
    if (present(frameDID)) dimensionID(1) = frameDID;
    TempVID = NC_define_var(ncid,NCTEMPERATURE,NF90_DOUBLE,NDIM-2,dimensionID)
    if (TempVID .eq. -1) return !  'NetCDF error on defining temperature'
    if (NC_error(nf90_put_att(ncid,TempVID,"units","kelvin"),&
                   'defining temperature units')) return
  endif
  ! Set fill mode
  if (NC_error(nf90_set_fill(ncid, NF90_NOFILL, ierr),&
                 'setting fill value.')) return
  ! End netcdf definitions
  if (NC_error(nf90_enddef(ncid), 'ending definitions')) return

  ! Specify spatial dimension labels
  if (NC_error(nf90_put_var(ncid, spatialVID, &
                    (/ 'x', 'y', 'z' /), start = (/ 1 /), count = (/ 3 /)), &
                    'write spatial variable')) return
  if (hasBox) then
    if (NC_error(nf90_put_var(ncid, cell_spatialVID, &
                   (/ 'a','b','c' /), start=(/ 1 /), count=(/ 3 /)), &
                   'write cell spatial variable')) return
    if (NC_error(nf90_put_var(ncid, cell_angularVID, &
                   (/ 'alpha','beta ','gamma' /), &
                   start=(/ 1, 1 /), count=(/ NCLABELLEN, 3 /)), &
                   'write spatial variable')) return
  endif

  ! All is well, exit false.
  NC_create=.false.
end function NC_create

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPMDCRD
!> @brief Setup a previously opened netcdf mdcrd/mdvel file
logical function NC_setupMdcrd(ncid, title, nframes, ncatom, &
                               coordVID, velocityVID, timeVID, &
                               cellLengthVID, cellAngleVID, TempVID)
  use netcdf
  implicit none
  integer, intent(in)       :: ncid
  character(*), intent(out) :: title
  integer, intent(out)      :: nframes
  integer, intent(out)      :: ncatom, coordVID, velocityVID, timeVID, &
                               cellLengthVID, cellAngleVID, TempVID
  NC_setupMdcrd=.true.
  ! Make sure this is a netcdf trajectory
  if (NC_checkConventions( ncid, .false. )) return 
  ! Get title
  if (NC_GetTitle( ncid, title )) return
  ! Get Frame Info
  if (NC_setupFrame( ncid, nframes )) return
  ! Setup Coordinates/Velocities, check atom count
  if (NC_setupCoordsVelo( ncid, ncatom, coordVID, velocityVID)) return
  ! Return a error if no coords
  if ( coordVID .eq. -1 ) then
    write(mdout,'(a)') 'Error: NetCDF file has no coordinates.'
    return
  endif
  ! Setup Time - Allowed to fail, these values are not needed for traj.
  if (NC_setupTime( ncid, timeVID )) then
    write(mdout,'(a)'), 'Warning: NetCDF trajectory has no time values.'
    timeVID = -1
  endif
  ! Setup box: If no box info cellLengthVID/cellAngleVID will be set to
  !            -1. It is the responsibility of the calling routine to
  !            decide if this is OK.
  call NC_setupBox(ncid, cellLengthVID, cellAngleVID)
  ! Setup replica temperature
  TempVID = NC_setupTemperature( ncid )
  ! All is well
  NC_setupMdcrd=.false.
end function NC_setupMdcrd

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPRESERVOIR
!> @brief Setup a previously opened netcdf structure reservoir file.
logical function NC_setupReservoir(ncid, reservoirsize, restemp0, numatoms,&
                                   coordVID, velocityVID, eptotVID, binsVID,&
                                   iseed)
  use netcdf
  implicit none
  integer, intent(in)           :: ncid
  integer, intent(out)          :: reservoirsize
  double precision, intent(out) :: restemp0
  integer, intent(out)          :: numatoms, coordVID, velocityVID, eptotVID
  integer, intent(out)          :: binsVID, iseed
  ! local vars
  character(80) title
  integer ierr, TempVID

  NC_setupReservoir=.true.
  ! Set up trajectory
  if (NC_setupMdcrd(ncid, title, reservoirsize, numatoms, coordVID, &
                    velocityVID, ierr, ierr, ierr, TempVID)) return
  ! Set up EPtot
  if ( NC_error( nf90_inq_varid(ncid, NCEPTOT, eptotVID),&
                   'Getting Netcdf reservoir eptot VID') ) return
  ! Set up Bins. Allowed to fail silently
  binsVID = GetVarInfo( ncid, NCBINS )
  ! Get reservoir temp
  if (TempVID .eq. -1) then
    write(mdout,'(a)') 'Error: netcdf reservoir does not contain temperature.'
    return
  endif
  if (NC_error(nf90_get_var(ncid, TempVID, restemp0),&
               'Getting netcdf reservoir temperature')) return
  ! Get random seed
  if (NC_error(nf90_get_att(ncid, NF90_GLOBAL, "iseed", iseed),&
               'Getting netcdf reservoir random seed')) return
  ! All is well
  NC_setupReservoir=.false.
end function NC_setupReservoir

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_READRESERVOIR
!> @brief Read energies and cluster numbers (if present) from netcdf reservoir
logical function NC_readReservoir(ncid, reservoirsize, eptotVID, binsVID, &
                                  saveene, clusternum)
  use netcdf
  implicit none
  integer, intent(in) :: ncid, reservoirsize, eptotVID, binsVID
  double precision, dimension(:), intent(out) :: saveene
  integer, dimension(:), intent(out) :: clusternum
  
  NC_readReservoir=.true.
  if (NC_error(nf90_get_var(ncid, eptotVID, saveene, &
                            start = (/ 1 /),&
                            count = (/ reservoirsize /)),&
               'Getting reservoir energies')) return
  if (binsVID.ne.-1) then
    if (NC_error(nf90_get_var(ncid, binsVID, clusternum, &
                              start = (/ 1 /), &
                              count = (/ reservoirsize /)),&
                 'Getting reservoir cluster numbers')) return
  endif
  ! All is well
  NC_readReservoir=.false.
end function NC_readReservoir

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_DEFINEREMDINDICES
!> @brief Create M-REMD index info in traj/restart file.
logical function NC_defineRemdIndices(ncid, remd_dimension, indicesVID, &
                                     remd_types, isRestart, frameDID,&
                                     remd_groupsVID)
  use netcdf
  implicit none
  integer, intent(in)               :: ncid, remd_dimension
  integer, intent(out)              :: indicesVID
  integer, dimension(:), intent(in) :: remd_types
  logical, intent(in)               :: isRestart
  integer, optional, intent(in)     :: frameDID
  integer, optional, intent(out)    :: remd_groupsVID
  ! Local vars
  integer, dimension(NF90_MAX_DIMS) :: dimensionID
  integer :: remd_dim_id, remd_dimtype_var_id, old_mode

  NC_defineRemdIndices=.true.
  ! If ncid is < 0 we are not writing coords, exit cleanly
  if (ncid < 0) then
    NC_defineRemdIndices=.false.
    return
  endif
  ! Place the netcdf file back into define mode
  if (NC_error(nf90_redef(ncid), &
                   'putting netcdf file back into define mode')) return
  ! Define number of replica dimensions
  if (NC_error(nf90_def_dim(ncid, NCREMD_DIMENSION, remd_dimension, &
                              remd_dim_id),&
                 'define replica dimension')) return
  ! For each dimension, need to know the type
  dimensionID(1) = remd_dim_id
  remd_dimtype_var_id = NC_define_var(ncid, NCREMD_DIMTYPE, NF90_INT, &
                                      1, dimensionID)
  if (remd_dimtype_var_id.eq.-1) return ! 'define replica type for each dimension'
  ! Need to store the indices of replica in each dimension each frame
  if (isRestart) then
    indicesVID = NC_define_var(ncid, NCREMD_INDICES, NF90_INT, &
                                        1, dimensionID)
    ! For restarts only, store replica groups. FIXME: Necessary?
    if (.not.present(remd_groupsVID)) then
      write(mdout,'(a)') 'Internal Error: NC_defineRemdIndices called without remd_groupsVID'
      return
    endif
    remd_groupsVID = NC_define_var(ncid, NCREMD_GROUPS, NF90_INT, 1, dimensionID)
  else
    if (.not.present(frameDID)) then
      write(mdout,'(a)') 'Internal Error: NC_defineRemdIndices called without frameDID'
      return
    endif
    dimensionID(2) = frameDID 
    indicesVID = NC_define_var(ncid, NCREMD_INDICES, NF90_INT, &
                                        2, dimensionID)
  endif
  if (indicesVID.eq.-1) return ! 'define replica indices'
  ! Eventually store temperature tables for each T dimension
  
  ! Set NoFill and end definition mode
  if (NC_error(nf90_set_fill(ncid, NF90_NOFILL, old_mode), &
                    'set no fill')) return
  if (NC_error(nf90_enddef(ncid), 'end define')) return

  ! Store the type of each replica dimension
  if (NC_error(nf90_put_var(ncid, remd_dimtype_var_id, remd_types(:), &
                              start=(/ 1 /), count = (/ remd_dimension /)),&
                 'write replica type for each dimension')) return
  ! All is well
  NC_defineRemdIndices=.false.
end function NC_defineRemdIndices

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPRESTART
logical function NC_setupRestart(ncid, title, ncatom, coordVID, &
                                 velocityVID, TempVID, Time)
  use netcdf
  implicit none
  integer, intent(in)       :: ncid
  character(*), intent(out) :: title
  integer, intent(out)      :: ncatom
  integer, intent(out)      :: coordVID, velocityVID, TempVID
  double precision, intent(out) :: Time
  ! local
  integer timeVID

  NC_setupRestart=.true.
  ! Check this is actually a restart
  if (NC_checkConventions(ncid, .true.)) return
  ! Read title
  if (NC_GetTitle(ncid, title)) return 
  ! number of atoms, coordVID and velocityVID
  if (NC_setupCoordsVelo(ncid, ncatom, coordVID, velocityVID)) return
  ! Return a error if no coords
  if ( coordVID .eq. -1 ) then
    write(mdout,'(a)') 'Error: NetCDF file has no coordinates.'
    return
  endif 
  ! Setup time. Allowed to fail since not all restarts have this.
  if (NC_setupTime(ncid, timeVID)) then
    write(mdout,'(a)') 'Warning: NetCDF restart has no time value.'
    Time = 0.d0
    timeVID = -1
  else
    if (NC_error(nf90_get_var(ncid, timeVID, Time),&
                  "read_nc_restart(): Getting restart time")) return
  endif
  ! Temperature
  TempVID = NC_setupTemperature(ncid) 
  ! All is well
  NC_setupRestart=.false.
end function NC_setupRestart

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_READRESTARTBOX
!> @brief Read box information from previously opened netcdf restart.
logical function NC_readRestartBox(ncid,a,b,c,alpha,beta,gamma)
  use netcdf
  implicit none
  integer, intent(in)           :: ncid 
  double precision, intent(out) :: a,b,c,alpha,beta,gamma
  ! Local vars
  double precision, dimension(3) :: box
  integer :: cellLengthVID, cellAngleVID
 
  NC_readRestartBox=.true.
  call NC_setupBox(ncid, cellLengthVID, cellAngleVID)
  if (cellLengthVID.eq.-1 .or. cellAngleVID.eq.-1) return
  ! Cell lengths
  if (NC_error(nf90_get_var(ncid, cellLengthVID, box(1:3), &
                            start = (/ 1 /), count = (/ 3 /)),&
               'Getting box lengths')) return
  a = box(1)
  b = box(2)
  c = box(3)
  ! Cell angles
  if (NC_error(nf90_get_var(ncid, cellAngleVID, box(1:3), &
                            start = (/ 1 /), count = (/ 3 /)),&
               'Getting box angles')) return
  alpha = box(1)
  beta  = box(2)
  gamma = box(3)
  ! All is well 
  NC_readRestartBox=.false.
end function NC_readRestartBox

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_SETUPMULTID
!> @brief Determine if Netcdf file contains multi-D REMD info.
!> @return 0 if VIDs successfully read, 1 if no indexes present,
!!         -1 if error.
integer function NC_setupMultiD(ncid, expected_dim, indicesVID, groupsVID)
  use netcdf
  implicit none
  integer, intent(in)  :: ncid, expected_dim
  integer, intent(out) :: indicesVID, groupsVID
  ! local
  integer dimensionDID, remd_dimension, dimtypeVID

  NC_setupMultiD=1
  if ( nf90_inq_dimid(ncid, NCREMD_DIMENSION, dimensionDID) .ne. NF90_NOERR ) &
       return
  NC_setupMultiD=-1
  ! Although this is a second call to dimid, makes for easier code
  dimensionDID = GetDimInfo(ncid, NCREMD_DIMENSION, remd_dimension)
  if (dimensionDID.eq.-1) return
  ! Ensure valid # dimensions
  if (remd_dimension.ne.expected_dim) then
    write(mdout,'(2(a,i6),a)') 'Warning: # of dimensions in restart file (',&
                               remd_dimension, &
                               ') does not match current # (', expected_dim, ')'
    return
  endif
  ! Get dimension types
  if ( NC_error(nf90_inq_varid(ncid, NCREMD_DIMTYPE, dimtypeVID),&
                'Getting dimension type VID') ) return
  ! TODO: CHECK DIMENSION TYPES
  ! Get VID for replica indices
  if ( NC_error(nf90_inq_varid(ncid, NCREMD_INDICES, indicesVID),&
                'Getting indices VID') ) return
  ! Get VID for replica groups. FIXME: Is this necessary?
  if ( NC_error(nf90_inq_varid(ncid, NCREMD_GROUPS, groupsVID),&
                'Getting groups VID') ) return
  ! All is well
  NC_setupMultiD=0
end function NC_setupMultiD 

! ===== NO BINTRAJ SECTION =====================================================
#else /* NO BINTRAJ */
  public NC_NoNetcdfError, NC_checkRestart, NC_readRestartIndices, NC_checkTraj
contains
!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_NONETCDFERROR
!> @brief Should never be called if BINTRAJ defined.
subroutine NC_NoNetcdfError(mdout_unit)
  integer, intent(in) :: mdout_unit
  write(mdout_unit,'(a)') 'No binary trajectory support in this version.'
  write(mdout_unit,'(a)') 'Recompile using the -DBINTRAJ flag.'
end subroutine NC_NoNetcdfError
#endif /* BINTRAJ */

! ===== SHARED BINTRAJ / NO BINTRAJ ============================================
!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_CHECKRESTART
!> @return true if file is an amber netcdf restart.
!> @return false if not netcdf restart or when no BINTRAJ
logical function NC_checkRestart(filename)
# ifdef BINTRAJ
  use netcdf
# endif
  implicit none
  character(*), intent(in) :: filename
  ! local
  integer ncid, ierr
  NC_checkRestart=.false.
# ifdef BINTRAJ
  ierr = nf90_open( filename, NF90_NOWRITE, ncid )
  if (ierr .eq. NF90_NOERR) then
    if (.not.NC_checkConventions(ncid, .true.)) NC_checkRestart=.true.
    call NC_close(ncid)
  endif
# endif
end function NC_checkRestart

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_CHECKTRAJ
!> @return true if file is an amber netcdf trajectory.
!> @return false if not netcdf trajectory or when no BINTRAJ
logical function NC_checkTraj(filename)
# ifdef BINTRAJ
  use netcdf
# endif
  implicit none
  character(*), intent(in) :: filename
  ! local
  integer ncid, ierr
  NC_checkTraj=.false.
# ifdef BINTRAJ
  ierr = nf90_open( filename, NF90_NOWRITE, ncid )
  if (ierr .eq. NF90_NOERR) then
    if (.not.NC_checkConventions(ncid, .false.)) NC_checkTraj=.true.
    call NC_close(ncid)
  endif
# endif
end function NC_checkTraj

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_READRESTARTINDICES
!> @brief Read indices/groups from NC restart
!> @param filename Restart file to read indices from.
!> @param replica_indexes Array of replica indexes
!> @param group_num Array of REMD communicators I belong to
!> @param remd_dimension Dimension of replica_indexes
!> @return 0 if indices/groups successfully read, 1 if no indexes 
!!         present, -1 if error or no BINTRAJ.
integer function NC_readRestartIndices(filename, replica_indexes, group_num,&
                                       remd_dimension)
# ifdef BINTRAJ
  use netcdf
# endif
  implicit none
  character(*), intent(in)           :: filename
  integer, dimension(:), intent(out) :: replica_indexes
  integer, dimension(:), intent(out) :: group_num
  integer, intent(in)                :: remd_dimension 
  ! local
  integer ncid, indicesVID, groupsVID, ierr
  
  NC_readRestartIndices=-1
# ifdef BINTRAJ
  ! Open file - exit cleanly if not netcdf
  if (nf90_open( filename, NF90_NOWRITE, ncid ) .ne. NF90_NOERR) return
  ! Make sure its NetCDF restart.
  if (NC_checkConventions(ncid, .true.)) then
    call NC_close(ncid)
    return
  endif
  ! Check Dimension, setup indices and groups VID
  ierr = NC_setupMultiD(ncid, remd_dimension, indicesVID, groupsVID)
  if (ierr.eq.0) then
    ! Get indices
    if (NC_error(nf90_get_var(ncid, indicesVID, &
                              replica_indexes(1:remd_dimension), &
                              start = (/ 1 /), count = (/ remd_dimension /)),&
                 'Getting replica indices.')) return
    ! Get groups
    if (NC_error(nf90_get_var(ncid, groupsVID, &
                             group_num(1:remd_dimension), &
                             start = (/ 1 /), count = (/ remd_dimension /)),&
                'Getting group nums.')) return
    ! All is well
    NC_readRestartIndices=0
  else
    ! No indices present
    NC_readRestartIndices=1
    write(mdout,'(a)') '| Warning: NetCDF restart does not contain replica indices.'
  endif
  call NC_close(ncid)
# endif
end function NC_readRestartIndices

end module AmberNetcdf_mod
