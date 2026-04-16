!=========================================================================================
! This subroutine reads in the input file. In the input file all other input files are
! specified.
!=========================================================================================
subroutine Initialize()
  use GlobalSetup
  use Constants
  use InOut
  implicit none
  integer :: ncla        ! number of command line arguments
  character(len=1000) :: readline, inputfile
  character(len=100) :: key, value

  call SetDefaults()

  call getarg(1, inputfile)
  open (unit=20, file=inputfile, RECL=1000)

  ncla = -1

  nmol = 0

  idum = -42

20 ncla = ncla + 1
10 continue
  if (ncla .eq. 0) then
    call ignorestar(20)
    read (20, '(a1000)', end=20, err=20) readline
  else
    call getarg(1 + ncla, readline)
    if (readline(1:2) .eq. '-s') then
      ncla = ncla + 1
      call getarg(1 + ncla, readline)
      call output("Command line argument: "//trim(readline))
      ncla = ncla + 1
    else
      if (readline .ne. ' ') then
!                                try to read another command line argument
        ncla = ncla + 1
        goto 10
      else
!                                all arguments are read
        goto 30
      end if
    end if
  end if

  if (readline .eq. ' ') goto 10

  call get_key_value(readline, key, value)

  select case (key)
  case ("mstar")
    call output("Mstar keyword no longer supported")
    call output("value is read from the fits file")
    stop
  case ("rstar")
    read (value, *) Rstar
    call output("Rstar keyword no longer supported")
    call output("value is read from the fits file")
    stop
  case ("structfile")
    call output("Please use forFLiTs file")
    stop
  case ("linefile")
    nmol = nmol + 1
    linefile(nmol) = value
  case ("popfile")
    call output("Please use forFLiTs file")
    stop
  case ("flitsfile")
    FLiTsfile = value
  case ("outputfile")
    outputFile = value
  case ("outputlinefluxfile")
    output_lineFluxFile = value
  case ("lte")
    read (value, *) LTE
  case ("blend")
    read (value, *) doblend
  case ("cylindrical")
    read (value, *) cylindrical
  case ("lmin")
    read (value, *) lmin
  case ("lmax")
    read (value, *) lmax
  case ("vres")        ! given in cm/s
    read (value, *) vresolution
  case ("vres_mult")
    call output("vres_mult no longer supported!")
    stop
    read (value, *) vres_mult
  case ("vres_profile")        ! given in cm/s
    read (value, *) vres_profile
  case ("tau_max")        ! given in cm/s
    read (value, *) tau_max
  case ("inc")        ! with respect to pole on
    read (value, *) inc
  case ("accuracy")
    read (value, *) accuracy
  case ("imagecube", "imcube")
    read (value, *) imagecube
  case("imagecube_npix","imcube_npix")
    read(value,*) npix
  case ("idum", "seed")
    read (value, *) idum
  case ("ngrids")  ! number of grids to use for the line RT, Default is 5
    read (value, *) ngrids
  case ("regular_grid")
    read (value, *) regular_grid
  case default
    call output("Unknown keyword: "//trim(key))
    stop
  end select

! read another command, so go back
  goto 10

30 continue
  close (unit=20)

  vres_mult = vresolution/vres_profile

  if (abs(inc) .lt. 1d0) then
    call output("Increasing inclination to 1 degrees")
    inc = 1
  end if

  if (regular_grid) then 
    call output("Using a regular grid and set ngrids=1")
    ngrids = 1
  end if

  call output("==================================================================")

  return
end subroutine Initialize

!=========================================================================================
! This subroutine just seperates the key and value component of a string given
! key=value syntax. Key is transformed to lowercase.
!=========================================================================================
subroutine get_key_value(line, key, value)
  implicit none
  character(len=1000) :: line
  character(len=100) :: key, value
  integer :: i

  key = line(1:index(line, '=') - 1)
  value = line(index(line, '=') + 1:len_trim(line))
  if (value(1:1) .eq. '"' .or. value(1:1) .eq. "'") then
    value = value(2:len_trim(value) - 1)
  end if
  do i = 1, len_trim(key)
    if (iachar(key(i:i)) .ge. 65 .and. iachar(key(i:i)) .le. 90) then
      key(i:i) = achar(iachar(key(i:i)) + 32)
    end if
  end do

  return
end subroutine get_key_value

!=========================================================================================
! This subroutine sets the default values for the global variables
!=========================================================================================
subroutine SetDefaults()
  use GlobalSetup
  implicit none

  Mstar = 1d0
  Rstar = 1d0
  FLiTsfile = 'ProDiMoForFLiTs.fits'
  outputFile = 'specFLiTs.out'
  output_lineFluxFile = 'lineFlux_FLiTs.out'
  lmin = 5
  lmax = 50
  inc = 30d0
  vresolution = 1d5        ! given in cm/s
  vres_profile = 1d4
  tau_max = 15d0
  idum = -42
  accuracy = 1
  LTE = .false.        ! default is non-LTE
  cylindrical = .true.
  doblend = .true.
  imagecube = .false.
  regular_grid = .false.
  npix=0 ! FIXME: no check yet for the case of imagecube=true
  ngrids = 5
  return
end subroutine SetDefaults

subroutine WriteSysinfo()
! Write some info about the system, openmp etc. into the log file
  use iso_fortran_env, ONLY: compiler_version, compiler_options
  use InOut
  implicit none
  character(len=100) :: host, folder
  integer :: numPE, OK
  integer, external :: OMP_GET_NUM_THREADS, omp_get_proc_bind
  integer, external :: omp_get_num_places
  integer :: ompversion, OMP_PROC_BIND, OMP_PLACES
  real*4 :: cfitsioversion

  write (*, *)
  write (*, *) "WRITE_SYSINFO: ..."
  call GET_ENVIRONMENT_VARIABLE(NAME="HOSTNAME", VALUE=host, STATUS=OK)

  if (ok .eq. 0) then ! does not always exist
    call output("running on host "//trim(host))
  end if

  call GET_ENVIRONMENT_VARIABLE(NAME="PWD", VALUE=folder, STATUS=OK)
  call output("running in folder "//trim(folder))

!$omp parallel
#ifdef _OPENMP
  numPE = OMP_GET_NUM_THREADS()
  ompversion = _OPENMP
  OMP_PROC_BIND = omp_get_proc_bind()
  OMP_PLACES = omp_get_num_places()
#else
  numPE = 1
  ompversion = 0
  OMP_PROC_BIND = -1
  OMP_PLACES = -1
#endif
!$omp end parallel
  
  call output(" using" // trim(int2string(numPE,'(i3)')) // " processors, with PROC_BIND=" // trim(int2string(OMP_PROC_BIND,'(i3)')) // " and PLACES=" // trim(int2string(OMP_PLACES,'(i3)')))
  if (ompversion > 0) then
    call output(" using OMP version: "//trim(int2string(ompversion, '(i7)')))
  end if

  call output(" compiler version: "//trim(compiler_version()))
  call output(" compiler options: "//trim(compiler_options()))

  ! check for required cfitsio version
  call ftvers(cfitsioversion)
  call output(" using cfitsio version: "//trim(dbl2string(cfitsioversion*1.d0, '(f7.5)')))

  ! This should not happen, because the compilation should fail anyway
  ! Hoever, keep it in case somebody is copying prodimo binaries around
  if (cfitsioversion .lt. 4.0199) then ! cfitsio return 4.01999 for 4.2.0
    call output("*** ERROR: FLiTs requires cfitsio version>=4.2.0")
    stop
  end if

end subroutine WriteSysinfo

