!=========================================================================================
! module containing the physical constants in cgs
!=========================================================================================
module Constants
  implicit none
  real(kind=8) :: pi, G, Msun, AU, clight, Rsun, mp, kb, hplanck, parsec
  parameter(pi=3.14159265358979323846264338328d0)
  parameter(clight=2.9979245800d10) !cm/s
  parameter(AU=1.4959787d13)
  parameter(parsec=3.08568025d18)
  parameter(Rsun=6.9599d10)
  parameter(Msun=1.9889225d33)
  !   parameter(Lsun=3.827e33)
  parameter(kb=1.38065812d-16)
  !   parameter(sigma=5.6704d-5)
  parameter(mp=1.67262178d-24)      ! proton mass
  parameter(G=6.6725985d-8)         ! in cm^3/g/s^2
  parameter(hplanck=6.62607554d-27) ! cm^2 g/s

end module Constants

!=========================================================================================
! global setup for the FLiTs
!=========================================================================================
module GlobalSetup
  implicit none
  ! stellar parameters
  real(kind=8) :: Mstar, Rstar, distance
  ! wavelength grid and resolution
  real(kind=8) :: lmin, lmax, rlines, vresolution, vres_mult, tau_max, vres_profile
  integer :: nvprofile, accuracy  ! accuracy can be 1,2,3 (increasing spatial sampling)
  logical :: doblend, cylindrical, LTE

  ! string converting functions
  !   character(len=20) :: int2string, dbl2string
  !   external int2string, dbl2string

  ! input files
  character(len=500) :: linefile(100)
  character(len=500) :: FLiTsfile, outputFile, output_lineFluxFile
  ! type of structure. 1=MCMax(LTE,homogeneous abundance), 2=ProDiMo
  integer :: structtype

  ! the grid setup. Note that we store cos(theta) in theta, but real theta in theta_av
  real(kind=8), allocatable :: R(:), theta(:, :), R_av(:), theta_av(:, :)
  real(kind=8), allocatable :: R_sphere(:), R_av_sphere(:)
  integer :: nR, nTheta, nlamCont, ilam1Cont, ilam2Cont, nmol, nlines, nblends, nspec
  real(kind=8) :: Rin, Rout, inc, Fstar_l
  real(kind=8), allocatable :: lam_cont(:), Fstar(:)

  ! the image grid
  integer :: nImR, ngrids, npoints_temp
  integer, allocatable :: nImPhi(:), npoints(:)
  real(kind=8) :: vmax

  ! the image cube
  real(kind=8), allocatable :: x_im(:, :, :), y_im(:, :, :), im_coord(:)
  ! FIXME: single precision should do it, to safe memory
  real(kind=8), allocatable :: imcube(:, :, :)
  integer npix, nvim, nlam_cube
  real(kind=8) :: dlam_cube ! the width of the channels in the image cube
  real(kind=8), allocatable :: lam_cube(:) ! the center wl, of the cube channels  
  ! produce an image cube
  logical imagecube
  ! use a more regular grid in the image planet (i.e. no random sampling)
  logical regular_grid
  character(len=100) :: imagecube_filename


  ! store all the blackbodies
  integer :: MAXT
  parameter(MAXT=5000)
  !   real(kind=8), allocatable :: BB(:,:)  ! dimensions nlam,MAXT

  ! random number generator
  real(kind=8) :: ran1
  external ran1
  integer :: idum
!$omp threadprivate(idum)
  ! line profile
  real(kind=8), allocatable :: profile(:)
  logical, allocatable :: profile_nz(:)

  character(len=10), allocatable :: mol_name0(:)

  type Line
    integer :: jup, jlow, imol
    real(kind=8) :: Aul, Blu, Bul, freq, lam, Eup
  end type Line

  type Blend
    type(Line), pointer :: L(:)
    real(kind=8), allocatable :: v(:)
    integer :: n, nvmin, nvmax
    type(Blend), pointer :: next
    integer, allocatable :: ib0(:), nb0(:)
    real(kind=8) :: lmin, lmax, lam, F
    logical :: computing, done
  end type Blend

  type Molecule
    real(kind=8), allocatable :: E(:), g(:) ! dimension is number of levels
    integer :: nlines, nlevels
    ! total mass of the molecule
    real(kind=8) :: M
    type(Line), allocatable :: L(:) ! dimension is number of lines
    character(len=10) :: name
    logical :: LTE
  end type Molecule

  type poplevels
    real(kind=8), allocatable :: N(:)
  end type poplevels

  ! cell structure
  type Cell
    ! Temperature and total gas density
    real(kind=8) :: Tdust, Tgas, dens, v
    ! properties of the molecule
    real(kind=8), allocatable :: line_width(:), N(:)
    type(poplevels), allocatable :: npop(:)
    integer, allocatable :: npopmax(:)
    !   real(kind=8), allocatable :: npop(:,:) ! dimension is nmol, number of levels
    real(kind=8), allocatable :: line_emis(:), line_abs(:)
    real(kind=8) :: kext_l, therm_l, scat_l
    integer :: iT
    ! Opacities and local radiation field. Opacities are given in units of tau/cm.
    real(kind=8), allocatable :: kabs(:), albedo(:), kext(:), LRF(:) ! dimension is wavelength
    real(kind=8), allocatable :: S(:) ! dimension is wavelength
  end type Cell

  type(Cell), allocatable, target :: C(:, :)  ! dimension nR,nTheta

  type Path
    ! minimum and maximum velocity encountered in this path
    real(kind=8), allocatable :: vmin(:), vmax(:)
    ! number of elements
    integer :: n
    integer, allocatable :: npopmax(:)
    ! surface area of this path in the image and its coordinates
    real(kind=8) :: A, R1, R2, phi1, phi2, R, phi
    real(kind=8) :: vx, vy, vz, x, y

    real(kind=8), allocatable :: v(:), d(:), v1(:), v2(:)
    integer, allocatable :: i(:), j(:)
    integer*2, allocatable :: im_ixy(:, :) ! first axis =2 (x,y), second axis can theoretically go to npix*2 (all pixels),
    ! but is likely much smaller FIXME don't know to estimate that, take npix for now
    integer*2 :: im_npix    ! nubmer of pixels associated to that path

    real(kind=8), allocatable :: flux_cont(:)  !continuum contribution at each wavelength
    real(kind=8), allocatable :: cont_contr(:) !continuum contribution at each path element
    real(kind=8), allocatable :: exptau_dust(:) !dust optical depth at each path element
    real(kind=8), allocatable :: S_dust(:) !dust source function at each path element
  end type Path

  type Tracer
    real(kind=8) :: x, y, z, vx, vy, vz
    integer :: edgeNr, i, j
    logical :: onEdge
  end type Tracer

  !==============================

  type(Path), allocatable, target :: P(:, :)
  type(Path), target :: path2star
  type(Molecule), allocatable :: Mol(:)
  type(Line), allocatable :: Lines(:)
  type(Blend), target :: Blends

end module GlobalSetup

!=========================================================================================

