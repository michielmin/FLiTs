c=========================================================================================
c module containing the physical constants in cgs
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,G,Msun,AU,clight,Rsun,mp,kb,hplanck,parsec
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(clight=2.9979245800d10) !cm/s
	parameter(AU=1.4959787d13)
	parameter(parsec=3.08568025d18)
	parameter(Rsun=6.9599d10)
	parameter(Msun=1.9889225d33)
c	parameter(Lsun=3.827e33)
	parameter(kb=1.38065812d-16)
c	parameter(sigma=5.6704d-5)
	parameter(mp=1.67262178d-24)  	  ! proton mass
	parameter(G=6.6725985d-8)         ! in cm^3/g/s^2
	parameter(hplanck=6.62607554d-27) ! cm^2 g/s

	end module Constants

c=========================================================================================
c global setup for the FLiTs
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
c stellar parameters
	real*8 Mstar,Rstar,distance
c wavelength grid and resolution
	real*8 lmin,lmax,rlines,vresolution,vres_mult,tau_max,vres_profile
	integer nvprofile,accuracy	! accuracy can be 1,2,3 (increasing spatial sampling)
	logical doblend,cylindrical,LTE

c string converting functions
c       character*20 int2string,dbl2string
c       external int2string,dbl2string

c input files
	character*500 linefile(100)
	character*500 FLiTsfile
c type of structure. 1=MCMax(LTE,homogeneous abundance), 2=ProDiMo
	integer structtype

c the grid setup. Note that we store cos(theta) in theta, but real theta in theta_av
	real*8,allocatable :: R(:),theta(:,:),R_av(:),theta_av(:,:)
	real*8,allocatable :: R_sphere(:),R_av_sphere(:)
	integer nR,nTheta,nlam,ilam1,ilam2,nmol,nlines,nblends,nspec
	real*8 Rin,Rout,inc,Fstar_l
	real*8,allocatable :: lam_cont(:),Fstar(:)

c the image grid
	integer nImR,ngrids,npoints_temp
	integer,allocatable :: nImPhi(:),npoints(:)
	real*8 vmax

c the image cube
	real*8,allocatable :: x_im(:,:,:),y_im(:,:,:)
	real*8,allocatable :: imcube(:,:,:)
	integer nint,npix,nvim
	logical imagecube

c store all the blackbodies
	integer MAXT
	parameter(MAXT=5000)
c	real*8,allocatable :: BB(:,:)	! dimensions nlam,MAXT

c random number generator
	real*8 ran1
	external ran1
	integer idum
!$omp threadprivate(idum)
c line profile
	real*8,allocatable :: profile(:)
	logical,allocatable :: profile_nz(:)

	character*10,allocatable :: mol_name0(:)

	type Line
		integer jup,jlow,imol
		real*8 Aul,Blu,Bul,freq,lam,Eup
	end type Line

	type Blend
		type(Line),pointer :: L(:)
		real*8,allocatable :: v(:)
		integer n,nvmin,nvmax
		type(Blend),pointer :: next
		integer,allocatable :: ib0(:),nb0(:)
		real*8 lmin,lmax,lam,F
		logical computing,done
	end type Blend

	type Molecule
		real*8,allocatable :: E(:),g(:) ! dimension is number of levels
		integer nlines,nlevels
c total mass of the molecule
		real*8 M
		type(Line),allocatable :: L(:) ! dimension is number of lines
		character*10 name
		logical LTE
	end type Molecule

	type poplevels
		real*8,allocatable :: N(:)
	end type poplevels
	integer,allocatable :: npop0(:)

c cell structure
	type Cell
c	Temperature and total gas density
		real*8 Tdust,Tgas,dens,v
c	properties of the molecule
		real*8,allocatable :: line_width(:),N(:)  ! dimension nmol
		real*8,allocatable :: N0(:),line_width0(:)
		type(poplevels),allocatable :: npop0(:)
		type(poplevels),allocatable :: npop(:)
		integer,allocatable :: npopmax(:)
c		real*8,allocatable :: npop(:,:) ! dimension is nmol, number of levels
		real*8,allocatable :: line_emis(:),line_abs(:)
		real*8 kext_l,therm_l,scat_l
		integer iT
c	Opacities and local radiation field. Opacities are given in units of tau/cm.
		real*8,allocatable :: kabs(:),albedo(:),kext(:),LRF(:) ! dimension is wavelength
		real*8,allocatable :: S(:) ! dimension is wavelength
	end type Cell

	type(Cell),allocatable,target :: C(:,:)	! dimension nR,nTheta

	type Path
c minimum and maximum velocity encountered in this path
		real*8,allocatable :: vmin(:),vmax(:)
c number of elements
		integer n
		integer,allocatable :: npopmax(:)
c surface area of this path in the image and its coordinates
		real*8 A,R1,R2,phi1,phi2,R,phi
		real*8 vx,vy,vz,x,y

		real*8,allocatable :: v(:),d(:),v1(:),v2(:)
		integer,allocatable :: i(:),j(:)

		real*8,allocatable :: flux_cont(:) !continuum contribution at each wavelength
		real*8,allocatable :: cont_contr(:) !continuum contribution at each path element
		real*8,allocatable :: exptau_dust(:) !dust optical depth at each path element
		real*8,allocatable :: S_dust(:) !dust source function at each path element
	end type Path

	type Tracer
		real*8 x,y,z,vx,vy,vz
		integer edgeNr,i,j
		logical onEdge
	end type Tracer


c==============================

	type(Path),allocatable,target :: P(:,:)
	type(Path),target :: path2star
	type(Molecule),allocatable :: Mol(:)
	type(Line),allocatable :: Lines(:)
	type(Blend),target :: Blends

	end module GlobalSetup

c=========================================================================================

