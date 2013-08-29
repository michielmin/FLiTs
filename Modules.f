c=========================================================================================
c module containing the physical constants in cgs
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,G,Msun,AU,clight,Rsun,mp,kb,hplanck,parsec
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(clight=2.9979245800d10) !cm/s
	parameter(AU=1.49598e13)
	parameter(parsec=3.08568025e18)
	parameter(Rsun=6.955e10)
	parameter(Msun=1.98892e33)
c	parameter(Lsun=3.827e33)
	parameter(kb=1.3806503d-16)
c	parameter(sigma=5.6704d-5)
	parameter(mp=1.67262178d-24)	!proton mass
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	parameter(hplanck=6.626068e-27) ! cm^2 g/s
	
	end module Constants

c=========================================================================================
c global setup for the FLiTs
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
c stellar parameters
	real*8 Mstar,Rstar,distance
c wavelength grid and resolution
	real*8 lmin,lmax,rlines,vresolution,vres_mult,tau_max
	integer nvprofile
	
c string converting functions
	character*20 int2string,dbl2string
	external int2string,dbl2string
	
c input files
	character*500 linefile
	character*500 structfile
	character*500 popfile
c type of structure. 1=MCMax(LTE,homogeneous abundance), 2=ProDiMo
	integer structtype,LTE

c the grid setup. Note that we store cos(theta) in theta, but real theta in theta_av
	real*8,allocatable :: R(:),theta(:),R_av(:),theta_av(:)
	integer nR,nTheta,nlam,ilam1,ilam2
	real*8 Rin,Rout,inc,Fstar_l
	real*8,allocatable :: lam_cont(:),Fstar(:)
	
c the image grid
	integer nImR,nImPhi
	real*8 vmax
	
c store all the blackbodies
	integer MAXT
	parameter(MAXT=5000)
	real*8,allocatable :: BB(:,:)	! dimensions nlam,MAXT

	type Line
		integer jup,jlow
		real*8 Aul,Blu,Bul,freq,lam,Eup
	end type Line

	type Molecule
		real*8,allocatable :: E(:),g(:) ! dimension is number of levels
		integer nlines,nlevels
c total mass of the molecule
		real*8 M
		type(Line),allocatable :: L(:) ! dimension is number of lines
		character*10 name
	end type Molecule
		
c cell structure
	type Cell
c	Temperature and total gas density
		real*8 Tdust,Tgas,dens,v,N,abun
c	properties of the molecule
		real*8 line_width
		real*8,allocatable :: profile(:)
		logical,allocatable :: profile_nz(:)
		real*8,allocatable :: npop(:) ! dimension is number of levels
		real*8 line_emis,line_abs
		real*8 kext_l,albedo_l,BB_l,LRF_l
		integer iT
c	Opacities and local radiation field. Opacities are given in units of tau/cm.
		real*8,allocatable :: kabs(:),albedo(:),kext(:),LRF(:) ! dimension is wavelength
	end type Cell

	type(Cell),allocatable,target :: C(:,:)	! dimension nR,nTheta
	
	type Path
c minimum and maximum velocity encountered in this path
		real*8 vmin,vmax
c number of elements
		integer n,npopmax
c surface area of this path in the image and its coordinates
		real*8 A,R1,R2,phi1,phi2,R,phi
		real*8 vx,vy,vz,x,y,z

		real*8,allocatable :: v(:),d(:),v1(:),v2(:)
		integer,allocatable :: i(:),j(:)
		
		real*8,allocatable :: flux_cont(:)
	end type Path
	
	type Tracer
		real*8 x,y,z,vx,vy,vz
		integer edgeNr,i,j
		logical onEdge
	end type Tracer

	
c==============================
	
	type(Path),allocatable,target :: P(:,:)
	type(Path),target :: path2star
	type(Molecule) Mol
	
	end module GlobalSetup

c=========================================================================================
	
