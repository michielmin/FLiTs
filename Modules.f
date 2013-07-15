c=========================================================================================
c module containing the physical constants in cgs
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,G,Msun,AU,clight,Rsun
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(clight=2.9979245800d10) !cm/s
	parameter(AU=1.49598e13)
c	parameter(parsec=3.08568025e18)
	parameter(Rsun=6.955e10)
	parameter(Msun=1.98892e33)
c	parameter(Lsun=3.827e33)
c	parameter(kb=1.3806503d-16,sigma=5.6704d-5)
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	
	end module Constants

c=========================================================================================
c global setup for the FLiTs
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
c number of species
	integer MAXSPECIES,nspecies
	parameter(MAXSPECIES=20)
c stellar parameters
	real*8 Mstar,Rstar
c wavelength grid and resolution
	real*8 lmin,lmax,rlines,vresolution
	
c string converting functions
	character*20 int2string,dbl2string
	external int2string,dbl2string
	
c input files
	character*500 linefile(MAXSPECIES)
	character*500 structfile
c homogenous abundance
	real*8 abun_hom(MAXSPECIES)
c type of structure. 1=MCMax(LTE,homogeneous abundance), 2=ProDiMo
	integer structtype
c molecule names
	character*10 mol_name(MAXSPECIES)

c the grid setup. Note that we store cos(theta) in theta, but real theta in theta_av
	real*8,allocatable :: R(:),theta(:),R_av(:),theta_av(:)
	integer nR,nTheta,nlam
	real*8 Rin,Rout,inc
	real*8,allocatable :: lam_cont(:)
	
c the image grid
	integer nImR,nImPhi
	
c store all the blackbodies
	integer MAXT
	parameter(MAXT=5000)
	real*8,allocatable :: BB(:,:)	! dimensions nlam,MAXT

	type Line
		integer jup,jlow
		real*8 Aul,Blu,Bul
	end type Line

	type Molecule
		real*8,allocatable :: E(:),g(:) ! dimension is number of levels
c total mass of the molecule
		real*8 M
		type(Line),allocatable :: L(:) ! dimension is number of lines
	end type Molecule

c this type refers to a molecule with given T
	type MoleculeT
		real*8 dens,line_width
		real*8,allocatable :: npop(:) ! dimension is number of levels
		real*8,allocatable :: line_strength(:) ! dimension is number of lines
		integer nlev,nlines
	end type MoleculeT
		
c cell structure
	type Cell
		type(MoleculeT),allocatable :: Mol(:) ! dimension is number of species
c	Temperature and total gas density
		real*8 T,dens,v
		integer iT
c	Opacities and local radiation field. Opacities are given in units of tau/cm.
		real*8,allocatable :: kabs(:),albedo(:),kext(:),LRF(:) ! dimension is wavelength
	end type Cell

	type(Cell),allocatable,target :: C(:,:)	! dimension nR,nTheta
	
	type PathElement
		real*8 v,d,v1,v2
		type(Cell),pointer :: C
		type(PathElement),pointer :: next
		integer i,j
	end type PathElement

	type Path
c minimum and maximum velocity encountered in this path
		real*8 vmin,vmax
c number of elements
		integer n
c starting element
		type(PathElement),pointer :: start
c surface area of this path in the image and its coordinates
		real*8 A,R1,R2,phi1,phi2,R,phi
		real*8 vx,vy,vz,x,y,z
		real*8,allocatable :: flux_cont(:)
	end type Path
	
	type Tracer
		real*8 x,y,z,vx,vy,vz
		integer edgeNr,i,j
		logical onEdge
	end type Tracer

	
c==============================
	type(Molecule),allocatable :: Mol(:)
	
	type(Path),allocatable :: P(:,:)
	
	end module GlobalSetup

c=========================================================================================
	
