c=========================================================================================
c module containing the physical constants in cgs
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,G,Msun,AU,clight
	parameter(pi=3.1415926536)
	parameter(clight=2.9979245800d10) !cm/s
	
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
	real*8 Mstar
c wavelength grid and resolution
	real*8 lmin,lmax,rcont,rlines
	
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
		real*8 T,dens
c	Opacities and local radiation field. Opacities are given in units of tau/cm.
		real*8,allocatable :: kabs(:),ksca(:),kext(:),LRF(:) ! dimension is wavelength
	end type Cell
	
	type PathElement
		real*8 v,d
		type(Cell),pointer :: C
		type(Path),pointer :: next
	end type PathElement

	type Path
c minimum and maximum velocity encountered in this path
		real*8 vmin,vmax
c number of elements
		integer n
		type(PathElement),pointer :: start
	end type Path
	
	
c==============================
	type(Molecule),allocatable :: Mol(:)
	
	
	end module GlobalSetup

c=========================================================================================
	
