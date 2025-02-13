c=========================================================================================
c This subroutine reads in the input file. In the input file all other input files are
c specified.
c=========================================================================================
	subroutine Initialize()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ncla	! number of command line arguments
	character*1000 readline,inputfile
	character*100 key,value
	integer i

	call SetDefaults()
	
	call getarg(1,inputfile)
	open(unit=20,file=inputfile,RECL=1000)
	
	ncla=-1

	nmol=0
	
	idum=-42
	
20	ncla=ncla+1
10	continue
	if(ncla.eq.0) then
		call ignorestar(20)
		read(20,'(a1000)',end=20,err=20) readline
	else
		call getarg(1+ncla,readline)
		if(readline(1:2).eq.'-s') then
			ncla=ncla+1
			call getarg(1+ncla,readline)
			call output("Command line argument: " // trim(readline))
			ncla=ncla+1
		else
			if(readline.ne.' ') then
c				try to read another command line argument
				ncla=ncla+1
				goto 10
			else
c				all arguments are read
				goto 30
			endif
		endif
	endif

	if(readline.eq.' ') goto 10

	call get_key_value(readline,key,value)

	select case(key)
		case("mstar")
			call output("Mstar keyword no longer supported")
			call output("value is read from the fits file")
			stop
		case("rstar")
			read(value,*) Rstar
			call output("Rstar keyword no longer supported")
			call output("value is read from the fits file")
			stop
		case("structfile")
			call output("Please use forFLiTs file")
			stop
		case("linefile")
			nmol=nmol+1
			linefile(nmol)=value
		case("popfile")
			call output("Please use forFLiTs file")
			stop
		case("flitsfile")
			FLiTsfile=value
		case("lte")
			read(value,*) LTE
		case("blend")
			read(value,*) doblend
		case("cylindrical")
			read(value,*) cylindrical
		case("lmin")
			read(value,*) lmin
		case("lmax")
			read(value,*) lmax
		case("vres")	! given in cm/s
			read(value,*) vresolution
		case("vres_mult")
			call output("vres_mult no longer supported!")
			stop
			read(value,*) vres_mult
		case("vres_profile")	! given in cm/s
			read(value,*) vres_profile
		case("tau_max")	! given in cm/s
			read(value,*) tau_max
		case("inc")	! with respect to pole on
			read(value,*) inc
		case("accuracy")
			read(value,*) accuracy
		case("imagecube","imcube")
			read(value,*) imagecube
		case("imagecube_npix","imcube_npix")
			read(value,*) npix
		case("idum","seed")
			read(value,*) idum
		case default
			call output("Unknown keyword: " // trim(key))
			stop
	end select

c read another command, so go back
	goto 10

30	continue
	close(unit=20)

	vres_mult=vresolution/vres_profile

	if(abs(inc).lt.1d0) then
		call output("Increasing inclination to 1 degrees")
		inc=1
	endif

	call output("==================================================================")
	
	return
	end


c=========================================================================================
c This subroutine just seperates the key and value component of a string given
c key=value syntax. Key is transformed to lowercase.
c=========================================================================================
	subroutine get_key_value(line,key,value)
	IMPLICIT NONE
	character*1000 line
	character*100 key,value
	integer i
	
	key=line(1:index(line,'=')-1)
	value=line(index(line,'=')+1:len_trim(line))
	if(value(1:1).eq.'"'.or.value(1:1).eq."'") then
		value=value(2:len_trim(value)-1)
	endif
	do i=1,len_trim(key)
		if(iachar(key(i:i)).ge.65.and.iachar(key(i:i)).le.90) then
			key(i:i)=achar(iachar(key(i:i))+32)
		endif
	enddo

	return
	end
	
c=========================================================================================
c This subroutine sets the default values for the global variables
c=========================================================================================
	subroutine SetDefaults()
	use GlobalSetup
	IMPLICIT NONE
	
	Mstar=1d0
	Rstar=1d0
	FLiTsfile='ProDiMoForFLiTs.fits.gz'
	lmin=5
	lmax=50
	inc=30d0
	vresolution=1d5	! given in cm/s
	vres_profile=1d4
	tau_max=15d0

	idum=-42
	
	accuracy=1
		
	LTE=.false.	! default is non-LTE

	cylindrical=.true.
	
	doblend=.true.
	
	! imagecube false by default
	imagecube=.false.
	npix=0 ! FIXME: not check yet for the case of imagecube=true
	
	return
	end
	
