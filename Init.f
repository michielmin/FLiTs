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
			read(value,*) Mstar
		case("rstar")
			read(value,*) Rstar
		case("structfile")
			structfile=value
		case("structtype")
			read(value,*) structtype
		case("linefile")
			read(value,*) linefile
		case("mass_mol")
			read(value,*) Mol%M
		case("lte")
			read(value,*) LTE
		case("lmin")
			read(value,*) lmin
		case("lmax")
			read(value,*) lmax
		case("vres")	! given in cm/s
			read(value,*) vresolution
		case("vres_mult")	! given in cm/s
			read(value,*) vres_mult
		case("tau_max")	! given in cm/s
			read(value,*) tau_max
		case("inc")	! with respect to pole on
			read(value,*) inc
		case default
			call output("Unknown keyword: " // trim(key))
			stop
	end select

c read another command, so go back
	goto 10

30	continue
	close(unit=20)

c everything is read in now
c output the setup to the screen and the log file
	call output("==================================================================")

	rlines=1d0/(sqrt((1d0+vresolution/clight)/(1d0-vresolution/clight))-1d0)

	call output("Mass of the star:   "//trim(dbl2string(Mstar,'(f13.4)'))//" Msun")
	call output("Minimum wavelength: "//trim(dbl2string(lmin,'(f13.4)'))//" micron")
	call output("Maximum wavelength: "//trim(dbl2string(lmax,'(f13.4)'))//" micron")
	call output("Resolution lines:   "//trim(dbl2string(rlines,'(f13.4)'))//" (dlam/lam)")
	call output("Velocity resolution:"//trim(dbl2string(vresolution/1d5,'(f13.4)'))//" km/s")
	call output("Inclination angle:  "//trim(dbl2string(inc,'(f13.4)'))//" degrees")

	call output("==================================================================")
	call output("Line file: "//trim(linefile))

	
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
	structfile='forFLiTs.fits.gz'
	structtype=1
	lmin=5
	lmax=50
	inc=35d0
	vresolution=1d5	! given in cm/s
	vres_mult=10d0
	tau_max=15d0
		
	Mol%M=28	! default is CO
	LTE=.true.	! default is LTE for now
	
	return
	end
	
