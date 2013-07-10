c=========================================================================================
c This subroutine reads in the input file. In the input file all other input files are
c specified.
c=========================================================================================
	subroutine Initialize()
	use GlobalSetup
	IMPLICIT NONE
	integer ncla	! number of command line arguments
	character*1000 line,inputfile
	character*100 key,value

	call SetDefaults()
	
c set the number of species to 0
	nspecies=0

	call getarg(1,inputfile)
	open(unit=20,file=inputfile,RECL=1000)
	
	ncla=1
10	continue
	if(ncla.eq.1) then
		call ignorestar(20)
		read(20,'(a1000)',end=20,err=20) line
	else
20		call getarg(1+ncla,line)
		if(line(1:2).eq.'-s') then
			ncla=ncla+1
			call getarg(2+ncla,line)
			call output("Command line argument: " // trim(line))
			ncla=ncla+1
		else
			if(line.ne.' ') then
c				try to read another command line argument
				goto 20
			else
c				all arguments are read
				goto 30
			endif
		endif
	endif

	call get_key_value(line,key,value)

	select case(key)
		case("mstar")
			read(value,*) mstar
		case("structfile")
			structfile=value
		case("structtype")
			read(value,*) structtype
		case("linefile")
			nspecies=nspecies+1
			read(value,*) linefile(nspecies)
		case("lmin")
			read(value,*) lmin
		case("lmax")
			read(value,*) lmax
		case("rcont")
			read(value,*) rcont
		case("rlines")
			read(value,*) rlines
		case default
			call output("Unknown keyword: " // trim(key))
			stop
	end select

c read another command, so go back
	goto 10

30	continue
	close(unit=20)

	
	return
	end

c=========================================================================================
c This subroutine just seperates the key and value component of a string given
c key=value syntax. Key is transformed to lowercase.
c=========================================================================================
	subroutine get_key_value(line,key,value)
	IMPLICIT NONE
	character line*(*)
	character*100 key,value
	
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

