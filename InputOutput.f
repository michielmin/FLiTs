	subroutine output(string)
	IMPLICIT NONE
	character string*(*)
	
	write(*,'(a)') trim(string)
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine ignorestar(un)
	IMPLICIT NONE
	integer un
	character c
1	read(un,fmt=3,end=2) c
	if(c.eq.'*') goto 1
	backspace(unit=un)
2	continue
3	format(a1)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
