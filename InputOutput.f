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


	character*20 function int2string(i,form)
	IMPLICIT NONE
	integer i
	character,intent(in),optional :: form*(*)
	
	if(form.ne.' ') then
		write(int2string,form) i
	else
		write(int2string,*) i
	endif
	
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	character*20 function dbl2string(x,form)
	IMPLICIT NONE
	real*8 x
	character,intent(in),optional :: form*(*)
	
	if(form.ne.' ') then
		write(dbl2string,form) x
	else
		write(dbl2string,*) x
	endif
	
	return
	end

