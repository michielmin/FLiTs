	subroutine output(string)
	IMPLICIT NONE
	character string(*)
	
	write(*,'(a)') trim(string)
	
	return
	end
	