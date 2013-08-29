	subroutine PrepareStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	
	do i=0,nR
		do j=1,nTheta
			if(popfile.eq.' ') C(i,j)%abun=1d-4
			C(i,j)%N=C(i,j)%dens*C(i,j)%abun/(Mol%M*mp)
		enddo
	enddo
	
	return
	end

