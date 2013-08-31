	subroutine PrepareStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,imol
	
	do i=0,nR
		do j=1,nTheta
			if(popfile.eq.' ') C(i,j)%abun(1:nmol)=1d-4
			allocate(C(i,j)%N(nmol))
			do imol=1,nmol
				C(i,j)%N(imol)=C(i,j)%dens*C(i,j)%abun(imol)/(Mol(imol)%M*mp)
			enddo
		enddo
	enddo
	
	return
	end

