	subroutine PrepareStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,imol
	
	do i=0,nR
		C(i,0)%dens=1d-50
	enddo
	do j=0,nTheta
		C(0,j)%dens=1d-50
		if(cylindrical) C(nR,j)%dens=1d-50
	enddo
	do i=0,nR
		do j=0,nTheta
			allocate(C(i,j)%N(nmol))
			do imol=1,nmol
				C(i,j)%N(imol)=C(i,j)%dens*C(i,j)%abun(imol)/(Mol(imol)%M*mp)
			enddo
		enddo
	enddo
	
	
	
	return
	end

