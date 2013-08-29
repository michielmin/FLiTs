	subroutine ComputeLTE()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilines,k,i_low,i_up
	real*8 UT,T
	
	call output("Computing LTE level populations")
	
	do i=0,nR
	do j=1,nTheta
		T=C(i,j)%Tgas
		if(T.lt.3d0) T=3d0
		UT=0d0
		do k=1,Mol%nlevels
			UT=UT+Mol%g(k)*exp(-Mol%E(k)/T)
		enddo
		
		if(.not.allocated(C(i,j)%npop)) allocate(C(i,j)%npop(Mol%nlevels))

		do k=1,Mol%nlevels
			C(i,j)%npop(k)=Mol%g(k)*exp(-Mol%E(k)/T)/UT
		enddo
	enddo
	enddo
	
	return
	end
	

	