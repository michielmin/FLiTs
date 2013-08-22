	subroutine ComputeLTE()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilines,k,i_low,i_up
	real*8 UT,T
	
	do i=0,nR
	do j=1,nTheta
		T=C(i,j)%Tgas
		if(T.lt.3d0) T=3d0
		UT=0d0
		do k=1,Mol%nlevels
			UT=UT+Mol%g(k)*exp(-Mol%E(k)/T)
		enddo
		
		allocate(C(i,j)%npop(Mol%nlevels))

		do k=1,Mol%nlevels
			C(i,j)%npop(k)=Mol%g(k)*exp(-Mol%E(k)/T)/UT
		enddo

		do ilines=1,Mol%nlines
			i_low=Mol%L(ilines)%jlow
			i_up=Mol%L(ilines)%jup
	
			Mol%L(ilines)%Bul=Mol%L(ilines)%Aul*2d0*hplanck*Mol%L(ilines)%freq**3/(clight**2)
			Mol%L(ilines)%Blu=Mol%L(ilines)%Bul*Mol%g(i_up)/Mol%g(i_low)
		enddo
	enddo
	enddo
	
	return
	end
	

	