	subroutine ComputeLTE()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilines,k,i_low,i_up,imol
	real*8 UT,T
	
	call output("Computing LTE level populations")
	
	do i=0,nR
	call tellertje(i+1,nR+1)
	do j=1,nTheta
		T=C(i,j)%Tgas
		if(T.lt.3d0) T=3d0
		do imol=1,nmol
			if(Mol(imol)%LTE) then
				UT=0d0
				do k=1,Mol(imol)%nlevels
					UT=UT+Mol(imol)%g(k)*exp(-Mol(imol)%E(k)/T)
				enddo
		
				do k=1,Mol(imol)%nlevels
					C(i,j)%npop(imol,k)=Mol(imol)%g(k)*exp(-Mol(imol)%E(k)/T)/UT
				enddo
			endif
		enddo
	enddo
	enddo
	
	return
	end
	

	