	subroutine RaytraceContinuum()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam
	real*8 flux,T,Planck
	
	allocate(BB(nlam,MAXT))
	do ilam=1,nlam
		do i=1,MAXT
			T=real(i)
			BB(ilam,i)=Planck(T,lam_cont(ilam))
		enddo
	enddo
	
	do i=1,nImR
	do j=1,nImPhi
		allocate(P(i,j)%flux_cont(nlam))
	enddo
	enddo
	
	open(unit=20,file='out.dat')
	do ilam=1,nlam
		flux=0d0
		do i=1,nImR
			do j=1,nImPhi
				call TraceFluxCont(i,j,ilam,P(i,j)%flux_cont(ilam))
				flux=flux+P(i,j)%flux_cont(ilam)*P(i,j)%A
			enddo
		enddo
		write(20,*) lam_cont(ilam),flux
	enddo
	close(unit=20)
	
	return
	end
	
	
	subroutine TraceFluxCont(i,j,ilam,flux)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 lam0,tau,exptau,flux,wl1,wl2,fact
	type(PathElement),pointer :: current

	current => P(i,j)%start
	fact=1d0
	flux=0d0
	do k=1,P(i,j)%n-1
		tau=current%C%kext(ilam)*current%d
		if(tau.gt.1d-6) then
			exptau=exp(-tau)
		else
			exptau=1d0-tau
		endif
		flux=flux+BB(ilam,current%C%iT)*(1d0-current%C%albedo(ilam))*(1d0-exptau)*fact
		fact=fact*exptau
		if(fact.lt.1d-6) exit
		current => current%next
	enddo
	
	return
	end
	
