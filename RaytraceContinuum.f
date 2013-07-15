	subroutine RaytraceContinuum()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k
	real*8 flux,lam0,T,Planck,wl1,wl2
	
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
				call TraceFlux(P(i,j),ilam,P(i,j)%flux_cont(ilam))
				flux=flux+P(i,j)%flux_cont(ilam)*P(i,j)%A
			enddo
		enddo
		write(20,*) lam_cont(ilam),flux
	enddo
	close(unit=20)
	
	return
	end
	
	
	subroutine TraceFlux(p0,ilam,flux)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 tau,exptau,flux,fact
	type(Path) p0

	fact=1d0
	flux=0d0
	do k=1,p0%n
		i=p0%i(k)
		j=p0%j(k)
		tau=C(i,j)%kext(ilam)*p0%d(k)
		if(tau.gt.1d-6) then
			exptau=exp(-tau)
		else
			exptau=1d0-tau
		endif
		flux=flux+BB(ilam,C(i,j)%iT)*(1d0-C(i,j)%albedo(ilam))*(1d0-exptau)*fact
		fact=fact*exptau
		if(fact.lt.1d-6) exit
	enddo
	
	return
	end
	
