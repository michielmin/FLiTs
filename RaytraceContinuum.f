	subroutine RaytraceContinuum()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k
	real*8 flux,lam0,T,Planck,wl1,wl2,flux0
	
	allocate(BB(nlam,MAXT))
	do ilam=1,nlam
		do i=1,MAXT
			T=real(i)
			BB(ilam,i)=Planck(T,lam_cont(ilam))
		enddo
	enddo
	
	open(unit=20,file='out.dat')
	lam0=lmin
	ilam=1
	do while(lam0.lt.lmax)
		wl1=(lam_cont(ilam+1)-lam0)/(lam_cont(ilam+1)-lam_cont(ilam))
		wl2=1d0-wl1
		flux=0d0
		do i=1,nImR
			do j=1,nImPhi
				call TraceFlux(P(i,j),lam0,ilam,wl1,wl2,flux0)
				flux=flux+flux0*P(i,j)%A
			enddo
		enddo
		write(20,*) lam0,flux
		lam0=lam0+lam0/rlines
		do while(lam0.gt.lam_cont(ilam+1))
			ilam=ilam+1
		enddo
	enddo
	close(unit=20)
	
	return
	end
	
	
	subroutine TraceFlux(p0,lam0,ilam,wl1,wl2,flux)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 lam0,tau,exptau,flux,wl1,wl2,fact
	type(Path) p0

	fact=1d0
	flux=0d0
	do k=1,p0%n
		i=p0%i(k)
		j=p0%j(k)
		tau=(wl1*C(i,j)%kext(ilam)+wl2*C(i,j)%kext(ilam+1))*p0%d(k)
		if(tau.gt.1d-6) then
			exptau=exp(-tau)
		else
			exptau=1d0-tau
		endif
		flux=flux+(wl1*BB(ilam,C(i,j)%iT)*(1d0-C(i,j)%albedo(ilam))+
     &		wl2*BB(ilam+1,C(i,j)%iT)*(1d0-C(i,j)%albedo(ilam+1)))*(1d0-exptau)*fact
		fact=fact*exptau
		if(fact.lt.1d-6) exit
	enddo
	
	return
	end
	
