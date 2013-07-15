	subroutine RaytraceContinuum()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam
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
				call TraceFlux(i,j,lam0,ilam,wl1,wl2,flux0)
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
	
	
	subroutine TraceFlux(i,j,lam0,ilam,wl1,wl2,flux)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 lam0,tau,exptau,flux,wl1,wl2,fact
	type(PathElement),pointer :: current

	current => P(i,j)%start
	fact=1d0
	flux=0d0
	do k=1,P(i,j)%n
		tau=(wl1*current%C%kext(ilam)+wl2*current%C%kext(ilam+1))*current%d
		if(tau.gt.1d-6) then
			exptau=exp(-tau)
		else
			exptau=1d0-tau
		endif
		flux=flux+(wl1*BB(ilam,current%C%iT)*(1d0-current%C%albedo(ilam))+
     &		wl2*BB(ilam+1,current%C%iT)*(1d0-current%C%albedo(ilam+1)))*(1d0-exptau)*fact
		fact=fact*exptau
		if(fact.lt.1d-6) exit
		current => current%next
	enddo
	
	return
	end
	
