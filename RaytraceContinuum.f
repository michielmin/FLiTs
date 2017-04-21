	subroutine RaytraceContinuum()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k
	real*8 flux,lam0,T,Planck,wl1,wl2
	logical doit

	call output("==================================================================")
	call output("Raytracing the continuum")

	
	allocate(BB(nlam,MAXT))
	do ilam=1,nlam
		do i=1,MAXT
			T=real(i)
			BB(ilam,i)=Planck(T,lam_cont(ilam))
		enddo
	enddo
	
	do i=1,ngrids
	do j=1,npoints(i)
		allocate(P(i,j)%flux_cont(nlam))
	enddo
	enddo
	allocate(path2star%flux_cont(nlam))
	
	do ilam=1,nlam
		call tellertje(ilam-1,nlam-2)
		doit=.false.
		if(ilam.eq.1) then
			if(lam_cont(ilam+1).gt.lmin) doit=.true.
		else if(ilam.eq.nlam) then
			if(lam_cont(ilam-1).lt.lmax) doit=.true.
		else if(lam_cont(ilam-1).lt.lmax.and.lam_cont(ilam+1).gt.lmin) then
			doit=.true.
		endif
		if(doit) then
		flux=0d0
		do i=1,ngrids
		do j=1,npoints(i)
			call TraceFluxCont(P(i,j),ilam,P(i,j)%flux_cont(ilam))
			flux=flux+P(i,j)%flux_cont(ilam)*P(i,j)%A/real(ngrids)
		enddo
		enddo
		call Trace2StarCont(path2star,ilam,path2star%flux_cont(ilam))
		flux=flux+path2star%flux_cont(ilam)*path2star%A
		endif
	enddo
	
	return
	end
	
	
	subroutine TraceFluxCont(p0,ilam,flux)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 tau,exptau,flux,fact,S
	type(Path) p0

	fact=1d0
	flux=0d0
	do k=1,p0%n
		i=p0%i(k)
		j=p0%j(k)
		tau=C(i,j)%kext(ilam)*p0%d(k)
c		S=BB(ilam,C(i,j)%iT)*(1d0-C(i,j)%albedo(ilam))
c		S=S+C(i,j)%LRF(ilam)*C(i,j)%albedo(ilam)
c Using the source function for now.
		S=C(i,j)%S(ilam)
		if(tau.gt.1d-6) then
			exptau=exp(-tau)
		else
			exptau=1d0-tau
		endif
		flux=flux+S*(1d0-exptau)*fact
		fact=fact*exptau
		if(fact.lt.1d-6) exit
	enddo
	
	return
	end
	

	
	subroutine Trace2StarCont(p0,ilam,flux)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 tau,flux
	type(Path) p0

	tau=0d0
	do k=1,p0%n
		i=p0%i(k)
		if(i.eq.0) exit
		j=p0%j(k)
		tau=tau+C(i,j)%kext(ilam)*p0%d(k)
	enddo

	flux=Fstar(ilam)*exp(-tau)
	
	return
	end
	
