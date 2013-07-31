	subroutine RaytraceLines()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k,n1,n2,nl
	real*8 flux,lam,T,Planck,wl1,wl2,v,flux0
	
	open(unit=20,file='out.dat')

	lam=lmin
	ilam=1
	do while(lam.gt.lam_cont(ilam+1))
		ilam=ilam+1
	enddo
	v=-70d5
	nl=0
	do while(lam.lt.lmax)
		flux=0d0
		n1=0
		n2=0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,flux0)
!$OMP& SHARED(ilam,P,flux,nImR,nImPhi,n1,n2,v)
!$OMP DO
		do i=1,nImR
			do j=1,nImPhi
				if(v.gt.P(i,j)%vmin.and.v.lt.P(i,j)%vmax) then
					call TraceFluxLines(P(i,j),ilam,flux0,1d0)
					n1=n1+1
				else
					flux0=P(i,j)%flux_cont(ilam)
					n2=n2+1
				endif
				flux=flux+flux0*P(i,j)%A/2d0
				if(v.lt.-P(i,j)%vmin.and.v.gt.-P(i,j)%vmax) then
					call TraceFluxLines(P(i,j),ilam,flux0,-1d0)
					n1=n1+1
				else
					flux0=P(i,j)%flux_cont(ilam)
					n2=n2+1
				endif
				flux=flux+flux0*P(i,j)%A/2d0
			enddo
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		write(20,*) lam,flux,real(n1)/real(n1+n2)
		lam=lam*(1d0+1d0/rlines)
		do while(lam.gt.lam_cont(ilam+1))
			ilam=ilam+1
		enddo
		v=v+vresolution
		if(v.gt.70d5) then
			v=-70d5
			call output("another line" // int2string(nl,'(i4)'))
			nl=nl+1
		endif
	enddo
	close(unit=20)
	
	return
	end
	
	
	subroutine TraceFluxLines(p0,ilam,flux,vmult)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 tau,exptau,flux,fact,vmult
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
	
