	subroutine RaytraceLines()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k,ilines,vmult,iv,nv,nl
	real*8 lam,T,Planck,wl1,wl2,v,flux0,starttime,stoptime,tot,fact
	real*8,allocatable :: flux(:)
	type(Path),pointer :: PP
	type(Line) :: LL
	
	call output("==================================================================")
	call output("Tracing the lines")

	call cpu_time(starttime)

	open(unit=20,file='out.dat')

	nv=int(vmax*1.5/vresolution)

	allocate(flux(-nv:nv))

	lam=lmin
	ilam=1
	do while(lam.gt.lam_cont(ilam+1))
		ilam=ilam+1
	enddo
	do i=0,nR
		do j=1,nTheta
			allocate(C(i,j)%profile(-nv:nv))
		enddo
	enddo

	do i=0,nR
		do j=1,nTheta
			do k=-nv,nv
				C(i,j)%profile(k)=((real(k)*vresolution)/C(i,j)%line_width)**2
			enddo
			C(i,j)%profile=clight*exp(-C(i,j)%profile)/(C(i,j)%line_width*sqrt(pi))
		enddo
	enddo

	nl=0
	do ilines=1,Mol%nlines
		call tellertje(ilines,Mol%nlines)
		flux=0d0

		LL = Mol%L(ilines)		
		lam=clight*1d4/(LL%freq)
		if(lam.gt.lmin.and.lam.lt.lmax) then
		nl=nl+1

		ilam=ilam1
		do while(lam.gt.lam_cont(ilam+1))
			ilam=ilam+1
		enddo
		wl1=(lam_cont(ilam+1)-lam)/(lam_cont(ilam+1)-lam_cont(ilam))
		wl2=1d0-wl1

		do i=0,nR
			do j=1,nTheta
				C(i,j)%kext_l=wl1*C(i,j)%kext(ilam)+wl2*C(i,j)%kext(ilam+1)
				C(i,j)%albedo_l=wl1*C(i,j)%albedo(ilam)+wl2*C(i,j)%albedo(ilam+1)
				C(i,j)%BB_l=wl1*BB(ilam,C(i,j)%iT)+wl2*BB(ilam+1,C(i,j)%iT)

				if(C(i,j)%npop(LL%jup).gt.0d0) then
					C(i,j)%line_emis=C(i,j)%npop(LL%jup)*LL%Aul/(C(i,j)%npop(LL%jlow)*LL%Blu-C(i,j)%npop(LL%jup)*LL%Bul)
				else
					C(i,j)%line_emis=0d0
				endif
				fact=hplanck*C(i,j)%N/(4d0*pi)
				C(i,j)%line_abs=fact*(C(i,j)%npop(LL%jlow)*LL%Blu-C(i,j)%npop(LL%jup)*LL%Bul)
			enddo
		enddo

		do i=1,nImR
!$OMP PARALLEL IF(.false.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(flux0,iv,j,PP,vmult)
!$OMP& SHARED(i,ilam,P,flux,nImR,nImPhi,nv,vresolution,wl1,wl2)
!$OMP DO
			do j=1,nImPhi
				PP => P(i,j)
				vmult=1
				do iv=-nv,nv
					if(real(iv*vresolution).gt.PP%vmax.or.real(iv*vresolution).lt.PP%vmin) then
						flux0=wl1*PP%flux_cont(ilam)+wl2*PP%flux_cont(ilam+1)
					else
						call TraceFluxLines(PP,flux0,iv,nv,vmult)
					endif
					do vmult=-1,1,2
						flux(iv*vmult)=flux(iv*vmult)+flux0*PP%A/2d0
					enddo
				enddo
			enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		enddo
		do i=-nv,nv
			write(20,*) lam*sqrt((1d0+real(i)*vresolution/clight)/(1d0-real(i)*vresolution/clight)),flux(i),real(i)*vresolution/1d5
		enddo
		
		endif

	enddo
	close(unit=20)
	
	call cpu_time(stoptime)

	call output("Time used for the lines:"//trim(dbl2string(stoptime-starttime,'(f8.2)'))
     &			//" s")
	call output("Time used per line:     "//trim(dbl2string((stoptime-starttime)/real(nl),'(f8.2)'))
     &			//" s")


	return
	end
	
	
	subroutine TraceFluxLines(p0,flux,ii,nv,vmult)
	use GlobalSetup
	use Constants
	integer i,j,k,vmult,iv,ii,nv
	real*8 tau,exptau,flux,fact,profile,S,tau_gas,tau_dust,tau_d
	type(Path) p0
	type(Cell),pointer :: CC

	fact=1d0
	flux=0d0

	do k=1,p0%n
		i=p0%i(k)
		if(i.ne.0) then
			j=p0%j(k)
			CC => C(i,j)
			iv=int(p0%v(k)/vresolution)
			jj=ii-iv
			if(jj.lt.-nv) jj=-nv
			if(jj.gt.nv) jj=nv
			profile=CC%profile(jj)

			tau_dust=CC%kext_l
			tau_gas=profile*CC%line_abs
			tau=tau_dust+tau_gas

c	dust thermal source function
			S=CC%BB_l*(1d0-CC%albedo_l)*tau_dust
c	gas source function
			S=S+CC%line_emis*tau_gas

			tau_d=tau*p0%d(k)
			if(tau_d.gt.1d-4) then
				exptau=exp(-tau_d)
				flux=flux+S*(1d0-exptau)*fact/tau
			else
				exptau=1d0-tau_d
				flux=flux+S*p0%d(k)*fact
			endif

			fact=fact*exptau
			if(fact.lt.1d-6) exit
		endif
	enddo
	
	return
	end
	
