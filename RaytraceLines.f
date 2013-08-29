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

	open(unit=20,file='out.dat',RECL=1000)

	nv=int(vmax*1.1/vresolution)+1
	nvprofile=int(vmax*vres_mult/vresolution)

	allocate(flux(-nv:nv))

	lam=lmin
	ilam=1
	do while(lam.gt.lam_cont(ilam+1).and.ilam.lt.nlam)
		ilam=ilam+1
	enddo
	do i=0,nR
		do j=1,nTheta
			allocate(C(i,j)%profile(-nvprofile:nvprofile))
			allocate(C(i,j)%profile_nz(-nvprofile:nvprofile))
		enddo
	enddo

	do i=0,nR
		do j=1,nTheta
			do k=-nvprofile,nvprofile
				C(i,j)%profile(k)=((real(k)*vresolution/vres_mult)/C(i,j)%line_width)**2
				if(abs(C(i,j)%profile(k)).lt.3d0) then
					C(i,j)%profile_nz(k)=.true.
				else
					C(i,j)%profile_nz(k)=.false.
				endif
			enddo
			C(i,j)%profile=clight*exp(-C(i,j)%profile)/(C(i,j)%line_width*sqrt(pi))
		enddo
	enddo

	nl=0
	
	print*,Mol%nlines
	
	do ilines=1,Mol%nlines
		call tellertje(ilines,Mol%nlines)
		flux=0d0

		LL = Mol%L(ilines)		
		lam=clight*1d4/(LL%freq)
		if(lam.gt.lmin.and.lam.lt.lmax) then
		nl=nl+1

		ilam=ilam1
		do while(lam.gt.lam_cont(ilam+1).and.ilam.lt.nlam)
			ilam=ilam+1
		enddo
		wl1=(lam_cont(ilam+1)-lam)/(lam_cont(ilam+1)-lam_cont(ilam))
		wl2=1d0-wl1

		do i=0,nR
			do j=1,nTheta
				C(i,j)%kext_l=wl1*C(i,j)%kext(ilam)+wl2*C(i,j)%kext(ilam+1)
				C(i,j)%albedo_l=wl1*C(i,j)%albedo(ilam)+wl2*C(i,j)%albedo(ilam+1)
				C(i,j)%BB_l=wl1*BB(ilam,C(i,j)%iT)+wl2*BB(ilam+1,C(i,j)%iT)
				C(i,j)%LRF_l=wl1*C(i,j)%LRF(ilam)+wl2*C(i,j)%LRF(ilam+1)

				fact=hplanck*C(i,j)%N/(4d0*pi)
				C(i,j)%line_abs=fact*(C(i,j)%npop(LL%jlow)*LL%Blu-C(i,j)%npop(LL%jup)*LL%Bul)

				C(i,j)%line_emis=fact*C(i,j)%npop(LL%jup)*LL%Aul
			enddo
		enddo
		Fstar_l=wl1*Fstar(ilam)+wl2*Fstar(ilam+1)

		do i=1,nImR
			do j=1,nImPhi
				PP => P(i,j)
				if(PP%npopmax.gt.LL%jlow) then
					do iv=-nv,nv
						vmult=1
						if(real(iv*vresolution).gt.PP%vmax.or.real(iv*vresolution).lt.PP%vmin) then
							flux0=wl1*PP%flux_cont(ilam)+wl2*PP%flux_cont(ilam+1)
						else
							call TraceFluxLines(PP,flux0,iv,vmult,1)
						endif
						do vmult=-1,1,2
							flux(iv*vmult)=flux(iv*vmult)+flux0*PP%A/2d0
						enddo
					enddo
				else
					flux0=wl1*PP%flux_cont(ilam)+wl2*PP%flux_cont(ilam+1)
					flux(-nv:nv)=flux(-nv:nv)+flux0*PP%A
				endif
			enddo
		enddo
		PP => path2star
		do iv=-nv,nv
			if(real(iv*vresolution).gt.PP%vmax.or.real(iv*vresolution).lt.PP%vmin) then
				flux0=wl1*PP%flux_cont(ilam)+wl2*PP%flux_cont(ilam+1)
			else
				call Trace2StarLines(PP,flux0,iv,1)
			endif
			do vmult=-1,1,2
				flux(iv*vmult)=flux(iv*vmult)+flux0*PP%A/2d0
			enddo
		enddo

		do i=-nv,nv
			write(20,*) lam*sqrt((1d0+real(i)*vresolution/clight)/(1d0-real(i)*vresolution/clight)),
     &					flux(i)*1e23/(distance*parsec)**2,
     &					real(i)*vresolution/1d5,
     &					LL%jlow,LL%jup
		enddo
		
		endif

	enddo

	do ilam=1,nlam
		if(lam_cont(ilam-1).lt.lmax.and.lam_cont(ilam+1).gt.lmin) then
			flux0=0d0
			do i=1,nImR
				do j=1,nImPhi
					PP => P(i,j)
					flux0=flux0+PP%flux_cont(ilam)*PP%A
				enddo
			enddo
			flux0=flux0+path2star%flux_cont(ilam)*path2star%A
			write(20,*) lam_cont(ilam),flux0*1e23/(distance*parsec)**2
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
	
	
	subroutine TraceFluxLines(p0,flux,ii,vmult,nn)
	use GlobalSetup
	use Constants
	integer i,j,k,vmult,iv,ii,nv,nn
	real*8 tau,exptau,flux,fact,profile,S,tau_gas,tau_dust,tau_d,tau_tot
	type(Path) p0
	type(Cell),pointer :: CC

	fact=1d0
	flux=0d0
	tau_tot=0d0

	do k=1,p0%n
		i=p0%i(k)
		if(i.ne.0) then
			j=p0%j(k)
			CC => C(i,j)
			tau_dust=CC%kext_l
c	dust thermal source function
			S=CC%BB_l*(1d0-CC%albedo_l)*tau_dust
c	dust scattering source function
			S=S+CC%LRF_l*CC%albedo_l*tau_dust

			tau=tau_dust
			do il=1,nn
				jj=int(p0%v(k)*vres_mult/vresolution-real(ii)*vres_mult)
				if(jj.lt.-nvprofile) jj=-nvprofile
				if(jj.gt.nvprofile) jj=nvprofile
				if(CC%profile_nz(jj)) then
					profile=CC%profile(jj)
					tau_gas=profile*CC%line_abs
c	gas source function
					S=S+CC%line_emis*profile
					tau=tau+tau_gas
				endif
			enddo

			tau_d=tau*p0%d(k)
			if(tau_d.gt.1d-4) then
				exptau=exp(-tau_d)
				flux=flux+S*(1d0-exptau)*fact/tau
			else
				exptau=1d0-tau_d
				flux=flux+S*p0%d(k)*fact
			endif

			fact=fact*exptau
			tau_tot=tau_tot+tau_d
			if(tau_tot.gt.tau_max) exit
		endif
	enddo
	
	return
	end
	



	subroutine Trace2StarLines(p0,flux,ii,nn)
	use GlobalSetup
	use Constants
	integer i,j,k,iv,ii,nv,nn
	real*8 tau,flux,profile
	type(Path) p0
	type(Cell),pointer :: CC

	tau=0d0

	do k=1,p0%n
		i=p0%i(k)
		if(i.eq.0) exit
		j=p0%j(k)
		CC => C(i,j)
		tau=tau+CC%kext_l

		do il=1,nn
			jj=int(p0%v(k)*vres_mult/vresolution-real(ii)*vres_mult)
			if(jj.lt.-nvprofile) jj=-nvprofile
			if(jj.gt.nvprofile) jj=nvprofile
			profile=CC%profile(jj)
			tau=tau+profile*CC%line_abs
		enddo
	enddo

	flux=Fstar_l*exp(-tau)
	
	return
	end
	
