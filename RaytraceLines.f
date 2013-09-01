	subroutine RaytraceLines()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k,ilines,vmult,iv,nv,nl,imol
	real*8 lam,T,Planck,wl1,wl2,v,flux0,starttime,stoptime,tot,fact
	real*8,allocatable :: flux(:)
	type(Path),pointer :: PP
	type(Line) :: LL
	
	call output("==================================================================")
	call output("Preparing the profiles")

	do i=1,nImR
		do j=1,nImPhi
			allocate(P(i,j)%cont_contr(P(i,j)%n))
			allocate(P(i,j)%exptau_dust(P(i,j)%n))
		enddo
	enddo

	open(unit=20,file='out.dat',RECL=1000)

	nv=int(vmax*1.1/vresolution)+1
	nvprofile=int(vmax*vres_mult/vresolution)

	allocate(flux(-nv:nv))

	lam=lmin
	ilam1=1
	do while(lam.gt.lam_cont(ilam1+1).and.ilam1.lt.nlam)
		ilam1=ilam1+1
	enddo
	do i=0,nR
		do j=1,nTheta
			allocate(C(i,j)%profile(nmol,-nvprofile:nvprofile))
			allocate(C(i,j)%profile_nz(nmol,-nvprofile:nvprofile))
		enddo
	enddo

	do i=0,nR
		call tellertje(i+1,nR+1)
		do j=1,nTheta
			do imol=1,nmol
				do k=-nvprofile,nvprofile
					C(i,j)%profile(imol,k)=((real(k)*vresolution/vres_mult)/C(i,j)%line_width(imol))**2
					if(abs(C(i,j)%profile(imol,k)).lt.3d0) then
						C(i,j)%profile_nz(imol,k)=.true.
					else
						C(i,j)%profile_nz(imol,k)=.false.
					endif
				enddo
				C(i,j)%profile(imol,:)=clight*exp(-C(i,j)%profile(imol,:))
     &					/(C(i,j)%line_width(imol)*sqrt(pi))
			enddo
		enddo
	enddo

	nl=0
	
	call output("==================================================================")

	call cpu_time(starttime)

	call output("Tracing " // trim(int2string(nlines,'(i5)')) // " lines")

	do ilines=1,nlines
		call tellertje(ilines,nlines)
		flux=0d0

		LL = Lines(ilines)
		lam=clight*1d4/(LL%freq)
		if(lam.gt.lmin.and.lam.lt.lmax) then
		nl=nl+1

		ilam=ilam1
		do while(lam.gt.lam_cont(ilam+1).and.ilam.lt.nlam)
			ilam=ilam+1
		enddo
		if(ilam.gt.ilam1) then
			do k=ilam1+1,ilam
				flux0=0d0
				do i=1,nImR
					do j=1,nImPhi
						PP => P(i,j)
						flux0=flux0+PP%flux_cont(k)*PP%A
					enddo
				enddo
				flux0=flux0+path2star%flux_cont(k)*path2star%A
				write(20,*) lam_cont(k),flux0*1e23/(distance*parsec)**2
			enddo
			ilam1=ilam
		endif
			
		wl1=(lam_cont(ilam+1)-lam)/(lam_cont(ilam+1)-lam_cont(ilam))
		wl2=1d0-wl1

		do i=0,nR
			do j=1,nTheta
				C(i,j)%kext_l=wl1*C(i,j)%kext(ilam)+wl2*C(i,j)%kext(ilam+1)
				C(i,j)%albedo_l=wl1*C(i,j)%albedo(ilam)+wl2*C(i,j)%albedo(ilam+1)
				C(i,j)%BB_l=wl1*BB(ilam,C(i,j)%iT)+wl2*BB(ilam+1,C(i,j)%iT)
				C(i,j)%LRF_l=wl1*C(i,j)%LRF(ilam)+wl2*C(i,j)%LRF(ilam+1)

				fact=hplanck*C(i,j)%N(LL%imol)/(4d0*pi)
				C(i,j)%line_abs=fact*(C(i,j)%npop(LL%imol,LL%jlow)*LL%Blu-C(i,j)%npop(LL%imol,LL%jup)*LL%Bul)

				C(i,j)%line_emis=fact*C(i,j)%npop(LL%imol,LL%jup)*LL%Aul
			enddo
		enddo
		Fstar_l=wl1*Fstar(ilam)+wl2*Fstar(ilam+1)

		do i=1,nImR
			do j=1,nImPhi
				PP => P(i,j)
				if(PP%npopmax(LL%imol).gt.LL%jlow) then
					call ContContrPath(PP)
					do iv=-nv,nv
						vmult=1
						if(real(iv*vresolution).gt.PP%vmax(LL%imol)
     &					.or.real(iv*vresolution).lt.PP%vmin(LL%imol)) then
							flux0=wl1*PP%flux_cont(ilam)+wl2*PP%flux_cont(ilam+1)
						else
							call TraceFluxLines(PP,flux0,iv,LL%imol)
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
			if(real(iv*vresolution).gt.PP%vmax(LL%imol)
     &			.or.real(iv*vresolution).lt.PP%vmin(LL%imol)) then
				flux0=wl1*PP%flux_cont(ilam)+wl2*PP%flux_cont(ilam+1)
			else
				call Trace2StarLines(PP,flux0,iv,LL%imol)
			endif
			do vmult=-1,1,2
				flux(iv*vmult)=flux(iv*vmult)+flux0*PP%A/2d0
			enddo
		enddo

		do i=-nv,nv
			write(20,*) lam*sqrt((1d0+real(i)*vresolution/clight)/(1d0-real(i)*vresolution/clight)),
     &					flux(i)*1e23/(distance*parsec)**2,
     &					real(i)*vresolution/1d5,
     &					trim(Mol(LL%imol)%name),
     &					LL%jlow,LL%jup
		enddo
		
		endif
	enddo

	ilam=ilam+1
	do while(ilam.le.nlam)
		if(lam_cont(ilam).lt.lmax) then
			flux0=0d0
			do i=1,nImR
				do j=1,nImPhi
					PP => P(i,j)
					flux0=flux0+PP%flux_cont(ilam)*PP%A
				enddo
			enddo
			flux0=flux0+path2star%flux_cont(ilam)*path2star%A
			write(20,*) lam_cont(ilam),flux0*1e23/(distance*parsec)**2
			ilam=ilam+1
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
	
	
	subroutine ContContrPath(p0)
	use GlobalSetup
	use Constants
	integer i,j,k,vmult,iv,ii,nv,nn
	real*8 tau,exptau,flux,fact,profile,S,tau_gas,tau_dust,tau_d,tau_tot
	type(Path) p0
	type(Cell),pointer :: CC

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

			tau_d=tau*p0%d(k)
			if(tau_d.gt.1d-4) then
				p0%exptau_dust(k)=exp(-tau_d)
				p0%cont_contr(k)=S*(1d0-p0%exptau_dust(k))/tau
			else
				p0%exptau_dust(k)=1d0-tau_d
				p0%cont_contr(k)=S*p0%d(k)
			endif

			tau_tot=tau_tot+tau_d
			if(tau_tot.gt.tau_max) exit
		endif
	enddo
	
	return
	end
	


	subroutine TraceFluxLines(p0,flux,ii,imol)
	use GlobalSetup
	use Constants
	integer i,j,k,iv,ii,nv,ilines,imol
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

			tau=0d0
			S=0d0
			jj=int(p0%v(k)*vres_mult/vresolution-real(ii)*vres_mult)
			if(jj.lt.-nvprofile) jj=-nvprofile
			if(jj.gt.nvprofile) jj=nvprofile
			if(CC%profile_nz(imol,jj)) then
				profile=CC%profile(imol,jj)
				tau_gas=profile*CC%line_abs
c	gas source function
				S=S+CC%line_emis*profile
				tau=tau+tau_gas

				tau_dust=CC%kext_l
c	dust thermal source function
				S=S+CC%BB_l*(1d0-CC%albedo_l)*tau_dust
c	dust scattering source function
				S=S+CC%LRF_l*CC%albedo_l*tau_dust

				tau=tau+tau_dust

				tau_d=tau*p0%d(k)
				if(tau_d.gt.1d-4) then
					exptau=exp(-tau_d)
					flux=flux+S*(1d0-exptau)*fact/tau
				else
					exptau=1d0-tau_d
					flux=flux+S*p0%d(k)*fact
				endif
			else
				tau_d=CC%kext_l*p0%d(k)
				exptau=p0%exptau_dust(k)
				flux=flux+p0%cont_contr(k)*fact
			endif

			fact=fact*exptau
			tau_tot=tau_tot+tau_d
			if(tau_tot.gt.tau_max) exit
		endif
	enddo
	
	return
	end
	



	subroutine Trace2StarLines(p0,flux,ii,imol)
	use GlobalSetup
	use Constants
	integer i,j,k,iv,ii,nv,ilines,imol
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

		jj=int(p0%v(k)*vres_mult/vresolution-real(ii)*vres_mult)
		if(jj.lt.-nvprofile) jj=-nvprofile
		if(jj.gt.nvprofile) jj=nvprofile
		profile=CC%profile(imol,jj)
		tau=tau+profile*CC%line_abs
	enddo

	flux=Fstar_l*exp(-tau)
	
	return
	end
	
