	subroutine RaytraceLines()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k,iblends,vmult,iv,nv,nl,imol,maxblend,ilines,nvmax,nb,ib0,nb0,ib
	integer,allocatable :: imol_blend(:),count(:)
	real*8,allocatable :: v_blend(:)
	real*8 lam,T,Planck,wl1,wl2,v,flux0,starttime,stoptime,tot,fact,lcmin
	real*8,allocatable :: flux(:),flux_cont(:)
	type(Path),pointer :: PP
	type(Line) :: LL
	type(Blend),pointer :: Bl
	logical gas
	real*8 flux1,flux2,flux3,fc,f
	real*8 wl11,wl21,wl12,wl22,wl13,wl23,flux_l1,flux_l2,flux_c
	character*1000 comment
	
	idum=-42
	
	call output("==================================================================")
	call output("Preparing the profiles")

	do i=1,nImR
		do j=1,nImPhi(i)
			allocate(P(i,j)%cont_contr(P(i,j)%n))
			allocate(P(i,j)%exptau_dust(P(i,j)%n))
			allocate(P(i,j)%S_dust(P(i,j)%n))
		enddo
	enddo

	open(unit=20,file='specFLiTs.out',RECL=1500)

	nv=int(vmax*1.1/vresolution)+1
	nvprofile=int(vmax*vres_mult/vresolution)

	call DetermineBlends(nv,maxblend)

	allocate(flux(-nv:nv+nv*2*(maxblend)))
	allocate(flux_cont(-nv:nv+nv*2*(maxblend)))
	allocate(imol_blend(maxblend))
	allocate(v_blend(maxblend))
	allocate(count(nmol))

	lam=lmin
	ilam1=1
	do while(lam.gt.lam_cont(ilam1+1).and.ilam1.lt.nlam)
		ilam1=ilam1+1
	enddo
	allocate(profile(-nvprofile:nvprofile))
	allocate(profile_nz(-nvprofile:nvprofile))
	do i=0,nR
		do j=0,nTheta
			allocate(C(i,j)%line_abs(maxblend))
			allocate(C(i,j)%line_emis(maxblend))
		enddo
	enddo

	do k=-nvprofile,nvprofile
		profile(k)=((real(k)*vresolution/vres_mult)/1d5)**2
		if(abs(profile(k)).lt.3d0) then
			profile_nz(k)=.true.
		else
			profile_nz(k)=.false.
		endif
	enddo
	profile(:)=exp(-profile(:))

	nl=0
	
	call output("==================================================================")

	call cpu_time(starttime)

	call output("Tracing " // trim(int2string(nlines,'(i5)')) // " lines")

	lcmin=lmin

	Bl => Blends
	do iblends=1,nblends
		flux=0d0
		
		nb=Bl%n
		do i=1,nb
			imol_blend(i)=Bl%L(i)%imol
			v_blend(i)=Bl%v(i)
		enddo

		LL = Bl%L(1)
		lam=LL%lam

		if(lam.gt.lmin.and.lam.lt.lmax) then
		nl=nl+nb

		ilam=ilam1
		do while(lam.gt.lam_cont(ilam+1).and.ilam.lt.nlam)
			ilam=ilam+1
			if(ilam.eq.nlam) exit
		enddo

		if(ilam.gt.ilam1) then
			do k=ilam1+1,ilam
				if(lam_cont(k).gt.lcmin.and.lam_cont(k).lt.Bl%lmin) then
					flux0=0d0
					do i=1,nImR
						do j=1,nImPhi(i)
							PP => P(i,j)
							flux0=flux0+PP%flux_cont(k)*PP%A
						enddo
					enddo
					flux0=flux0+path2star%flux_cont(k)*path2star%A
					write(20,*) lam_cont(k),flux0*1e23/(distance*parsec)**2
				endif
			enddo
			ilam1=ilam
		endif
		lcmin=Bl%lmax

c		wl1=(lam_cont(ilam+1)-lam)/(lam_cont(ilam+1)-lam_cont(ilam))
c		wl2=1d0-wl1
		wl1=log10(lam_cont(ilam+1)/lam)/log10(lam_cont(ilam+1)/lam_cont(ilam))
		wl2=1d0-wl1

		do i=0,nR
			do j=1,nTheta
				C(i,j)%kext_l=10d0**(wl1*log10(C(i,j)%kext(ilam))+wl2*log10(C(i,j)%kext(ilam+1)))
				C(i,j)%albedo_l=10d0**(wl1*log10(C(i,j)%albedo(ilam))+wl2*log10(C(i,j)%albedo(ilam+1)))
				C(i,j)%BB_l=10d0**(wl1*log10(BB(ilam,C(i,j)%iT))+wl2*log10(BB(ilam+1,C(i,j)%iT)))
				C(i,j)%LRF_l=10d0**(wl1*log10(C(i,j)%LRF(ilam))+wl2*log10(C(i,j)%LRF(ilam+1)))

				do ilines=1,Bl%n
					LL = Bl%L(ilines)
					fact=clight*hplanck*C(i,j)%N(LL%imol)/(4d0*pi*C(i,j)%line_width(LL%imol)*sqrt(pi))
					C(i,j)%line_abs(ilines)=fact*(C(i,j)%npop(LL%imol,LL%jlow)*LL%Blu-C(i,j)%npop(LL%imol,LL%jup)*LL%Bul)
					C(i,j)%line_emis(ilines)=fact*C(i,j)%npop(LL%imol,LL%jup)*LL%Aul
				enddo
			enddo
		enddo
		Fstar_l=10d0**(wl1*log10(Fstar(ilam))+wl2*log10(Fstar(ilam+1)))

		nvmax=nv+int(v_blend(nb)/vresolution)

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,PP,iv,vmult,ib0,nb0,gas,ib,imol,flux0)
!$OMP& SHARED(P,nImPhi,nImR,nb,LL,nv,nvmax,Bl,vresolution,wl1,wl2,imol_blend,v_blend,flux,ilam)
!$OMP DO
		do i=1,nImR
			do j=1,nImPhi(i)
				PP => P(i,j)
				call ContContrPath(PP,flux_c)
				if(nb.gt.1.or.PP%npopmax(LL%imol).gt.LL%jlow) then
					do iv=-nv,nvmax
						vmult=1
						if(nb.gt.1) then
							ib0=Bl%ib0(iv)
							nb0=Bl%nb0(iv)
							vmult=-1
							gas=.false.
							do ib=ib0,ib0+nb0-1
								imol=Bl%L(ib)%imol
								if((real(iv*vresolution)-Bl%v(ib)).lt.-PP%vmin(imol)
     &								.and.(real(iv*vresolution)-Bl%v(ib)).gt.-PP%vmax(imol)) then
	   								gas=.true.
	   								exit
	   							endif
							enddo
							if(gas) then
								call TraceFluxLines(PP,flux0,iv,vmult,imol_blend(ib0),v_blend(ib0),nb0,ib0)
							else
								flux0=flux_c
							endif
							flux(iv)=flux(iv)+flux0*PP%A/2d0

							vmult=1
							gas=.false.
							do ib=ib0,ib0+nb0-1
								imol=Bl%L(ib)%imol
								if((real(iv*vresolution)-Bl%v(ib)).lt.PP%vmax(imol)
     &								.and.(real(iv*vresolution)-Bl%v(ib)).gt.PP%vmin(imol)) then
	   								gas=.true.
	   								exit
	   							endif
							enddo
							if(gas) then
								call TraceFluxLines(PP,flux0,iv,vmult,imol_blend(ib0),v_blend(ib0),nb0,ib0)
							else
								flux0=flux_c
							endif
							flux(iv)=flux(iv)+flux0*PP%A/2d0
						else if(real(iv*vresolution).gt.PP%vmax(LL%imol)
     &					.or.real(iv*vresolution).lt.PP%vmin(LL%imol)) then
							flux0=flux_c
							do vmult=-1,1,2
								flux(iv*vmult)=flux(iv*vmult)+flux0*PP%A/2d0
							enddo
						else
							call TraceFluxLines(PP,flux0,iv,vmult,imol_blend,v_blend,nb,1)
							do vmult=-1,1,2
								flux(iv*vmult)=flux(iv*vmult)+flux0*PP%A/2d0
							enddo
						endif
					enddo
				else
					flux0=flux_c
					flux(-nv:nvmax)=flux(-nv:nvmax)+flux0*PP%A
				endif
			enddo
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		PP => path2star
		do iv=-nv,nvmax
			if(nb.gt.1) then
				call Trace2StarLines(PP,flux0,iv,imol_blend,v_blend,nb)
				flux(iv)=flux(iv)+flux0*PP%A
			else
				call Trace2StarLines(PP,flux0,iv,imol_blend,v_blend,nb)
				flux(iv)=flux(iv)+flux0*PP%A
			endif
		enddo

		flux1=0d0
		flux2=0d0
		flux3=0d0

c		f=sqrt((1d0+real(-nv)*vresolution/clight)/(1d0-real(-nv)*vresolution/clight))
c		wl11=(lam_cont(ilam+1)-(lam*f))/(lam_cont(ilam+1)-lam_cont(ilam))
c		wl21=1d0-wl11
c		f=1d0
c		wl12=(lam_cont(ilam+1)-(lam*f))/(lam_cont(ilam+1)-lam_cont(ilam))
c		wl22=1d0-wl12
c		f=sqrt((1d0+real(nvmax)*vresolution/clight)/(1d0-real(nvmax)*vresolution/clight))
c		wl13=(lam_cont(ilam+1)-(lam*f))/(lam_cont(ilam+1)-lam_cont(ilam))
c		wl23=1d0-wl13

		f=sqrt((1d0+real(-nv)*vresolution/clight)/(1d0-real(-nv)*vresolution/clight))
		wl11=log10(lam_cont(ilam+1)/(lam*f))/log10(lam_cont(ilam+1)/lam_cont(ilam))
		wl21=1d0-wl11
		f=1d0
		wl12=log10(lam_cont(ilam+1)/(lam*f))/log10(lam_cont(ilam+1)/lam_cont(ilam))
		wl22=1d0-wl12
		f=sqrt((1d0+real(nvmax)*vresolution/clight)/(1d0-real(nvmax)*vresolution/clight))
		wl13=log10(lam_cont(ilam+1)/(lam*f))/log10(lam_cont(ilam+1)/lam_cont(ilam))
		wl23=1d0-wl13

		flux_l1=0d0
		flux_l2=0d0
		do i=1,nImR
			do j=1,nImPhi(i)
				PP => P(i,j)

				flux_l1=flux_l1+PP%flux_cont(ilam)*PP%A
				flux_l2=flux_l2+PP%flux_cont(ilam+1)*PP%A
			enddo
		enddo		
		PP => path2star

		flux_l1=flux_l1+PP%flux_cont(ilam)*PP%A
		flux_l2=flux_l2+PP%flux_cont(ilam+1)*PP%A

		flux1=10d0**(wl11*log10(flux_l1)+wl21*log10(flux_l2))
		flux2=10d0**(wl12*log10(flux_l1)+wl22*log10(flux_l2))
		flux3=10d0**(wl13*log10(flux_l1)+wl23*log10(flux_l2))

		do i=-nv,nvmax
			fc=-flux2+flux1+(flux3-flux1)*real(i+nv)/real(nvmax+nv)
			flux(i)=flux(i)+fc
			flux_cont(i)=flux1+(flux3-flux1)*real(i+nv)/real(nvmax+nv)
		enddo

		if(Bl%n.eq.1) then
			comment = trim(Mol(Bl%L(1)%imol)%name) // "  up: " // trim(int2string(Bl%L(1)%jup,'(i5)'))
     &											  // "  low: " // trim(int2string(Bl%L(1)%jlow,'(i5)'))
		else
			count=0
			do i=1,Bl%n
				count(Bl%L(i)%imol)=count(Bl%L(i)%imol)+1
			enddo
			comment = "blend of "
			do i=1,nmol
				if(count(i).gt.0) then
					comment=trim(comment) // trim(int2string(count(i),'(i3)')) // " " // trim(Mol(i)%name)
				endif
			enddo
		endif

		do i=-nv,nvmax
			write(20,*) lam*sqrt((1d0+real(i)*vresolution/clight)/(1d0-real(i)*vresolution/clight)),
     &					flux(i)*1e23/(distance*parsec)**2,
     &					real(i)*vresolution/1d5,
     &					flux_cont(i)*1e23/(distance*parsec)**2,
     &					trim(comment)
		enddo
		
		endif

		if(iblends.lt.nblends) Bl => Bl%next
				
c		call cpu_time(stoptime)
c		call output("Time used per line:     "//trim(dbl2string((stoptime-starttime)/real(nl),'(f8.2)'))//" s")
		call tellertje_time(iblends,nblends,nl,nlines,starttime)
	enddo

	ilam=ilam+1
	do while(ilam.le.nlam)
		if(lam_cont(ilam).lt.lmax) then
			flux0=0d0
			do i=1,nImR
				do j=1,nImPhi(i)
					PP => P(i,j)
					flux0=flux0+PP%flux_cont(ilam)*PP%A
				enddo
			enddo
			flux0=flux0+path2star%flux_cont(ilam)*path2star%A
			write(20,*) lam_cont(ilam),flux0*1e23/(distance*parsec)**2
		endif
		ilam=ilam+1
	enddo

	close(unit=20)
	
	call cpu_time(stoptime)

	call output("Time used for the lines:"//trim(dbl2string(stoptime-starttime,'(f8.2)'))
     &			//" s")
	call output("Time used per line:     "//trim(dbl2string((stoptime-starttime)/real(nl),'(f8.2)'))
     &			//" s")


	return
	end
	
	
	subroutine ContContrPath(p0,flux)
	use GlobalSetup
	use Constants
	integer i,j,k,vmult,iv,ii,nv,nn
	real*8 tau,exptau,flux,fact,prof,S,tau_gas,tau_dust,tau_d,tau_tot
	type(Path) p0
	type(Cell),pointer :: CC

	fact=1d0
	flux=0d0
	tau_tot=0d0

	do k=1,p0%n
		i=p0%i(k)
		j=p0%j(k)
		if(i.gt.0.and.i.lt.nR.and.j.gt.0) then
			CC => C(i,j)
			tau_dust=CC%kext_l
c	dust thermal source function
			S=CC%BB_l*(1d0-CC%albedo_l)*tau_dust
c	dust scattering source function
			S=S+CC%LRF_l*CC%albedo_l*tau_dust

			p0%S_dust(k)=S

			tau=tau_dust

			tau_d=tau*p0%d(k)
			if(tau_d.gt.1d-4) then
				p0%exptau_dust(k)=exp(-tau_d)
				p0%cont_contr(k)=S*(1d0-p0%exptau_dust(k))/tau
			else
				p0%exptau_dust(k)=1d0-tau_d
				p0%cont_contr(k)=S*p0%d(k)
			endif
			flux=flux+p0%cont_contr(k)*fact

			fact=fact*p0%exptau_dust(k)
			tau_tot=tau_tot+tau_d
			if(tau_tot.gt.tau_max) return
		endif
	enddo
	
	return
	end
	


	subroutine TraceFluxLines(p0,flux,ii,vmult,imol_blend,v_blend,nb,ib0)
	use GlobalSetup
	use Constants
	integer i,j,k,iv,ii,nv,ilines,imol,vmult,nb,imol_blend(nb)
	real*8 tau,exptau,flux,fact,prof,S,tau_gas,tau_dust,tau_d,tau_tot,v_blend(nb)
	type(Path) p0
	type(Cell),pointer :: CC
	type(Blend) Bl
	logical gas
	
	fact=1d0
	flux=0d0
	tau_tot=0d0

	v=real(ii)+ran1(idum)-0.5d0
	do k=1,p0%n
		i=p0%i(k)
		j=p0%j(k)
		if(i.gt.0.and.i.lt.nR.and.j.gt.0) then
			CC => C(i,j)

			tau=0d0
			S=0d0
			gas=.false.

			do ib=1,nb
				imol=imol_blend(ib)
				jj=int(((real(vmult)*p0%v(k)+v_blend(ib))*vres_mult/vresolution-v*vres_mult)
     &					*1d5/CC%line_width(imol))
				if(jj.lt.-nvprofile) jj=-nvprofile
				if(jj.gt.nvprofile) jj=nvprofile
				if(profile_nz(jj)) then
					prof=profile(jj)
					tau_gas=prof*CC%line_abs(ib+ib0-1)
c	gas source function
					S=S+prof*CC%line_emis(ib+ib0-1)
					tau=tau+tau_gas
					gas=.true.
				endif
			enddo

			if(gas) then
				tau_dust=CC%kext_l
				S=S+p0%S_dust(k)

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

			if(tau_tot.gt.tau_max) return
		endif
	enddo
	
	return
	end
	



	subroutine Trace2StarLines(p0,flux,ii,imol_blend,v_blend,nb)
	use GlobalSetup
	use Constants
	integer i,j,k,iv,ii,nv,ilines,imol,nb,imol_blend(nb)
	real*8 tau,flux,prof,v_blend(nb)
	type(Path) p0
	type(Cell),pointer :: CC
	type(Blend) Bl

	tau=0d0

	do k=1,p0%n
		i=p0%i(k)
		if(i.eq.0) exit
		j=p0%j(k)
		if(i.ne.0.and.i.ne.nR.and.j.ne.0) then
			CC => C(i,j)
			tau=tau+CC%kext_l

			do ib=1,nb
				imol=imol_blend(ib)
				jj=int(((p0%v(k)+v_blend(ib))*vres_mult/vresolution-real(ii)*vres_mult)*1d5/CC%line_width(imol))
				if(jj.lt.-nvprofile) jj=-nvprofile
				if(jj.gt.nvprofile) jj=nvprofile
				prof=profile(jj)
				tau=tau+prof*CC%line_abs(ib)
			enddo
		endif
	enddo

	flux=Fstar_l*exp(-tau)
	
	return
	end
	

	subroutine DetermineBlends(nv,maxblend)
	use GlobalSetup
	use Constants
	integer ilines,ilines0,i,maxblend,nv,nvmax,iv
	real*8 maxvshift,maxmult,v
	type(Blend),pointer :: Bl

	maxvshift=2d0*real(nv)*vresolution
	maxmult=sqrt((1d0+maxvshift/clight)/(1d0-maxvshift/clight))
	
	Bl => Blends

	ilines=2
	ilines0=1
	nblends=0
	maxblend=0
	do while(ilines0.le.nlines)
		Bl%n=1
		if(doblend) then
			do while((Lines(ilines)%lam/Lines(ilines-1)%lam).lt.maxmult.and.ilines.le.nlines)
				Bl%n=Bl%n+1
				ilines=ilines+1
				if(ilines.gt.nlines) exit
			enddo
		endif
		allocate(Bl%L(Bl%n))
		allocate(Bl%v(Bl%n))
		if(Bl%n.gt.maxblend) maxblend=Bl%n
		do i=1,Bl%n
			Bl%L(i) = Lines(ilines0+i-1)
			if(i.eq.1) then
				Bl%v(i)=0d0
			else
				f=(Bl%L(i)%lam/Bl%L(1)%lam)**2
				Bl%v(i)=clight*(f-1d0)/(f+1d0)
			endif
		enddo
		nvmax=nv+int(Bl%v(Bl%n)/vresolution)
		allocate(Bl%ib0(-nv:nvmax))
		allocate(Bl%nb0(-nv:nvmax))
		do iv=-nv,nvmax
			do i=1,Bl%n
				if(real(iv-nv)*vresolution.le.Bl%v(i)) exit
			enddo
			Bl%ib0(iv)=i
			do i=Bl%n,1,-1
				if(real(iv+nv)*vresolution.ge.Bl%v(i)) exit
			enddo
			Bl%nb0(iv)=i-Bl%ib0(iv)+1
		enddo

		v=-real(nv)*vresolution
		Bl%lmin=Bl%L(1)%lam*sqrt((1d0+v/clight)/(1d0-v/clight))
		v=real(nvmax)*vresolution
		Bl%lmax=Bl%L(1)%lam*sqrt((1d0+v/clight)/(1d0-v/clight))

		ilines0=ilines
		ilines=ilines+1
		allocate(Bl%next)
		Bl => Bl%next
		nblends=nblends+1
	enddo
	
	return
	end
		
	
	
	
	
	
	
	
	
	
	
	
	
	
