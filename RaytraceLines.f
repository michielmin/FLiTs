	subroutine RaytraceLines()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k,iblends,vmult,iv,nv,nl,imol,maxblend,ilines,nvmax,nvmin
	integer nb,ib0,nb0,ib,nltot
	integer*4 counts,count_rate,count_max
	integer,external :: OMP_GET_THREAD_NUM
	integer,allocatable :: imol_blend(:),count(:)
	real*8,allocatable :: v_blend(:),flux4(:)
	real*8 lam,T,Planck,wl1,wl2,v,flux0,starttime,stoptime,tot,fact,lcmin
	real*8 lmin_next,lam_velo,lam_velo_imcube
	real*8,allocatable :: flux(:),flux_cont(:)
	type(Path),pointer :: PP
	type(Line) :: LL
	type(Blend),pointer :: Bl
	logical gas,doit
	logical,allocatable :: doit_ib(:),doit_ib0(:)
	real*8 flux1,flux2,flux3,fc,f,dnu,lam_w,lam_w_min,lam_w_max
	real*8 wl11,wl21,wl12,wl22,wl13,wl23,flux_l1,flux_l2,flux_c
	character*1000 comment
	character*500 imcubename
	character*40 callerstr
	real*8 :: delpix
	interface
	  subroutine output(string)
	  IMPLICIT NONE
	  character string*(*)
	  end
	  character*20 function int2string(i,form)
	  IMPLICIT NONE
	  integer i
	  character,intent(in),optional :: form*(*)
	  end
	  character*20 function dbl2string(x,form)
	  IMPLICIT NONE
	  real*8 x
	  character,intent(in),optional :: form*(*)
	  end
	end interface
		
	call output("==================================================================")
	call output("Preparing the profiles")

	do i=1,ngrids
		do j=1,npoints(i)
			allocate(P(i,j)%cont_contr(P(i,j)%n))
			allocate(P(i,j)%exptau_dust(P(i,j)%n))
			allocate(P(i,j)%S_dust(P(i,j)%n))
		enddo
	enddo

	open(unit=20,file='specFLiTs.out',RECL=1500)
	open(unit=21,file='lineFlux_FLiTs.out',RECL=1500)
	write(21,'("# lam[mu]   Fline[W/m^2]    lmin[mu]    lmax[mu]    comment")')

	nv=int(vmax*1.1/vresolution)+1
	nvprofile=int(vmax*vres_mult/vresolution)

	! this should always give nvmin=-nv and nmvmax=nv
	call DetermineBlends(nv,maxblend,nvmin,nvmax)

	allocate(flux(nvmin:nvmax))
	allocate(flux_cont(nvmin:nvmax))
	allocate(imol_blend(maxblend))
	allocate(v_blend(maxblend))
	allocate(count(nmol))

	if(imagecube) then
		allocate(imcube(npix,npix,nvmin:nvmax))
		allocate(im_coord(npix))
		delpix=2.d0*Rout*AU/real(npix)	
		! create the center coordinates for the pixels		
		do i=1,npix ! center coordinates of pixels in cm
			im_coord(i)=delpix/2d0+delpix*(i-1)-Rout*AU
		enddo	
		!write(*,*) im_coord/AU
	endif

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
		if(abs(profile(k)).lt.10d0) then
			profile_nz(k)=.true.
		else
			profile_nz(k)=.false.
		endif
	enddo
	profile(:)=exp(-profile(:))

	nl=0
	
	call output("==================================================================")

	!call cpu_time(starttime)
	call SYSTEM_CLOCK(counts, count_rate, count_max)
        starttime = DBLE(counts)/DBLE(count_rate)

	call output("Tracing " // trim(int2string(nlines,'(i5)')) // " lines")

	lcmin=lmin
	lmin_next=0d0

	nltot=0
	Bl => Blends
	do iblends=1,nblends
		nltot=nltot+Bl%n
		if(iblends.lt.nblends) Bl => Bl%next
	enddo
	print*,nltot,nlines,nblends

	Bl => Blends
	do iblends=1,nblends
		if(imagecube) then
			imcube=0d0 ! FIXME: in the end I want only one cube, not for each blend, I think that shoulbe be moved outside of the blends loop (like flux)
		endif
		flux=0d0
		
		nb=Bl%n
		if(iblends.eq.1) call output("First line blend (" // trim(int2string(nb,'(i4)')) // " lines)")

		do ib=1,nb
			imol_blend(ib)=Bl%L(ib)%imol
			v_blend(ib)=Bl%v(ib)
		enddo

		LL = Bl%L(1)
		lam=Bl%lam

		if(lam.gt.lmin.and.lam.lt.lmax) then
		nl=nl+nb

		ilam=ilam1
		do while(lam.gt.lam_cont(ilam+1).and.ilam.lt.nlam)
			ilam=ilam+1
			if(ilam.eq.nlam) exit
		enddo

		! CHR: What is this doing ?
		if(ilam.gt.ilam1) then
			do k=ilam1+1,ilam
				if(lam_cont(k).gt.lcmin.and.lam_cont(k).lt.Bl%lmin) then
					flux0=0d0
					do i=1,ngrids
					do j=1,npoints(i)
						PP => P(i,j)
						flux0=flux0+PP%flux_cont(k)*PP%A/real(ngrids)
					enddo
					enddo
					flux0=flux0+path2star%flux_cont(k)*path2star%A
					write(20,*) lam_cont(k),flux0*1e23/(distance*parsec)**2
				endif
			enddo
			ilam1=ilam
		endif
		lcmin=Bl%lmax

		call InterpolateLam(lam,ilam)				
		do i=0,nR
			do j=1,nTheta
				do ilines=1,Bl%n
					LL = Bl%L(ilines)
					fact=clight*hplanck*C(i,j)%N(LL%imol)/(4d0*pi*C(i,j)%line_width(LL%imol)*sqrt(pi))
					C(i,j)%line_abs(ilines)=fact*(C(i,j)%npop(LL%imol)%N(LL%jlow)*LL%Blu-C(i,j)%npop(LL%imol)%N(LL%jup)*LL%Bul)
					C(i,j)%line_emis(ilines)=fact*C(i,j)%npop(LL%imol)%N(LL%jup)*LL%Aul
				enddo
			enddo
		enddo

		flux2=0d0

		i=ran1(idum)*real(ngrids)+1


		if (imagecube) call map_pixels_to_path(i,npoints(i))

		! FIXME: check parallel implementation for imagecube
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(j,PP,iv,vmult,ib0,nb0,gas,ib,imol,flux0,flux4,doit,flux_c,doit_ib,doit_ib0,lam_velo)
!$OMP& SHARED(i,P,ngrids,npoints,nb,LL,nv,Bl,vresolution,wl1,wl2,imol_blend,v_blend,flux,ilam)
!$OMP& SHARED(iblends,flux2,maxblend,lam,lmin_next,imagecube,nvmin,nvmax)
		allocate(doit_ib0(maxblend))
		allocate(doit_ib(maxblend))
		allocate(flux4(nvmin:nvmax))

#ifdef USE_OPENMP
		idum = -42-OMP_GET_THREAD_NUM()     ! setting the seed 
#endif

		flux4(:)=0d0
!$OMP DO SCHEDULE(dynamic,1)
		do j=1,npoints(i)
			if(iblends.eq.1) call tellertje(j,npoints(i))
			PP => P(i,j)

			!write(88,*) i,j,PP%x/AU,PP%y/AU

			call ContContrPath(PP,flux_c)
!$OMP ATOMIC
			flux2=flux2+flux_c*PP%A
			doit=.false.
			doit_ib0=.false.
			do ib=1,nb
				if(PP%npopmax(Bl%L(ib)%imol).gt.Bl%L(ib)%jlow) then
					doit=.true.
					doit_ib0(ib)=.true.
				endif
			enddo
			if(doit) then
				!write(*,*) "doit",i,j
				!write(*,*) Bl%nvmin,Bl%nvmax
				do iv=Bl%nvmin,Bl%nvmax
					lam_velo=lam*sqrt((1d0+real(iv)*vresolution/clight)/(1d0-real(iv)*vresolution/clight))
					if(lam_velo.gt.lmin_next) then
						vmult=1
						if(nb.gt.1) then
							ib0=Bl%ib0(iv)
							nb0=Bl%nb0(iv)
							vmult=-1
							doit_ib(1:nb)=doit_ib0(1:nb)
							gas=.false.
							do ib=ib0,ib0+nb0-1
								imol=Bl%L(ib)%imol
								if(((real(iv)-0.5d0)*vresolution-Bl%v(ib)).lt.-PP%vmin(imol)
     &                        .and.((real(iv)+0.5d0)*vresolution-Bl%v(ib)).gt.-PP%vmax(imol)) then
									gas=.true.
									exit
								endif
							enddo
							if(gas) then
								call TraceFluxLines(PP,flux0,iv,vmult,imol_blend(ib0),v_blend(ib0),doit_ib(ib0),nb0,ib0)
							else
								flux0=flux_c
							endif
							flux4(iv)=flux4(iv)+flux0*PP%A/2d0
							if(imagecube) call AddImage(iv,i,j,flux0*PP%A/2d0,PP,vmult,"call nb 1")
							vmult=1
							doit_ib(1:nb)=doit_ib0(1:nb)
							gas=.false.
							do ib=ib0,ib0+nb0-1
								imol=Bl%L(ib)%imol
								if(((real(iv)-0.5d0)*vresolution-Bl%v(ib)).lt.PP%vmax(imol)
     &                        .and.((real(iv)+0.5d0)*vresolution-Bl%v(ib)).gt.PP%vmin(imol)) then
									gas=.true.
									exit
								endif
							enddo
							if(gas) then
								call TraceFluxLines(PP,flux0,iv,vmult,imol_blend(ib0),v_blend(ib0),doit_ib(ib0),nb0,ib0)
							else
								flux0=flux_c
							endif
							flux4(iv)=flux4(iv)+flux0*PP%A/2d0
							! FIXME: I think here  times vmult is required, but then it should also be in the line above ? or maybe not
							! or maybe in that case vmult also needs to be ignored for the imagecube
							if(imagecube) call AddImage(iv,i,j,flux0*PP%A/2d0,PP,vmult,"call nb 2")

						else if(((real(iv)+0.5d0)*vresolution.gt.PP%vmax(LL%imol).and.
     &                           (real(iv)-0.5d0)*vresolution.gt.PP%vmax(LL%imol))
     &                      .or.((real(iv)+0.5d0)*vresolution.lt.PP%vmin(LL%imol).and.
     &                           (real(iv)-0.5d0)*vresolution.lt.PP%vmin(LL%imol))) then
							flux0=flux_c
							do vmult=-1,1,2
								flux4(iv*vmult)=flux4(iv*vmult)+flux0*PP%A/2d0								
								! Seems to be required here and is done twice
								! FIXME: callerstring is only for debugging, remove it 
								write(callerstr,"(A,' ' ,i5,' ',i2)") "call nb else if", iv,vmult	
								if(imagecube) then
									call AddImage(iv*vmult,i,j,flux0*PP%A/2d0,PP,vmult,callerstr)
								endif
							enddo							
						else
							doit_ib=.true.
							call TraceFluxLines(PP,flux0,iv,vmult,imol_blend,v_blend,doit_ib,nb,1)
							do vmult=-1,1,2
								flux4(iv*vmult)=flux4(iv*vmult)+flux0*PP%A/2d0								
								write(callerstr,"(A,' ' ,i5,' ',i2)") "call nb else", iv,vmult																
								if(imagecube) then 
									call AddImage(iv*vmult,i,j,flux0*PP%A/2d0,PP,vmult,callerstr)
								endif
							enddo							
						endif
					endif
				enddo
			else
				flux0=flux_c
				flux4(Bl%nvmin:Bl%nvmax)=flux4(Bl%nvmin:Bl%nvmax)+flux0*PP%A
				if(imagecube) then
					do iv=Bl%nvmin,Bl%nvmax
						call AddImage(iv,i,j,flux0*PP%A,PP,1,"call not doit")
					enddo
				endif
			endif
		enddo		
!$OMP END DO		
!$OMP CRITICAL
		flux(:)=flux(:)+flux4(:)
!$OMP END CRITICAL
		deallocate(doit_ib0)
		deallocate(doit_ib)
		deallocate(flux4)
!$OMP FLUSH
!$OMP END PARALLEL
		PP => path2star
		! FIXME: whatever is done here is not considered for the Images
		! I guess it is for the star ... need to check what are the coordinates, then one might just need to calls AddImage again
		if (.not.allocated(PP%im_ixy)) allocate(PP%im_ixy(2,1))
		do iv=Bl%nvmin,Bl%nvmax
			if(nb.gt.1) then  ! FIXME: this if seems to be unnecessary
				call Trace2StarLines(PP,flux0,iv,imol_blend,v_blend,nb)
				flux(iv)=flux(iv)+flux0*PP%A				
			else
				call Trace2StarLines(PP,flux0,iv,imol_blend,v_blend,nb)
				flux(iv)=flux(iv)+flux0*PP%A
			endif
			! FIXME: workaround assume that the star is at the center of the grid
			PP%im_ixy(:,1)=int(npix/2)+1
			PP%im_npix=1
			if(imagecube) call AddImage(iv,0,0,flux0*PP%A,PP,1,"trace star")
		enddo
		flux2=flux2+flux0*PP%A

		flux1=0d0
		flux3=0d0

		f=sqrt((1d0+real(Bl%nvmin)*vresolution/clight)/(1d0-real(Bl%nvmin)*vresolution/clight))
		wl11=log10(lam_cont(ilam+1)/(lam*f))/log10(lam_cont(ilam+1)/lam_cont(ilam))
		wl21=1d0-wl11

		f=sqrt((1d0+real(Bl%nvmax)*vresolution/clight)/(1d0-real(Bl%nvmax)*vresolution/clight))
		wl13=log10(lam_cont(ilam+1)/(lam*f))/log10(lam_cont(ilam+1)/lam_cont(ilam))
		wl23=1d0-wl13

		flux_l1=0d0
		flux_l2=0d0
		do k=1,ngrids
			do j=1,npoints(k)
				PP => P(k,j)
				flux_l1=flux_l1+PP%flux_cont(ilam)*PP%A/real(ngrids)
				flux_l2=flux_l2+PP%flux_cont(ilam+1)*PP%A/real(ngrids)
			enddo
		enddo
		PP => path2star

		flux_l1=flux_l1+PP%flux_cont(ilam)*PP%A
		flux_l2=flux_l2+PP%flux_cont(ilam+1)*PP%A

		flux1=flux_l1**wl11*flux_l2**wl21
		flux3=flux_l1**wl13*flux_l2**wl23

		do iv=Bl%nvmin,Bl%nvmax
			fc=flux1+(flux3-flux1)*real(i-Bl%nvmin)/real(Bl%nvmax-Bl%nvmin)
			flux(iv)=flux(iv)-flux2+fc
			flux_cont(iv)=fc
		enddo

		if(Bl%n.eq.1) then
			comment = trim(Mol(Bl%L(1)%imol)%name) // "  up: " // trim(int2string(Bl%L(1)%jup,'(i5)'))
     &											  // "  low: " // trim(int2string(Bl%L(1)%jlow,'(i5)'))
		else
			count=0
			do iv=1,Bl%n
				count(Bl%L(iv)%imol)=count(Bl%L(iv)%imol)+1
			enddo
			comment = "blend of "
			do iv=1,nmol
				if(count(iv).gt.0) then
					comment=trim(comment) // trim(int2string(count(iv),'(i3)')) // " " // trim(Mol(iv)%name)
				endif
			enddo

		endif

		Bl%F=0d0
		lam_velo=lam*sqrt((1d0+vresolution/clight)/(1d0-vresolution/clight))
		dnu=dabs(clight*1d4*(1d0/lam_velo-1d0/lam))
		lam_w=0d0
		lam_w_min=1d200
		lam_w_max=0d0
		do iv=Bl%nvmin,Bl%nvmax
			lam_velo=lam*sqrt((1d0+real(iv)*vresolution/clight)/(1d0-real(iv)*vresolution/clight))
			if(lam_velo.gt.lmin_next) then

				write(20,*) lam_velo,
     &					flux(iv)*1e23/(distance*parsec)**2,
     &					real(iv)*vresolution/1d5,
     &					flux_cont(iv)*1e23/(distance*parsec)**2,
     &					trim(comment)
				Bl%F=Bl%F+dnu*(flux(iv)-flux_cont(iv))
				lam_w=lam_w+lam_velo*dnu*(flux(iv)-flux_cont(iv))
				if(lam_velo.lt.lam_w_min) lam_w_min=lam_velo
				if(lam_velo.gt.lam_w_max) lam_w_max=lam_velo
			endif
		enddo
		if(Bl%F.ne.0d0) then     ! added PW, June 13, 2022      
			lam_w=lam_w/Bl%F
		endif
		Bl%F=Bl%F*1d-3/(distance*parsec)**2
		write(21,*) lam_w,Bl%F,lam_w_min,lam_w_max,trim(comment)

		! has to be here, (i.e. before lmin_next is set, and beofre BL%next is done)
		if(imagecube) then
			! do it similar to the flux output ... find the index where we actuall start (have data)
			do iv=Bl%nvmin,Bl%nvmax
				lam_velo_imcube=lam*sqrt((1d0+real(iv)*vresolution/clight)/(1d0-real(iv)*vresolution/clight))
				if(lam_velo_imcube.gt.lmin_next) exit
			enddo
			!write(*,*) iv,Bl%nvmax,nvmin,nvmax,lmin_next,lam_velo
			imcubename="imcube" // trim(int2string(iblends,'(i0.10)')) // ".fits.gz"
			write(*,*) "Writing image cube to file ",trim(imcubename)
			!write(*,*) imcube(12,60,-1),imcube(12,42,-1),imcube(12,60,1),imcube(12,42,1)
			! currently this is per blend 
			!call writefitsfile(imcubename,imcube*1e23/(distance*parsec)**2,nvmax-nvmin+1,npix)
			!lam_velo=lam*sqrt((1d0+real(Bl%nvmin)*vresolution/clight)/(1d0-real(Bl%nvmin)*vresolution/clight))
			call writefitsfile(imcubename,imcube(:,:,iv:nvmax)*1e23/(distance*parsec)**2,nvmax-iv+1,npix,Bl%lam,lam_velo_imcube)
			
			!imcubename="imcube_hit" // trim(int2string(iblends,'(i0.10)')) // ".fits.gz"
			!call writefitsfile(imcubename,imcube_hit,1,npix)
		endif		


		lmin_next=max(lmin_next,lam_velo)

		
		endif ! if(lam.gt.lmin.and.lam.lt.lmax)

		if(iblends.lt.nblends) Bl => Bl%next

		call tellertje_time(iblends,nblends,nl,nltot,starttime)
c		call tellertje_time(iblends,nblends,iblends,nblends,starttime)

		
	enddo

	ilam=ilam+1
	do while(ilam.le.nlam)
		if(lam_cont(ilam).lt.lmax) then
			flux0=0d0
			do j=1,npoints(i)
				PP => P(i,j)
				flux0=flux0+PP%flux_cont(ilam)*PP%A
			enddo
			flux0=flux0+path2star%flux_cont(ilam)*path2star%A
			write(20,*) lam_cont(ilam),flux0*1e23/(distance*parsec)**2
		endif
		ilam=ilam+1
	enddo

	close(unit=20)
	close(unit=21)
	
	!call cpu_time(stoptime)
	call SYSTEM_CLOCK(counts, count_rate, count_max)
        stoptime = DBLE(counts)/DBLE(count_rate)

	call output("Time used for the lines:"//trim(dbl2string(stoptime-starttime,'(f8.2)'))
     &			//" s")
	call output("Time used per line:     "//trim(dbl2string((stoptime-starttime)/real(nlines),'(f8.2)'))
     &			//" s")
c	call output("Time used per line:     "//trim(dbl2string((stoptime-starttime)/real(nblends),'(f8.2)'))
c     &			//" s")


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
			S=CC%therm_l*tau_dust
c	dust scattering source function
			S=S+CC%scat_l*tau_dust

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
	


	subroutine TraceFluxLines(p0,flux,ii,vmult,imol_blend,v_blend,doit,nb,ib0)
	use GlobalSetup
	use Constants
	integer i,j,k,iv,ii,nv,ilines,imol,vmult,nb,imol_blend(nb)
	real*8 tau,exptau,flux,fact,prof,S,tau_gas,tau_dust,tau_d,tau_tot,v_blend(nb)
	type(Path) p0
	type(Cell),pointer :: CC
	type(Blend) Bl
	logical gas,doit(nb)
	real*8 rj
	
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
				if(doit(ib)) then
				imol=imol_blend(ib)
				rj=((real(vmult)*p0%v(k)+v_blend(ib))*vres_mult/vresolution-v*vres_mult)
     &					*1d5/CC%line_width(imol)
				jj=int(rj)
				if(jj.gt.-nvprofile.and.jj.lt.nvprofile) then
					if(profile_nz(jj)) then
						prof=profile(jj)
						tau_gas=prof*CC%line_abs(ib+ib0-1)
c	gas source function
						S=S+prof*CC%line_emis(ib+ib0-1)
						tau=tau+tau_gas
						gas=.true.
					endif
				endif
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
	

	subroutine DetermineBlends(nv,maxblend,nvmin0,nvmax0)
	use GlobalSetup
	use Constants
	integer ilines,ilines0,i,maxblend,nv,nvmin0,nvmax0,iv,j
	real*8 maxvshift,maxmult,v
	type(Blend),pointer :: Bl

	maxvshift=2d0*real(nv)*vresolution
	maxmult=sqrt((1d0+maxvshift/clight)/(1d0-maxvshift/clight))
	
	Bl => Blends

	nblends=0
	maxblend=0
	nvmin0=-nv
	nvmax0=nv

	do i=1,nlines
		Bl%n=0
		Bl%nvmin=-nv
		Bl%nvmax=nv
		Bl%lam=Lines(i)%lam
		do j=1,nlines
			if((j.gt.i.and.(Lines(j)%lam/Lines(i)%lam).lt.maxmult).or.(j.lt.i.and.(Lines(i)%lam/Lines(j)%lam).lt.maxmult).or.j.eq.i) then
				Bl%n=Bl%n+1
			endif
		enddo
		allocate(Bl%L(Bl%n))
		allocate(Bl%v(Bl%n))
		if(Bl%n.gt.maxblend) maxblend=Bl%n
		Bl%n=0
		do j=1,nlines
			if((j.gt.i.and.(Lines(j)%lam/Lines(i)%lam).lt.maxmult).or.(j.lt.i.and.(Lines(i)%lam/Lines(j)%lam).lt.maxmult).or.j.eq.i) then
				Bl%n=Bl%n+1
				Bl%L(Bl%n) = Lines(j)
				if(i.eq.j) then
					Bl%v(Bl%n)=0d0
				else
					f=(Lines(j)%lam/Lines(i)%lam)**2
					Bl%v(Bl%n)=clight*(f-1d0)/(f+1d0)
				endif
			endif
		enddo
		allocate(Bl%ib0(-nv:nv))
		allocate(Bl%nb0(-nv:nv))
		do iv=-nv,nv
			do j=1,Bl%n
				if(real(iv-nv)*vresolution.le.Bl%v(j)) exit
			enddo
			Bl%ib0(iv)=j
			do j=Bl%n,1,-1
				if(real(iv+nv)*vresolution.ge.Bl%v(j)) exit
			enddo
			Bl%nb0(iv)=j-Bl%ib0(iv)+1
		enddo

		v=-real(nv)*vresolution
		Bl%lmin=Bl%L(1)%lam*sqrt((1d0+v/clight)/(1d0-v/clight))
		v=real(nv)*vresolution
		Bl%lmax=Bl%L(1)%lam*sqrt((1d0+v/clight)/(1d0-v/clight))

		allocate(Bl%next)
		Bl => Bl%next
		nblends=nblends+1
	enddo
	
	print*,nlines,nblends
	
	return
	end
		


	subroutine DetermineBlendsOld(nv,maxblend,nvmin0,nvmax0)
	use GlobalSetup
	use Constants
	integer ilines,ilines0,i,maxblend,nv,nvmax,iv,nvmin0,nvmax0
	real*8 maxvshift,maxmult,v
	type(Blend),pointer :: Bl

	maxvshift=2d0*real(nv)*vresolution
	maxmult=sqrt((1d0+maxvshift/clight)/(1d0-maxvshift/clight))
	
	Bl => Blends

	ilines=2
	ilines0=1
	nblends=0
	maxblend=0
	nvmin0=-nv
	nvmax0=0
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
		if(nvmax.gt.nvmax0) nvmax0=nvmax
		Bl%nvmin=-nv
		Bl%nvmax=nvmax
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

		Bl%lam=Bl%L(1)%lam

		ilines0=ilines
		ilines=ilines+1
		allocate(Bl%next)
		Bl => Bl%next
		nblends=nblends+1
	enddo
	
	return
	end
		

	
	subroutine InterpolateLam(lam0,ilam)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 wl1,wl2,w1,w2,lam0,x1,x2
	integer ilam,i,j

	w1=(lam_cont(ilam+1)-lam0)/(lam_cont(ilam+1)-lam_cont(ilam))
	w2=1d0-w1
	wl1=log10(lam_cont(ilam+1)/lam0)/log10(lam_cont(ilam+1)/lam_cont(ilam))
	wl2=1d0-wl1
	do i=0,nR
		do j=1,nTheta
			x1=C(i,j)%kext(ilam)
			x2=C(i,j)%kext(ilam+1)
			C(i,j)%kext_l=(x1*w1)+(x2*w2)
c			x1=C(i,j)%LRF(ilam)*C(i,j)%albedo(ilam)*C(i,j)%kext(ilam)
c			x2=C(i,j)%LRF(ilam+1)*C(i,j)%albedo(ilam+1)*C(i,j)%kext(ilam+1)
c			C(i,j)%scat_l=(x1*w1+x2*w2)/C(i,j)%kext_l
c			x1=BB(ilam,C(i,j)%iT)
c			x2=BB(ilam+1,C(i,j)%iT)
c			C(i,j)%therm_l=w1*x1+w2*x2
c			x1=(1d0-C(i,j)%albedo(ilam))*C(i,j)%kext(ilam)
c			x2=(1d0-C(i,j)%albedo(ilam+1))*C(i,j)%kext(ilam+1)
c			C(i,j)%therm_l=C(i,j)%therm_l*(w1*x1+w2*x2)/C(i,j)%kext_l

c Using the source function for now.
			x1=C(i,j)%S(ilam)
			x2=C(i,j)%S(ilam+1)
			C(i,j)%therm_l=(w1*x1+w2*x2)/2d0
			C(i,j)%scat_l=(w1*x1+w2*x2)/2d0
		enddo
	enddo
	Fstar_l=(Fstar(ilam)**wl1)*(Fstar(ilam+1)**wl2)

	return
	end

	

	! ! create a regular pixel coordinate system in the image (fits) plane 
	! ! with the center pixel being (0,0) x is the horizontal axis and y the vertical one
	! ! the x,y coordinates correspond to the center of the grid
	! subroutine Imcoords(npix,Rout,im_coord)
	! use GlobalSetup
	! use Constants
	! integer, intent(in) :: npix
	! real*8, intent(in) :: Rout
	! real*8, intent(inout) :: im_coord(npix,npix) ! should be already allocated
	! real*8 :: delpix

	! delpix=2.d0*Rout*AU/real(npix)	
	! ! create the coordinates for the pixels
	! do i=1,npix ! center coordinates of pixels in cm
	! 	im_coord(:,:)=delpix/2d0+delpix*(i-1)-Rout*AU
	! enddo	

	! end subroutine im_coords


	subroutine map_pixels_to_path(igrid,npath)
	use GlobalSetup
	use Constants
	implicit none
	integer, intent(in) :: igrid,npath
	real*8 :: px(npath),py(npath),dist(npath),im_x,im_y
	real*8 :: maxx,maxy,maxr,maxrxp,maxrxm
	type(Path),pointer :: p0

	integer :: i,j,ipath,iminpath,maxnpixpath


    ! If already allocated this was done already (same grid igrid is used again), no need to map things again
	if (allocated(P(igrid,1)%im_ixy)) return


	write(*,*) "Mapping pixels to path igrid,npath",igrid,npath

    
	! the paths do not sample the whole image (just the disk)
	! In theory this could be different (e.g. if not a disk)
	!maxr=maxval(abs(P(igrid,:)%y)) ! this should always be the max r, major axis
	!write(*,*) maxr/AU
	maxr=R(nr)
	!maxrxp=maxval(P(igrid,:)%x)
	!maxrxm=minval(P(igrid,:)%x) 
	!maxr=max(maxrxp,maxrxm)

	P(igrid,:)%im_npix=0
	px(:)=P(igrid,:)%x
	py(:)=P(igrid,:)%y

	maxnpixpath=npix*int(npix/10)
	do ipath=1,npath
		! FIXME: think about how large that can be
		allocate(P(igrid,ipath)%im_ixy(2,maxnpixpath))
	enddo

	do i=1,npix
		im_x=im_coord(i)
		! only care about the positive half
		do j=int(npix/2+1),npix
			im_y=im_coord(j)
			! the Path grid does not cover a circel (due to inclination)
			! along the x is the minor axis, which has differen Rmax depending on sign
			if (sqrt((im_x)**2+(im_y)**2)>maxr) cycle
			!if ((im_x<0).and.im_x<maxrxm) cycle FIXME: doesn't really work
			!if ((im_x>0).and.im_x>maxrxp) cycle
			dist=sqrt((px-im_x)**2+(py-im_y)**2)
			! don't need squareroot
			iminpath=minloc(dist,dim=1)
!			if (i==int(npix/2)) then 
!				write(*,*) "out: ",i,j,im_x/AU,im_y/AU,dist(iminpath)/AU,iminpath,P(igrid,iminpath)%x/AU,P(igrid,iminpath)%y/AU
!			endif
       			
			P(igrid,iminpath)%im_npix=P(igrid,iminpath)%im_npix+1

			if (P(igrid,iminpath)%im_npix>maxnpixpath) then 
				write(*,*)'Problem found to many pixels for path'
				write(*,*) P(igrid,iminpath)%x/AU,P(igrid,iminpath)%y/AU,P(igrid,iminpath)%im_npix
				stop
			endif
			P(igrid,iminpath)%im_ixy(1,P(igrid,iminpath)%im_npix)=i
			P(igrid,iminpath)%im_ixy(2,P(igrid,iminpath)%im_npix)=j
		enddo
	enddo

	! now there will be paths with no pixel assigned as in case for Paths being
	! close to each other the pixel is only assigned to the closest one
	! so simply go through the paths without pixels, and assign the closest one. 

	do ipath=1,npath
		p0 => P(igrid,ipath)
		if (p0%im_npix>0) cycle
		! find the closest pixel
		p0%im_npix=1
		p0%im_ixy(1,1)=minloc(abs(im_coord-p0%x),dim=1) 
		p0%im_ixy(2,1)=minloc(abs(im_coord-p0%y),dim=1)
	enddo 

	end subroutine map_pixels_to_path
	
	subroutine AddImage(iv,i,j,flux0,p0,vmult,caller) ! FIXME: don't need i, j anymore 
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 flux0,x,y
	integer, intent(in) :: iv
	integer, intent(in) :: vmult
	type(Path), intent(in) :: p0
	character(len=*), intent(in) :: caller	! just for logging, can be removed
	integer i,j,iint,ix,iy,iym,ipix
	real*8 :: delpix,fluxperpix
 
	! simply distribute the flux for the path over all pixels equally
	fluxperpix=flux0/p0%im_npix

	do ipix=1,p0%im_npix

		ix=p0%im_ixy(1,ipix)
		iy=p0%im_ixy(2,ipix)

		! assume here x and y = 0 at the center of the image 
		if (vmult < 0) then
			iy=npix-(iy-1) ! for mirroring the whole thing
		endif

		! if (iv==1.or.iv==-1) then 	
		! 	write(*,*) "Caller: ",trim(caller)
		! 	write(*,*) "AddImage",iv,i,j,flux0,p0%x/AU,p0%y/AU,ix,iy
		! 	!write(*,*) delpix/AU,2.d0*Rout,delpix*real(npix)/AU
		! endif	
		
		imcube(iy,ix,iv)=imcube(iy,ix,iv)+fluxperpix ! P%y (Vertical) seems to be along the major axis - make it x
	enddo
	return
	end
		
	
	
	
	
	
