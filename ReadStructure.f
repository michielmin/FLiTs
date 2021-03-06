	subroutine ReadStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	
	call ReadForProDiMo()
	
	do i=0,nR
		do j=0,nTheta
			C(i,j)%iT=C(i,j)%Tdust+0.5d0
			if(C(i,j)%iT.lt.1) C(i,j)%iT=1
			if(C(i,j)%iT.gt.MAXT) C(i,j)%iT=MAXT
		enddo
	enddo
	
	
c everything is read in now
c output the setup to the screen and the log file
	call output("==================================================================")

	rlines=1d0/(sqrt((1d0+vresolution/clight)/(1d0-vresolution/clight))-1d0)

	call output("Mass of the star:   "//trim(dbl2string(Mstar,'(f13.4)'))//" Msun")
	call output("Minimum wavelength: "//trim(dbl2string(lmin,'(f13.4)'))//" micron")
	call output("Maximum wavelength: "//trim(dbl2string(lmax,'(f13.4)'))//" micron")
	call output("Resolution lines:   "//trim(dbl2string(rlines,'(f13.4)'))//" (dlam/lam)")
	call output("Velocity resolution:"//trim(dbl2string(vresolution/1d5,'(f13.4)'))//" km/s")
	call output("Inclination angle:  "//trim(dbl2string(inc,'(f13.4)'))//" degrees")

	call output("==================================================================")
	do i=1,nmol
		call output("Line file: "//trim(linefile(i)))
	enddo

	return
	end


c=========================================================================================
c This subroutine reads in a forProDiMo.fits file generated by MCMax
c It fills the arrays:
c dens,T,kabs,kext,ksca,v,lam_cont,R,Theta
c=========================================================================================
	subroutine ReadForProDiMo()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,ivars,i,j,l,naxis,nhdu
	character*7 vars(10),hdu
	real,allocatable :: array(:,:,:,:)
	real*8,allocatable :: array_d(:,:,:,:)
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real  :: nullval
	real*8  :: nullval_d,xx,zz,rr
	logical*4 :: anynull
	integer*4, dimension(4) :: naxes
	character*80 comment,errmessage
	character*30 errtext

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	call ftopen(unit,structfile,readwrite,blocksize,status)
	if (status /= 0) then
		write(*,'("Density file not found")')
		write(9,'("Density file not found")')
		print*,trim(structfile)
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif
	group=1
	firstpix=1
	nullval=-999
	nullval_d=-999


	!------------------------------------------------------------------------
	! HDU0 : grid
	!------------------------------------------------------------------------
	! Check dimensions
	call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

	npixels=naxes(1)*naxes(2)*naxes(3)

	! Read model info

	call ftgkyd(unit,'Rin',Rin,comment,status)
	call ftgkyd(unit,'Rout',Rout,comment,status)

	call ftgkyd(unit,'Mstar',Mstar,comment,status)
	call ftgkyd(unit,'Rstar',Rstar,comment,status)

	call ftgkyd(unit,'distance',distance,comment,status)

	call ftgkyj(unit,'n_rad',nR,comment,status)
	nR=nR+1
	call ftgkyj(unit,'nz',nTheta,comment,status)
	
	allocate(C(0:nR,0:nTheta))
	allocate(R(0:nR+1))
	allocate(Theta(0:nTheta+1))

	! read_image
	allocate(array_d(nR-1,nTheta,2,1))
	allocate(R_av(0:nR))
	allocate(theta_av(0:nTheta))

	call ftgpvd(unit,group,firstpix,npixels,nullval_d,array_d,anynull,status)

	do j=1,nTheta
		theta_av(j)=acos(array_d(1,nTheta+1-j,2,1)/array_d(1,1,1,1))
	enddo

c in the theta grid we actually store cos(theta) for convenience
	Theta(0)=1d0
	do i=2,nTheta
		Theta(i)=cos((theta_av(i-1)+theta_av(i))/2d0)
	enddo
	Theta(nTheta+1)=0d0
	if(cylindrical) then
		Theta(1)=2d0*cos(theta_av(1))-Theta(2)
		if(Theta(1).gt.1d0) Theta(1)=1d0
		theta_av(0)=acos(Theta(1))/2d0
	else
		Theta(1)=1d0
		theta_av(0)=0d0
	endif

	do i=2,nR-2
		R_av(i)=array_d(i,1,1,1)*AU
	enddo

	Rin=array_d(1,1,1,1)
	call output("Adjusting Rin to:  "//trim(dbl2string(Rin,'(f8.3)')) //" AU")
	Rout=array_d(nR-1,1,1,1)
	call output("Adjusting Rout to: "//trim(dbl2string(Rout,'(f8.3)')) //" AU")

	R(0)=Rstar*Rsun
	R(1)=Rin*AU
	R_av(0)=(R(1)+R(0))/2d0
	R_av(1)=sqrt(R_av(2)*R(1))
	R(2)=sqrt(R_av(1)*R_av(2))
	do i=3,nR-2
		R(i)=10d0**((2d0*log10(R_av(i-1))-log10(R(i-1))))
		if(R_av(i).lt.R(i).or.R_av(i-1).gt.R(i)) then
			R(i)=(R_av(i-1)+R_av(i))/2d0
		endif
	enddo
	R(nR)=Rout*AU
	R(nR-1)=sqrt(R_av(nR-2)*R(nR))

	if(cylindrical) then
		if(Theta(1).lt.1d0) then
			Rout=Rout*1.0001/sin(acos(Theta(1)))
		else
			call output("Grid seems to be spherical")
			call output("SWITCHING TO SPHERICAL GRID")
			cylindrical=.false.
			Rout=Rout*1.0001
		endif
	else
		Rout=Rout*1.0001
	endif

	R(nR+1)=Rout*AU

	R_av(nR-1)=sqrt(R(nR-1)*R(nR))
	R_av(nR)=sqrt(R(nR)*R(nR+1))

	call sort(R(1:nR+1),nR+1)

	allocate(R_sphere(0:nR+1))
	allocate(R_av_sphere(0:nR+1))
	if(cylindrical) then
		do i=0,nR
			R_sphere(i)=R(i)/sin(acos(Theta(1)))
			R_av_sphere(i)=R_av(i)/sin(acos(Theta(1)))
		enddo
		R_sphere(nR+1)=R(nR+1)
		R_av_sphere(nR)=sqrt(R_sphere(nR)*R_sphere(nR+1))
	else
		do i=0,nR
			R_sphere(i)=R(i)
			R_av_sphere(i)=R_av(i)
		enddo
		R_sphere(nR+1)=R(nR+1)
	endif

	deallocate(array_d)

	!------------------------------------------------------------------------------
	! HDU 2: Temperature 
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	naxis=2

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

	do i=naxis+1,4
		naxes(i)=1
	enddo
	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! read_image
	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpve(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,nR-1
		do j=1,nTheta
			C(i,j)%Tdust=array(i,nTheta+1-j,1,1)
			C(i,j)%Tgas=C(i,j)%Tdust
		enddo
	enddo

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 3 : Longueurs d'onde
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	naxis=1

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

	nlam=naxes(1)
	do i=0,nR
		do j=0,nTheta
			allocate(C(i,j)%kabs(nlam))
			allocate(C(i,j)%albedo(nlam))
			allocate(C(i,j)%kext(nlam))
			allocate(C(i,j)%LRF(nlam))
		enddo
	enddo
 	allocate(lam_cont(nlam))
 	allocate(Fstar(nlam))
 
	do i=naxis+1,4
		naxes(i)=1
	enddo
	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! read_image
	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpve(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,nlam
		lam_cont(i)=array(i,1,1,1)
	enddo

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 4 : Spectre stellaire
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	naxis=1

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

	do i=naxis+1,4
		naxes(i)=1
	enddo
	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! read_image
	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpve(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,nlam
		Fstar(i)=array(i,1,1,1)
		Fstar(i)=Fstar(i)*lam_cont(i)*1d3*1d-4/clight
	enddo

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 5 : Spectre ISM (input)
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	!------------------------------------------------------------------------------
	! HDU 6 : Champ de radiation en W.m-2 (lambda.F_lambda)  
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	naxis=3

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

	do i=naxis+1,4
		naxes(i)=1
	enddo
	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! read_image
	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpve(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,nR-1
		do j=1,nTheta
			do l=1,nlam
				C(i,j)%LRF(l)=array(i,nTheta+1-j,l,1)
				C(i,j)%LRF(l)=C(i,j)%LRF(l)*lam_cont(l)*1d3*1d-4/clight
			enddo
		enddo
	enddo

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 7 : Statistique du champ de radiation (nombre de paquet)
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	!------------------------------------------------------------------------------
	! HDU 8 : Champ de radiation ISM en W.m-2 (lambda.F_lambda)  
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	!------------------------------------------------------------------------------
	! HDU 9 : Statistique du champ de radiation ISM (nombre de paquet)
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	!------------------------------------------------------------------------------
	! HDU 10 : Densite de gaz pour un rapport de masse de 100 / poussiere
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	naxis=2

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

	do i=naxis+1,4
		naxes(i)=1
	enddo
	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! read_image
	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpve(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,nR-1
		do j=1,nTheta
			C(i,j)%dens=array(i,nTheta+1-j,1,1)
			if(C(i,j)%dens.lt.1d-50) C(i,j)%dens=1d-60
		enddo
	enddo

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 11 : Opacites
	!------------------------------------------------------------------------------

	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)
	if(status.ne.0) then
		status=0
		goto 1
	endif

	naxis=4

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

	do i=naxis+1,4
		naxes(i)=1
	enddo
	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! read_image
	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpve(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,nR-1
		do j=1,nTheta
			do l=1,nlam
				C(i,j)%kext(l)=array(i,nTheta+1-j,1,l)/AU
				C(i,j)%kabs(l)=array(i,nTheta+1-j,2,l)/AU
				if(C(i,j)%kext(l).gt.1d-150) then
					C(i,j)%albedo(l)=(C(i,j)%kext(l)-C(i,j)%kabs(l))/C(i,j)%kext(l)
				else
					C(i,j)%albedo=0.5d0
				endif
			enddo
		enddo
	enddo

	deallocate(array)


1	continue

	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	!  Get the text string which describes the error
	if (status > 0) then
	   call ftgerr(status,errtext)
	   print *,'FITSIO Error Status =',status,': ',errtext

	   !  Read and print out all the error messages on the FITSIO stack
	   call ftgmsg(errmessage)
	   do while (errmessage .ne. ' ')
		  print *,errmessage
		  call ftgmsg(errmessage)
	   end do
	endif

	
	do i=1,nR-1
		do j=1,nTheta
c			C(i,j)%v=sqrt(G*Mstar*Msun*sin(theta_av(j))/R_av(i))
c			C(i,j)%v=sqrt(G*Mstar*Msun*sin(theta_av(j))**2/R_av(i))
			xx=R_av(i)
			zz=xx/tan(theta_av(j))
			rr=sqrt(xx*xx+zz*zz)
			C(i,j)%v=sqrt(G*Mstar*Msun*xx*xx/(rr*rr*rr))
		enddo
	enddo
	i=0
	do j=0,nTheta
		C(i,j)%kext(1:nlam)=1d-70
		C(i,j)%kabs(1:nlam)=1d-70
c		C(i,j)%v=sqrt(G*(Mstar*Msun)*sin(theta_av(j))/(sqrt(R_av(1)*Rstar*Rsun)))
c		C(i,j)%v=sqrt(G*Mstar*Msun*sin(theta_av(j))**2/(sqrt(R_av(1)*Rstar*Rsun)))
		xx=R_av(1)
		zz=xx/tan(theta_av(j))
		rr=sqrt(xx*xx+zz*zz)
		C(i,j)%v=sqrt(G*Mstar*Msun*xx*xx/(rr*rr*rr))
		C(i,j)%dens=1d-60
	enddo
	i=nR
	do j=0,nTheta
		C(i,j)%kext(1:nlam)=1d-70
		C(i,j)%kabs(1:nlam)=1d-70
c		C(i,j)%v=sqrt(G*(Mstar*Msun)*sin(theta_av(j))/(sqrt(R_av(1)*Rstar*Rsun)))
c		C(i,j)%v=sqrt(G*Mstar*Msun*sin(theta_av(j))**2/(sqrt(R_av(1)*Rstar*Rsun)))
		xx=R_av(i)
		zz=xx/tan(theta_av(j))
		rr=sqrt(xx*xx+zz*zz)
		C(i,j)%v=sqrt(G*Mstar*Msun*xx*xx/(rr*rr*rr))
		C(i,j)%dens=1d-60
	enddo

	j=0
	do i=1,nR
		C(i,j)%kext(1:nlam)=1d-70
		C(i,j)%kabs(1:nlam)=1d-70
c		C(i,j)%v=sqrt(G*(Mstar*Msun)*sin(theta_av(j))/(sqrt(R_av(1)*Rstar*Rsun)))
c		C(i,j)%v=sqrt(G*Mstar*Msun*sin(theta_av(j))**2/(sqrt(R_av(1)*Rstar*Rsun)))
		xx=R_av(i)
		zz=xx/tan(theta_av(j))
		rr=sqrt(xx*xx+zz*zz)
		C(i,j)%v=sqrt(G*Mstar*Msun*xx*xx/(rr*rr*rr))
		C(i,j)%dens=1d-60
	enddo
	
	return
	end
	





c=========================================================================================
c This subroutine reads in a forFLiTs.fits file generated by MCMax
c It fills the arrays:
c dens,T,kabs,kext,ksca,v,lam_cont,R,Theta
c=========================================================================================
	subroutine ReadMCMaxStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	
	call readfitsMCMax(structfile)
	
	do i=1,nR
		do j=1,nTheta
			C(i,j)%kext(1:nlam)=C(i,j)%kext(1:nlam)*C(i,j)%dens/100d0
			C(i,j)%kabs(1:nlam)=C(i,j)%kabs(1:nlam)*C(i,j)%dens/100d0
			C(i,j)%v=sqrt(G*(Mstar*Msun)*sin(theta_av(j))/R_av(i))
		enddo
	enddo
	i=0
	do j=1,nTheta
		C(i,j)%kext(1:nlam)=1d-70
		C(i,j)%kabs(1:nlam)=1d-70
		C(i,j)%v=sqrt(G*(Mstar*Msun)*sin(theta_av(j))/(sqrt(R_av(1)*Rstar*Rsun)))
	enddo
	
	return
	end
	
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine readfitsMCMax(filename)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,ivars,i,j,l,naxis,nhdu
	character*7 vars(10),hdu
	character*500 filename
	real*8,allocatable :: array(:,:,:,:)
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real*8  :: nullval
	logical*4 :: anynull
	integer*4, dimension(4) :: naxes
	character*80 comment,errmessage
	character*30 errtext

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	call ftopen(unit,filename,readwrite,blocksize,status)
	if (status /= 0) then
		write(*,'("Density file not found")')
		write(9,'("Density file not found")')
		print*,trim(filename)
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif
	group=1
	firstpix=1
	nullval=-999


	!------------------------------------------------------------------------
	! HDU0 : grid
	!------------------------------------------------------------------------
	! Check dimensions
	call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

	npixels=naxes(1)*naxes(2)*naxes(3)

	! Read model info

	call ftgkyd(unit,'Rin',Rin,comment,status)
	call ftgkyd(unit,'Rout',Rout,comment,status)

	call ftgkyj(unit,'nR',nR,comment,status)
	call ftgkyj(unit,'nlam',nlam,comment,status)
	call ftgkyj(unit,'nTheta',nTheta,comment,status)
	call ftgkyj(unit,'nHDU',nhdu,comment,status)
	if(status.ne.0) then
		nhdu=nvars
		status=0
	endif
	do i=1,nhdu
		write(hdu,'("HDU",i2)') i
		call ftgkys(unit,hdu,vars(i),comment,status)
	enddo
	
	allocate(C(0:nR,nTheta))
	allocate(R(nR+1))
	allocate(Theta(nTheta+1))
	do i=0,nR
		do j=1,nTheta
			allocate(C(i,j)%kabs(nlam))
			allocate(C(i,j)%albedo(nlam))
			allocate(C(i,j)%kext(nlam))
		enddo
	enddo
 	allocate(lam_cont(nlam))
 
	! read_image
	allocate(array(nR,nTheta,2,1))
	allocate(R_av(nR))
	allocate(theta_av(nTheta))

	call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,nR
		R_av(i)=array(i,1,1,1)*AU
	enddo

	R(1)=Rin*AU
	do i=2,nR
		R(i)=10d0**((2d0*log10(R_av(i-1))-log10(R(i-1))))
		if(R_av(i).lt.R(i).or.R_av(i-1).gt.R(i)) then
			R(i)=(R_av(i-1)+R_av(i))/2d0
		endif
	enddo
	R(nR+1)=Rout*AU
	call sort(R(1:nR+1),nR+1)

	do j=1,nTheta
		theta_av(j)=array(1,j,2,1)
	enddo


c in the theta grid we actually store cos(theta) for convenience
	Theta(1)=1d0
	do i=2,nTheta
		Theta(i)=cos((theta_av(i-1)+theta_av(i))/2d0)
	enddo
	Theta(nTheta+1)=0d0

	deallocate(array)

	do ivars=1,nhdu
		!  move to next hdu
		call ftmrhd(unit,1,hdutype,status)
		if(status.ne.0) then
			status=0
			goto 1
		endif

		select case (vars(ivars))
			case('LAM')
				naxis=1
			case('OPACITY')
				naxis=4
			case default
				naxis=2
		end select

		! Check dimensions
		call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

		do i=naxis+1,4
			naxes(i)=1
		enddo
		npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

		! read_image
		allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

		call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)

		select case (vars(ivars))
			case ('TEMP')
				do i=1,nR
					do j=1,nTheta
						C(i,j)%Tdust=array(i,j,1,1)
						C(i,j)%Tgas=C(i,j)%Tdust
					enddo
				enddo
			case ('GASDENS')
				do i=1,nR
					do j=1,nTheta
						C(i,j)%dens=array(i,j,1,1)
					enddo
				enddo
			case ('LAM')
				do i=1,nlam
					lam_cont(i)=array(i,1,1,1)
				enddo
			case ('OPACITY')
				do i=1,nR
					do j=1,nTheta
						do l=1,nlam
							C(i,j)%kext(l)=array(i,j,l,1)
							C(i,j)%kabs(l)=array(i,j,l,2)
							C(i,j)%albedo(l)=array(i,j,l,3)/C(i,j)%kext(l)
						enddo
					enddo
				enddo
		end select

		deallocate(array)
	enddo
	
1	continue

	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	!  Get the text string which describes the error
	if (status > 0) then
	   call ftgerr(status,errtext)
	   print *,'FITSIO Error Status =',status,': ',errtext

	   !  Read and print out all the error messages on the FITSIO stack
	   call ftgmsg(errmessage)
	   do while (errmessage .ne. ' ')
		  print *,errmessage
		  call ftgmsg(errmessage)
	   end do
	endif

	return
	end





