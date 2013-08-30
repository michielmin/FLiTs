	subroutine ReadLineData()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,i_low,i_up
	
	open(unit=80,file=linefile,RECL=6000)
	read(80,*)
	read(80,*) Mol%name
	read(80,*)
	read(80,*) Mol%M
	read(80,*)
	read(80,*) Mol%nlevels
	read(80,*)
	allocate(Mol%E(Mol%nlevels))
	allocate(Mol%g(Mol%nlevels))
	do i=1,Mol%nlevels
		read(80,*) j,Mol%E(i),Mol%g(i)
	enddo
	read(80,*)
	read(80,*) Mol%nlines
	read(80,*)
	allocate(Mol%L(Mol%nlines))
	do i=1,Mol%nlines
		read(80,*) j,Mol%L(i)%jup,Mol%L(i)%jlow,Mol%L(i)%Aul,Mol%L(i)%freq,Mol%E(Mol%L(i)%jup)	!Mol%L(i)%Eup
		Mol%L(i)%freq=Mol%L(i)%freq*1d9
		Mol%L(i)%lam=clight*1d4/(Mol%L(i)%freq)

		i_low=Mol%L(i)%jlow
		i_up=Mol%L(i)%jup
		Mol%L(i)%Bul=Mol%L(i)%Aul/(2d0*hplanck*Mol%L(i)%freq**3/clight**2)
		Mol%L(i)%Blu=Mol%L(i)%Bul*Mol%g(i_up)/Mol%g(i_low)
	enddo
	close(unit=80)


	if(popfile.ne.' ') 	call ReadPopData()

	if(LTE.or.popfile.eq.' ') call ComputeLTE()
	
	do i=0,nR
		do j=1,nTheta
			if(popfile.eq.' ') then
				C(i,j)%line_width=sqrt(2d0*kb*C(i,j)%Tgas/(mp*Mol%M))
				C(i,j)%line_width=C(i,j)%line_width+0.5d0*sqrt((7.0/5.0)*kb*C(i,j)%Tgas/(mp*2.3))
			endif
			if(C(i,j)%line_width.lt.vresolution/vres_mult) C(i,j)%line_width=vresolution/vres_mult
		enddo
	enddo
	
	return
	end
	
	
	
c=========================================================================================
c This subroutine reads in a forMCFOST.fits type of file generated by ProDiMo
c It fills the population levels
c=========================================================================================
	subroutine ReadPopData()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,ivars,i,j,l,k,naxis,npopname,ipop
	character*7 vars(10),hdu
	real,allocatable :: array(:,:,:,:)
	real*8,allocatable :: array_d(:,:,:,:)
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real  :: nullval
	real*8  :: nullval_d,tot
	logical*4 :: anynull
	integer*4, dimension(4) :: naxes
	character*80 comment,errmessage
	character*30 errtext,popname(20)

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	call ftopen(unit,popfile,readwrite,blocksize,status)
	if (status /= 0) then
		call output("Population file not found "//trim(popfile))
		call output("==================================================================")
		stop
	endif
	group=1
	firstpix=1
	nullval=-999
	nullval_d=-999

	call output("Reading level populations from: "//trim(popfile))

c set default names of the species
	popname(1) = "C+"
	popname(2) = "O"
	popname(3) = "CO"
	popname(4) = "o-H2O"
	popname(5) = "p-H2O"  
	npopname=5

	do i=1,npopname
		if(trim(Mol%name).eq.trim(popname(i))) ipop=i
	enddo

	!------------------------------------------------------------------------
	! HDU0 : grid
	!------------------------------------------------------------------------
	! Skip this, it is already done

	!------------------------------------------------------------------------------
	! HDU 2: Gas Temperature 
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
	allocate(array_d(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpvd(unit,group,firstpix,npixels,nullval_d,array_d,anynull,status)

	do i=1,nR
		do j=1,nTheta
			C(i,j)%Tgas=array_d(i,nTheta+1-j,1,1)
		enddo
	enddo

	deallocate(array_d)
	
	!------------------------------------------------------------------------------
	! HDU 3 : Molecular particle densities [1/cm^3]
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
	allocate(array_d(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpvd(unit,group,firstpix,npixels,nullval_d,array_d,anynull,status)

	do i=1,nR
		do j=1,nTheta
			if(C(i,j)%dens.gt.1d-50) then
				C(i,j)%abun=array_d(ipop,i,nTheta+1-j,1)*Mol%M*mp/C(i,j)%dens
			else
				C(i,j)%abun=1d-4
			endif
		enddo
	enddo

	deallocate(array_d)

	!------------------------------------------------------------------------------
	! HDU 4 : Line broadening parameter (should I read this in?)
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
	allocate(array_d(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpvd(unit,group,firstpix,npixels,nullval_d,array_d,anynull,status)

	do i=1,nR
		do j=1,nTheta
			C(i,j)%line_width=array_d(ipop,i,nTheta+1-j,1)*1d5
		enddo
	enddo

	deallocate(array_d)
	
	!------------------------------------------------------------------------------
	! HDU 5... : level populations
	!------------------------------------------------------------------------------

	do i=1,ipop
	!  move to next hdu
		call ftmrhd(unit,1,hdutype,status)
		if(status.ne.0) then
			status=0
			goto 1
		endif
	enddo

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

	if(naxes(1).lt.Mol%nlevels.and..not.LTE) then
		call output("Assuming levels above " // int2string(naxes(1),'(i4)') // "unpopulated")
	endif

	if(naxes(1).gt.Mol%nlevels) naxes(1)=Mol%nlevels
	do i=1,nR
		do j=1,nTheta
			allocate(C(i,j)%npop(Mol%nlevels))
			C(i,j)%npop(1:Mol%nlevels)=0d0
			tot=0d0
			do k=1,naxes(1)
				C(i,j)%npop(k)=array(k,i,nTheta+1-j,1)
			enddo
c			do k=2,naxes(1)
c				C(i,j)%npop(k)=C(i,j)%npop(k-1)*C(i,j)%npop(k)
c			enddo
c			C(i,j)%npop(1)=1d0-sum(C(i,j)%npop(2:naxes(1)))
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

	do j=1,nTheta
		allocate(C(0,j)%npop(Mol%nlevels))
		C(0,j)%npop(1:Mol%nlevels)=C(1,j)%npop(1:Mol%nlevels)
		C(0,j)%abun=C(1,j)%abun
	enddo
	
	return
	end
	


