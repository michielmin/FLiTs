	subroutine ReadStructure()
	use GlobalSetup
	IMPLICIT NONE
	integer i
	
	if(structtype.eq.1) then
		call ReadMCMaxStructure()
	else if(structtype.eq.2) then
c		call ReadProDiMoStructure()
	endif
	
	
	open(unit=20,file='output.dat',RECL=6000)
	do i=1,nR
		write(20,*) R(i),C(i,nTheta)%dens
	enddo
	close(unit=20)
	
	return
	end
	
	
	subroutine ReadMCMaxStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	
	call readfitsMCMax(structfile)
	
	
	return
	end
	
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine readfitsMCMax(filename)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,l,iopac,naxis,nhdu
	character*7 vars(10),hdu
	character*500 filename
	logical doalloc,truefalse
	real*8,allocatable :: array(:,:,:,:)
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real*8  :: nullval
	logical*4 :: anynull
	integer*4, dimension(4) :: naxes
	character*80 comment,errmessage
	character*30 errtext
	real*8,allocatable :: R_av(:),theta_av(:)

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
	
	allocate(C(nR,nTheta))
	allocate(R(nR+1))
	allocate(Theta(nTheta+1))
 
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
	deallocate(R_av)
	deallocate(theta_av)

	do ivars=1,nhdu
		!  move to next hdu
		call ftmrhd(unit,1,hdutype,status)
		if(status.ne.0) then
			status=0
			goto 1
		endif

		select case (vars(ivars))
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
						C(i,j)%T=array(i,j,1,1)
					enddo
				enddo
			case ('GASDENS')
				do i=1,nR
					do j=1,nTheta
						C(i,j)%dens=array(i,j,1,1)
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


