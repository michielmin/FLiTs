	subroutine SortLines()
	use GlobalSetup
	IMPLICIT NONE
	real*8 lam_min
	integer i,j,imin,imolmin,maxlines,imol
	logical found
	logical,allocatable :: used(:,:)
	
	call output("Sorting the lines")
	
	nlines=0
	maxlines=0
	do imol=1,nmol
		nlines=nlines+Mol(imol)%nlines
		if(Mol(imol)%nlines.gt.maxlines) maxlines=Mol(imol)%nlines
	enddo
	allocate(used(nmol,maxlines))
	used=.false.

	allocate(Lines(nlines))

	nlines=0

	found=.true.
	do while(found)
		found=.false.
		lam_min=lmax
		do imol=1,nmol
			do i=1,Mol(imol)%nlines
				if(Mol(imol)%L(i)%lam.le.lam_min.and.
     &				Mol(imol)%L(i)%lam.gt.lmin.and..not.used(imol,i)) then
					lam_min=Mol(imol)%L(i)%lam
					imin=i
					imolmin=imol
					found=.true.
				endif
			enddo
		enddo
		if(found) then
			nlines=nlines+1
			Lines(nlines)=Mol(imolmin)%L(imin)
			used(imolmin,imin)=.true.
		endif
	enddo

	return
	end

	subroutine LineList()
	use GlobalSetup
	IMPLICIT NONE
	integer i
	
	call output("Writing line list")
	
	open(unit=20,file='LineList.txt',RECL=1000)
	do i=1,nlines
		write(20,'(f16.8,"  ",a8,"  ",i8,i8)') Lines(i)%lam,trim(Mol(Lines(i)%imol)%name),Lines(i)%jup,Lines(i)%jlow
	enddo
	close(unit=20)
	
	return
	end
	



