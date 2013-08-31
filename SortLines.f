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
				if(Mol(imol)%L(i)%lam.lt.lam_min.and.
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
