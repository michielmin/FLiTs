	subroutine SortLines()
	use GlobalSetup
	IMPLICIT NONE
	type(Line),allocatable :: Line_temp(:)
	type(Line) Line_min
	real*8 lam_min
	integer i,j,imin,nlines_temp
	
	allocate(Line_temp(Mol%nlines))

	nlines_temp=0
	do j=1,Mol%nlines-1
		lam_min=Mol%L(j)%lam
		imin=j
		do i=j,Mol%nlines
			if(Mol%L(i)%lam.lt.lam_min) then
				lam_min=Mol%L(i)%lam
				imin=i
			endif
		enddo
		Line_min=Mol%L(j)
		Mol%L(j)=Mol%L(imin)
		Mol%L(imin)=Line_min
		if(Mol%L(j)%lam.gt.lmin.and.Mol%L(j)%lam.lt.lmax) then
			nlines_temp=nlines_temp+1
			Line_temp(nlines_temp)=Mol%L(j)
		endif
	enddo

	deallocate(Mol%L)

	Mol%nlines=nlines_temp
	allocate(Mol%L(Mol%nlines))
	Mol%L(1:Mol%nlines)=Line_temp(1:Mol%nlines)

	deallocate(Line_temp)

	return
	end
