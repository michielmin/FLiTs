	subroutine ReadLineData()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	
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
	enddo
	close(unit=80)

	if(LTE) then
		call ComputeLTE()
	endif
	
	do i=0,nR
		do j=1,nTheta
			C(i,j)%line_width=sqrt(2d0*kb*C(i,j)%Tgas/(mp*Mol%M))
			C(i,j)%line_width=C(i,j)%line_width+0.5d0*sqrt((7.0/5.0)*kb*C(i,j)%Tgas/(mp*2.3))
			if(C(i,j)%line_width.lt.vresolution/vres_mult) C(i,j)%line_width=vresolution/vres_mult
		enddo
	enddo
	
	return
	end
	