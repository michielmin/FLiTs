	subroutine ReadLambdaFiles()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,i_low,i_up,imol
	real*8 v1,v2
	
	allocate(Mol(nmol))

	do imol=1,nmol
		open(unit=80,file=linefile(imol),RECL=6000)
		read(80,*)
		read(80,*) Mol(imol)%name
		read(80,*)
		read(80,*) Mol(imol)%M
		read(80,*)
		read(80,*) Mol(imol)%nlevels
		read(80,*)
		allocate(Mol(imol)%E(Mol(imol)%nlevels))
		allocate(Mol(imol)%g(Mol(imol)%nlevels))
		do i=1,Mol(imol)%nlevels
			read(80,*) j,Mol(imol)%E(i),Mol(imol)%g(i)
		enddo
		read(80,*)
		read(80,*) Mol(imol)%nlines
		read(80,*)
		allocate(Mol(imol)%L(Mol(imol)%nlines))
		do i=1,Mol(imol)%nlines
			read(80,*) j,Mol(imol)%L(i)%jup,Mol(imol)%L(i)%jlow,Mol(imol)%L(i)%Aul,
     &				Mol(imol)%L(i)%freq,Mol(imol)%E(Mol(imol)%L(i)%jup)	!Mol%L(i)%Eup
			Mol(imol)%L(i)%freq=Mol(imol)%L(i)%freq*1d9
			Mol(imol)%L(i)%lam=clight*1d4/(Mol(imol)%L(i)%freq)
			Mol(imol)%L(i)%imol=imol

			i_low=Mol(imol)%L(i)%jlow
			i_up=Mol(imol)%L(i)%jup
			Mol(imol)%L(i)%Bul=Mol(imol)%L(i)%Aul/(2d0*hplanck*Mol(imol)%L(i)%freq**3/clight**2)
			Mol(imol)%L(i)%Blu=Mol(imol)%L(i)%Bul*Mol(imol)%g(i_up)/Mol(imol)%g(i_low)
		enddo
		close(unit=80)
	enddo

	return
	end
	
