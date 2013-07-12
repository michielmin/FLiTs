	subroutine SetupPaths()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ip,jp,i,j,ir
	real*8 inc_min
	real*8,allocatable :: imR(:),imPhi(:)
	
	inc_min=5d0
	
	call output("==================================================================")
	call output("Setup up the paths for raytracing")
	
	nImPhi=abs(sin(inc*pi/180d0))*(C(1,nTheta)%v/vresolution)*3d0
	if(nImPhi.lt.10) nImPhi=10
	if(nImPhi.gt.180) nImPhi=180
	
	if(inc.gt.inc_min) then
		nImR=nR*2+nTheta*2
	else
		nImR=nR
	endif

	nImR=nImR+abs(sin(inc*pi/180d0))*(C(1,nTheta)%v/vresolution)*3d0

	call output("Number of radial image points: "//trim(int2string(nImR,'(i5)')))
	call output("Number of phi image points:    "//trim(int2string(nImPhi,'(i5)')))
	
	allocate(imR(nImR))
	allocate(imPhi(nImPhi))
	
	ir=0
	do i=1,nR
		ir=ir+1
		imR(ir)=R_av(i)
	enddo
	if(inc.gt.inc_min) then
		do i=1,nR
			ir=ir+1
			imR(ir)=abs(R_av(i)*sin(inc*pi/180d0))
		enddo
		do i=1,nTheta
			ir=ir+1
			imR(ir)=abs(R(1)*sin(theta_av(i)-inc))
			ir=ir+1
			imR(ir)=abs(R(1)*sin(pi-theta_av(i)+inc))
		enddo
	endif

	j=nImR-ir
	do i=1,j
		ir=ir+1
		imR(ir)=10d0**(log10(Rstar*Rsun)+log10(R(nR+1)/(Rstar*Rsun))*real(i-1)/real(j-1))
	enddo
	
	call sort(imR,nImR)
	
	do i=1,nImPhi
		ImPhi(i)=pi*(real(i)-0.5)/real(nImPhi)
	enddo

	allocate(P(nImR,nImPhi))
	do i=1,nImR
		if(i.ne.1) P(i,1)%R1=sqrt(ImR(i-1)*ImR(i))
		if(i.ne.nImR) P(i,1)%R2=sqrt(ImR(i)*ImR(i+1))
	enddo
	P(1,1)%R1=ImR(1)**2/P(1,1)%R2
	P(nImR,1)%R2=ImR(nImR)**2/P(nImR,1)%R1

	do i=1,nImR
		do j=1,nImPhi
			P(i,j)%R=ImR(i)
			P(i,j)%Phi=ImPhi(j)
			P(i,j)%phi1=pi*real(j-1)/real(nImPhi)
			P(i,j)%phi2=pi*real(j)/real(nImPhi)
			P(i,j)%R1=P(i,1)%R1
			P(i,j)%R2=P(i,1)%R2
			P(i,j)%A=pi*(P(i,j)%R2**2-P(i,j)%R1**2)/real(2*nImPhi)
		enddo
	enddo

	return
	end
	