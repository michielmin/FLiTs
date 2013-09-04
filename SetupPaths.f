	subroutine SetupPaths()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ip,jp,i,j,k,ir,nRreduce,ilam,imol,ntheta_reduce,ninc
	real*8 inc_min,ct,res_inc,x,y
	real*8,allocatable :: imR(:),imPhi(:)
	type(Tracer) trac
	type(Path),pointer :: PP

	ilam1=1
	ilam2=nlam
	do ilam=1,nlam
		if(lam_cont(ilam).lt.lmin) ilam1=ilam
		if(lam_cont(nlam+1-ilam).gt.lmax) ilam2=nlam+1-ilam
	enddo
	
c for inclinations smaller don't use the additional radial points
	inc_min=1d0
c increase the resolution in velocity by this factor
	res_inc=1d0

	nrReduce=2
	ntheta_reduce=6
	ninc=2
	
	call output("==================================================================")
	call output("Setup up the paths for raytracing")
	
	nImPhi=abs(sin(inc*pi/180d0))*(C(1,nTheta)%v/vresolution)*res_inc
	if(nImPhi.lt.30) nImPhi=30
	if(nImPhi.gt.90) nImPhi=90
	
	if(inc.gt.inc_min) then
		nImR=nR*ninc/nRreduce+nTheta*2/ntheta_reduce
	else
		nImR=nR/nRreduce
	endif

	nImR=nImR+abs(sin(inc*pi/180d0))*(C(1,nTheta)%v/vresolution)*res_inc/(nImPhi/3)

	call output("Number of radial image points: "//trim(int2string(nImR,'(i5)')))
	call output("Number of phi image points:    "//trim(int2string(nImPhi,'(i5)')))
	
	allocate(imR(nImR))
	allocate(imPhi(nImPhi))
	
	ir=0
	if(inc.gt.inc_min) then
		do i=1,nR,nRreduce
			do j=1,ninc
				x=R_av(i)*cos(pi*real(j-1)/real(ninc-1)/2d0)
				y=R_av(i)*sin(pi*real(j-1)/real(ninc-1)/2d0)
				ir=ir+1
				imR(ir)=sqrt(abs(x**2+(y*sin(inc*pi/180d0))**2))
			enddo
		enddo
		do i=1,nTheta,ntheta_reduce
			ir=ir+1
			imR(ir)=abs(R(1)*sin(theta_av(i)+inc*pi/180d0))
			ir=ir+1
			imR(ir)=abs(R(1)*sin(theta_av(i)-inc*pi/180d0))
		enddo
	else
		do i=1,nR,nRreduce
			do j=1,ninc
				ir=ir+1
				imR(ir)=R_av(i)
			enddo
		enddo	
	endif

	j=nImR-ir
	do i=1,j
		ir=ir+1
		imR(ir)=10d0**(log10(Rstar*Rsun)+log10(R(nR+1)/(Rstar*Rsun))*(real(i)-0.1)/real(j))
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
			P(i,j)%A=pi*(P(i,j)%R2**2-P(i,j)%R1**2)/real(nImPhi)
		enddo
	enddo

	vmax=0d0
	do i=1,nImR
	do j=1,nImPhi
		PP => P(i,j)
		trac%x=P(i,j)%R*cos(P(i,j)%Phi)
		trac%y=P(i,j)%R*sin(P(i,j)%Phi)
		trac%z=sqrt(R(nR+1)**2-P(i,j)%R**2)
		trac%edgeNr=2
		trac%onEdge=.true.
		call rotate(trac%x,trac%y,trac%z,0d0,1d0,0d0,inc*pi/180d0)
		ct=abs(trac%z)/R(nR+1)
		trac%i=nR
		do k=1,nTheta
			if(ct.lt.Theta(k).and.ct.ge.Theta(k+1)) then
				trac%j=k
			endif
		enddo

		trac%vx=0d0
		trac%vy=0d0
		trac%vz=-1d0
		call rotate(trac%vx,trac%vy,trac%vz,0d0,1d0,0d0,inc*pi/180d0)

		allocate(P(i,j)%vmin(nmol))
		allocate(P(i,j)%vmax(nmol))
		allocate(P(i,j)%npopmax(nmol))
		P(i,j)%vmin(1:nmol)=1d30
		P(i,j)%vmax(1:nmol)=-1d30
		P(i,j)%npopmax(1:nmol)=0
		call tracepath(trac,PP)
		do imol=1,nmol
			if((P(i,j)%vmax(imol)).gt.vmax) vmax=abs(P(i,j)%vmax(imol))
			if((-P(i,j)%vmin(imol)).gt.vmax) vmax=abs(P(i,j)%vmin(imol))
		enddo
	enddo
	enddo
	call output("Maximum velocity encountered: "//trim(dbl2string(vmax/1d5,'(f6.1)'))
     &			//" km/s")

	path2star%R=0d0
	path2star%Phi=0d0
	path2star%phi1=0d0
	path2star%phi2=2d0*pi
	path2star%R1=0d0
	path2star%R2=Rstar*Rsun
	path2star%A=pi*(path2star%R2**2-path2star%R1**2)

	trac%x=0d0
	trac%y=0d0
	trac%z=R(nR+1)
	trac%edgeNr=2
	trac%onEdge=.true.
	call rotate(trac%x,trac%y,trac%z,0d0,1d0,0d0,inc*pi/180d0)
	ct=abs(trac%z)/R(nR+1)
	trac%i=nR
	do k=1,nTheta
		if(ct.lt.Theta(k).and.ct.ge.Theta(k+1)) then
			trac%j=k
		endif
	enddo

	trac%vx=0d0
	trac%vy=0d0
	trac%vz=-1d0
	call rotate(trac%vx,trac%vy,trac%vz,0d0,1d0,0d0,inc*pi/180d0)

	allocate(path2star%vmin(nmol))
	allocate(path2star%vmax(nmol))
	allocate(path2star%npopmax(nmol))
	path2star%vmin(1:nmol)=1d30
	path2star%vmax(1:nmol)=-1d30
	path2star%npopmax(1:nmol)=0
	PP => path2star
	call tracepath(trac,PP)

	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine tracepath(trac,PP)
	use GlobalSetup
	IMPLICIT NONE
	type(Tracer) trac
	real*8 x,y,z,phi,d,vtot,taumin
	integer inext,jnext,ntrace,i,j,k,ipop,imol
	logical hitstar
	type(Path) PP

	PP%n=0
	PP%x=trac%x
	PP%y=trac%y
	PP%z=trac%z
	PP%vx=trac%vx
	PP%vy=trac%vy
	PP%vz=trac%vz

	allocate(PP%v(1000))
	allocate(PP%v1(1000))
	allocate(PP%v2(1000))
	allocate(PP%d(1000))
	allocate(PP%i(1000))
	allocate(PP%j(1000))

	taumin=0d0

1	continue

	call Trace2edge(trac,d,inext,jnext)

	if(PP%n.eq.0) then
		PP%n=PP%n+1
	else if(C(PP%i(k),PP%j(k))%dens.gt.1d-50.and.PP%i(k).gt.0) then
		PP%n=PP%n+1
	endif
	k=PP%n

	PP%d(k)=d

	PP%i(k)=trac%i
	PP%j(k)=trac%j
	
	x=trac%x
	y=trac%y
	z=trac%z

	PP%v1(k)=C(trac%i,trac%j)%v*(trac%vx*y-trac%vy*x)/sqrt(x**2+y**2)

	x=trac%x+trac%vx*d/2d0
	y=trac%y+trac%vy*d/2d0
	z=trac%z+trac%vz*d/2d0

	PP%v(k)=C(trac%i,trac%j)%v*(trac%vy*x-trac%vx*y)/sqrt(x**2+y**2)

c here I still have to add the turbulent velocity widening of the line
c add this for all species to get the absolute max and min velocity contributing.
	if(C(trac%i,trac%j)%dens.gt.1d-50) then
		do imol=1,nmol
			vtot=abs(PP%v(k))+3d0*C(trac%i,trac%j)%line_width(imol)
			if(vtot.gt.PP%vmax(imol).and.trac%i.gt.0.and.taumin.lt.tau_max) PP%vmax(imol)=vtot
			vtot=abs(PP%v(k))-3d0*C(trac%i,trac%j)%line_width(imol)
			if(vtot.lt.PP%vmin(imol).and.trac%i.gt.0.and.taumin.lt.tau_max) PP%vmin(imol)=vtot
			do ipop=Mol(imol)%nlevels,1,-1
				if(C(trac%i,trac%j)%npop(imol,ipop).gt.1d-150) exit
			enddo
			if(ipop.gt.PP%npopmax(imol)) PP%npopmax(imol)=ipop
		enddo
	endif


	trac%x=trac%x+trac%vx*d
	trac%y=trac%y+trac%vy*d
	trac%z=trac%z+trac%vz*d

	x=trac%x
	y=trac%y
	z=trac%z

	PP%v2(k)=C(trac%i,trac%j)%v*(trac%vx*y-trac%vy*x)/sqrt(x**2+y**2)

	if(trac%i.gt.0) taumin=taumin+minval(C(trac%i,trac%j)%kext(ilam1:ilam2))*d

	if(inext.gt.nR) then
		return
	endif
	if(inext.lt.0) then
c		p%hitstar=.true.
		return
	endif

	trac%i=inext
	trac%j=jnext

	goto 1

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


!c-----------------------------------------------------------------------
!c This subroutine determines the distance a photon has to travel to
!c the next border of the cell. The output variables are v, the distance
!c inext, the indicator for the next radial cell, and jnext, for the
!c next theta cell.
!c-----------------------------------------------------------------------
	subroutine Trace2edge(trac,v,inext,jnext)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Tracer) trac
	real*8 b,v,rr,R1,R2,T1,T2,vR1,vR2,vT1,vT2
	integer inext,jnext
	logical hitR1,hitR2,hitR,hitT1,hitT2,hitT,hitTsame

	rr=trac%x**2+trac%y**2+trac%z**2
	if(trac%i.ne.0) then
		R1=R(trac%i)**2
	else
		R1=(Rstar*Rsun)**2
	endif
	R2=R(trac%i+1)**2
	T1=Theta(trac%j)**2
	T2=Theta(trac%j+1)**2

	b=2d0*(trac%x*trac%vx+trac%y*trac%vy+trac%z*trac%vz)

	if(.not.trac%onEdge) then
		hitR1=hitR(trac,R1,rr,b,vR1)
		hitR2=hitR(trac,R2,rr,b,vR2)
		hitT1=hitT(trac,T1,rr,b,vT1)
		hitT2=hitT(trac,T2,rr,b,vT2)
	else
	if(trac%edgeNr.eq.1) then
		hitR1=.false.
		vR1=1d200
		hitR2=hitR(trac,R2,rr,b,vR2)
		hitT1=hitT(trac,T1,rr,b,vT1)
		hitT2=hitT(trac,T2,rr,b,vT2)
	else if(trac%edgeNr.eq.2) then
		hitR1=hitR(trac,R1,rr,b,vR1)
		hitR2=.true.
		vR2=-b
		hitT1=hitT(trac,T1,rr,b,vT1)
		hitT2=hitT(trac,T2,rr,b,vT2)
	else if(trac%edgeNr.eq.3) then
		hitR1=hitR(trac,R1,rr,b,vR1)
		hitR2=hitR(trac,R2,rr,b,vR2)
		hitT1=.false.
		vT1=1d200
		if(trac%j.ne.(nTheta)) then
			hitT2=hitT(trac,T2,rr,b,vT2)
		else
			vT2=-trac%z/trac%vz
			if(vT2.gt.0d0) then
				hitT2=.true.
			else
				hitT2=.false.
			endif
		endif
	else if(trac%edgeNr.eq.4) then
		hitR1=hitR(trac,R1,rr,b,vR1)
		hitR2=hitR(trac,R2,rr,b,vR2)
		if(trac%j.ne.1) then
			hitT1=hitT(trac,T1,rr,b,vT1)
		else
			hitT1=.false.
			vT1=1d200
		endif
		if(trac%j.ne.(nTheta)) then
			hitT2=hitTsame(trac,T2,rr,b,vT2)
		else
			hitT2=.false.
			vT2=1d200
		endif
	endif
	endif

	if(.not.hitR2) then
		print*,'Cannot hit outer boundary...'
		stop
	endif

	if(.not.hitR1.and..not.hitR2.and..not.hitT1.and..not.hitT2) then
		print*,'nothing to hit!'
		stop
	endif

	v=1d200
	if(hitR1.and.vR1.lt.v) then
		v=vR1
		inext=trac%i-1
		jnext=trac%j
		trac%edgeNr=2
	endif
	if(hitR2.and.vR2.lt.v) then
		v=vR2
		inext=trac%i+1
		jnext=trac%j
		trac%edgeNr=1
	endif
	if(hitT1.and.vT1.lt.v) then
		v=vT1
		inext=trac%i
		jnext=trac%j-1
		trac%edgeNr=4
	endif
	if(hitT2.and.vT2.lt.v) then
		v=vT2
		inext=trac%i
		jnext=trac%j+1
		trac%edgeNr=3
	endif


	if(jnext.eq.nTheta+1) then
		jnext=nTheta
		trac%edgeNr=4
	endif
	if(jnext.eq.0) then
		jnext=1
		trac%edgeNr=3
	endif
	
	trac%onEdge=.true.

	return
	end	
	
	logical function hitR(trac,Rad,rr,b,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Tracer) trac
	real*8 Rad,rr,b,cc,discr,vr1,vr2,v,q
	
	hitR=.false.
	v=1d200

	cc=rr-Rad
	discr=(b**2-4d0*cc)
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			q=-0.5d0*(b+discr)
		else
			q=-0.5d0*(b-discr)
		endif
		vr1=q
		vr2=cc/q
		if(vr1.gt.0d0) then
			v=vr1
			hitR=.true.
		endif
		if(vr2.gt.0d0.and.vr2.lt.v) then
			v=vr2
			hitR=.true.
		endif
	endif
	return
	end

	logical function hitT(trac,Thet,rr,b,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Tracer) trac
	real*8 Thet,rr,b,at,bt,ct,discr,vt1,vt2,v,q

	hitT=.false.
	v=1d200

	at=Thet-trac%vz*trac%vz
	bt=Thet*b-2d0*trac%z*trac%vz
	ct=Thet*rr-trac%z*trac%z
	discr=bt*bt-4d0*at*ct
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(bt.gt.0d0) then
			q=-0.5d0*(bt+discr)
		else
			q=-0.5d0*(bt-discr)
		endif
		vt1=q/at
		vt2=ct/q
		if(vt1.gt.0d0) then
			v=vt1
			hitT=.true.
		endif
		if(vt2.gt.0d0.and.vt2.lt.v) then
			v=vt2
			hitT=.true.
		endif
	endif
	return
	end

	logical function hitTsame(trac,Thet,rr,b,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Tracer) trac
	real*8 Thet,rr,b,at,bt,v

	hitTsame=.true.
	v=1d200

	bt=Thet*b-2d0*trac%z*trac%vz
	at=Thet-trac%vz**2
	v=-bt/at
	if(v.le.0d0) hitTsame=.false.

	return
	end

!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------
	