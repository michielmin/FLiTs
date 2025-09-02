	subroutine SetupPaths()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ip,jp,i,j,k,ir,nRreduce,ilam,imol,nImPhi_max,nPhiMin,nPhiMax,iint,it
	integer,allocatable :: startphi(:)
	real*8 ct,res_inc,x,y,rr,tt,ph
	real*8,allocatable :: imR(:),imPhi(:,:)
	type(Tracer) trac,trac_count
	type(Path),pointer :: PP
	real*8 incfact,x11,x12,x21,x22,y11,y12,y21,y22,r1,r2,r3,s
	real*8 ComputeIncFact,maxRjump
	real*8 imx,imy,imz
	integer,allocatable :: t_node(:,:),t_neighbor(:,:)
	integer matri,icount,ncount
	real*8,allocatable :: xy(:,:)
	real*8 matrix(2,2)
	interface
	  subroutine output(string)
	  IMPLICIT NONE
	  character string*(*)
	  end
	  character*20 function int2string(i,form)
	  IMPLICIT NONE
	  integer i
	  character,intent(in),optional :: form*(*)
	  end
	  character*20 function dbl2string(x,form)
	  IMPLICIT NONE
	  real*8 x
	  character,intent(in),optional :: form*(*)
	  end
	end interface

	ilam1=1
	ilam2=nlam
	do ilam=1,nlam
		if(lam_cont(ilam).lt.lmin) ilam1=ilam
		if(lam_cont(nlam+1-ilam).gt.lmax) ilam2=nlam+1-ilam
	enddo
	
c increase the resolution in velocity by this factor

	if(accuracy.lt.0) then
		nrReduce=8
		res_inc=1d0
		nPhiMin=2
		nPhiMax=10
	else if(accuracy.eq.0) then
		nrReduce=4
		res_inc=1d0
		nPhiMin=15
		nPhiMax=60
	else if(accuracy.eq.1) then
		nrReduce=2
		res_inc=1d0
		nPhiMin=30
		nPhiMax=75
	else if(accuracy.eq.2) then
		nrReduce=1
		res_inc=2d0
		nPhiMin=30
		nPhiMax=90
	else if(accuracy.eq.3) then
		nrReduce=1
		res_inc=4d0
		nPhiMin=45
		nPhiMax=90
	else if(accuracy.eq.4) then
		nrReduce=1
		res_inc=8d0
		nPhiMin=45
		nPhiMax=180
	else
		nrReduce=1
		res_inc=10d0
		nPhiMin=120
		nPhiMax=270
	endif		
		
	call output("==================================================================")
	call output("Setup up the paths for raytracing")

	do while(nR/nRreduce.lt.40.and.nRreduce.gt.1)
		nRreduce=nRreduce-1
	enddo
	
	ir=0
	do i=nTheta,1,-1
		call trace2tau1(i,rr)
		if(rr.gt.R_sphere(1)) then
			ir=ir+3
		endif
	enddo
	do i=1,nR,nRreduce
		ir=ir+1
		if(cylindrical) then
			ir=ir+1
		endif
	enddo

	nImR=ir+int(abs(sin(inc*pi/180d0))*(C(1,nTheta)%v/vresolution)*res_inc/2d0)+(nTheta-1)*2+100

	allocate(imR(nImR))

	ir=0
	do i=nTheta,1,-1
		call trace2tau1(i,rr)
		if(rr.gt.R_sphere(1)) then
			ir=ir+1
			imR(ir)=rr
			ir=ir+1
			imR(ir)=rr*cos(inc*pi/180d0-(pi/2d0-theta_av(nR-1,i)))
			ir=ir+1
			imR(ir)=rr*cos(inc*pi/180d0+(pi/2d0-theta_av(nR-1,i)))
		endif
	enddo
	do i=1,nR,nRreduce
		ir=ir+1
		imR(ir)=R_av_sphere(i)
		if(cylindrical) then
			ir=ir+1
			imR(ir)=R_av(i)
		endif
	enddo

	call sort(imR,ir)
	do while(imR(ir).gt.R_sphere(nr+1))
		ir=ir-1
	enddo

	do i=2,nTheta
		if(cylindrical) then
			rr=R(1)/sin(theta_av(nR-1,i))
		else
			rr=R_sphere(1)
		endif
		ir=ir+1
		imR(ir)=rr*cos(inc*pi/180d0-(pi/2d0-theta_av(nR-1,i)))/ComputeIncFact(rr)
		imR(ir)=rr*cos(inc*pi/180d0-(pi/2d0-theta_av(nR-1,i)))/ComputeIncFact(imR(ir))
		imR(ir)=rr*cos(inc*pi/180d0-(pi/2d0-theta_av(nR-1,i)))/ComputeIncFact(imR(ir))
		ir=ir+1
		imR(ir)=rr*cos(inc*pi/180d0+(pi/2d0-theta_av(nR-1,i)))/ComputeIncFact(rr)
		imR(ir)=rr*cos(inc*pi/180d0+(pi/2d0-theta_av(nR-1,i)))/ComputeIncFact(imR(ir))
		imR(ir)=rr*cos(inc*pi/180d0+(pi/2d0-theta_av(nR-1,i)))/ComputeIncFact(imR(ir))
	enddo

	j=(nImR-ir)
	do i=1,j
		ir=ir+1
		imR(ir)=10d0**(log10(Rstar*Rsun)+log10(R_sphere(nR+1)/(Rstar*Rsun))*(real(i)-0.1)/real(j))
		imR(ir)=10d0**(log10(R(nR))+log10(R_sphere(nR+1)/R(nR))*(real(i)-0.1)/real(j))
	enddo

	imR=abs(imR)

	call sort(imR,nImR)

	open(unit=20,file='imagegrid.out',RECL=1000)
	do i=1,nImR
		write(20,*) imR(i)/AU,R_sphere(nR+1)/AU
	enddo
	close(unit=20)

	allocate(nImPhi(nImR))
	allocate(startphi(nImR+1))

	nImPhi_max=1
	startphi(1)=1

	npoints_temp=0
	do i=1,nImR
		do k=1,nR-1
			if(imR(i).gt.R(k).and.imR(i).lt.R(k+1)) exit
		enddo
		nImPhi(i)=abs(sin(inc*pi/180d0))*(C(k,nTheta)%v/vresolution)*res_inc
		if(nImPhi(i).lt.nPhiMin) nImPhi(i)=nPhiMin
		if(nImPhi(i).gt.nPhiMax) nImPhi(i)=nPhiMax
		if(startphi(i).eq.1) then
			nImPhi(i)=nImPhi(i)+1
			startphi(i+1)=0
		else
			startphi(i+1)=1
		endif
		if(nImPhi(i).gt.nImPhi_max) nImPhi_max=nImPhi(i)
		npoints_temp=npoints_temp+nImPhi(i)
	enddo
	ngrids=5

	allocate(imPhi(nImR,nImPhi_max))
	
	do i=1,nImR
		if(startphi(i).eq.1) then
			do j=1,nImPhi(i)
				ImPhi(i,j)=pi*(real(j)-0.5)/real(nImPhi(i))
			enddo
		else
			do j=1,nImPhi(i)
				ImPhi(i,j)=pi*real(j-1)/real(nImPhi(i)-1)
			enddo
		endif
	enddo
	
	allocate(xy(2,npoints_temp))
	allocate(npoints(ngrids))
	matri=2*npoints_temp-3
	allocate(P(ngrids,matri))
	allocate(t_node(3,matri))
	allocate(t_neighbor(3,matri))
	k=0
	do i=1,ngrids
		do j=1,npoints_temp
			ir=ran1(idum)*real(nR)
			it=ran1(idum)*real(nTheta)
			rr=ran1(idum)
			rr=sqrt(R(ir)**2*rr+R(ir+1)**2*(1d0-rr))
			tt=ran1(idum)
			tt=Theta(ir,it)*tt+Theta(ir,it+1)*(1d0-tt)
			tt=acos(tt)
			ph=ran1(idum)*pi*2d0
			imx=rr*cos(ph)*sin(tt)
			imy=rr*sin(ph)*sin(tt)
			imz=rr*cos(tt)
			! CHR: my guess is that this rotate is done to optimiye the path grid 
			! however, because of this rotate the outer edge of the disk is not sampled on 
			! the far side of the disk. To be backward compatible still do the rotation if 
			! imagecube is not used
			if (.not.imagecube) then 
				call rotate(imx,imy,imz,0d0,1d0,0d0,-inc*pi/180d0)
			endif
			xy(1,j)=imx
			xy(2,j)=imy
		enddo
		matrix(1,1)=1d0
		matrix(1,2)=0d0
		matrix(2,1)=0d0
		matrix(2,2)=1d0
		call dtris2_lmap (npoints_temp, xy, matrix, npoints(i), t_node, t_neighbor )
		do j=1,npoints(i)
			P(i,j)%x=(xy(1,t_node(1,j))+xy(1,t_node(2,j))+xy(1,t_node(3,j)))/3d0
			P(i,j)%y=(xy(2,t_node(1,j))+xy(2,t_node(2,j))+xy(2,t_node(3,j)))/3d0
			P(i,j)%y=abs(P(i,j)%y)
			r1=sqrt((xy(1,t_node(1,j))-xy(1,t_node(2,j)))**2+(xy(2,t_node(1,j))-xy(2,t_node(2,j)))**2)
			r2=sqrt((xy(1,t_node(1,j))-xy(1,t_node(3,j)))**2+(xy(2,t_node(1,j))-xy(2,t_node(3,j)))**2)
			r3=sqrt((xy(1,t_node(3,j))-xy(1,t_node(2,j)))**2+(xy(2,t_node(3,j))-xy(2,t_node(2,j)))**2)
			s=(r1+r2+r3)/2d0
			P(i,j)%A=sqrt(s*(s-r1)*(s-r2)*(s-r3))
		enddo
		k=k+npoints(i)
	enddo
	ncount=k
	k=k/ngrids

	call output("Number of image gridpoints: "//trim(int2string(k,'(i7)')))	

	vmax=0d0
	icount=0
	do i=1,ngrids
	do j=1,npoints(i)
		icount=icount+1
		call tellertje(icount,ncount)
		PP => P(i,j)
		trac%x=P(i,j)%x
		trac%y=P(i,j)%y

		rr=trac%x**2+trac%y**2

		trac%z=sqrt(R_sphere(nR+1)**2-rr)
		trac%edgeNr=2
		trac%onEdge=.true.
		call rotate(trac%x,trac%y,trac%z,0d0,1d0,0d0,inc*pi/180d0)
		ct=abs(trac%z)/R_sphere(nR+1)
		trac%i=nR
		do k=0,nTheta
			if(ct.lt.Theta(trac%i,k).and.ct.ge.Theta(trac%i,k+1)) then
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

		trac_count=trac
		call tracepath(trac_count,PP,.true.)
		call tracepath(trac,PP,.false.)

		do imol=1,nmol
			if((P(i,j)%vmax(imol)).gt.vmax) vmax=abs(P(i,j)%vmax(imol))
			if((-P(i,j)%vmin(imol)).gt.vmax) vmax=abs(P(i,j)%vmin(imol))
		enddo
	enddo
	enddo
	call output("Maximum velocity encountered: "//trim(dbl2string(vmax/1d5,'(f6.1)'))
     &			//" km/s")

	path2star%R=Rstar*Rsun/4d0
	path2star%Phi=pi/4d0
	path2star%phi1=0d0
	path2star%phi2=2d0*pi
	path2star%R1=0d0
	path2star%R2=Rstar*Rsun
	path2star%A=pi*(path2star%R2**2-path2star%R1**2)

	trac%x=path2star%R*cos(path2star%Phi)
	trac%y=path2star%R*sin(path2star%Phi)
	trac%z=sqrt(R_sphere(nR+1)**2-path2star%R**2)

	trac%edgeNr=2
	trac%onEdge=.true.
	call rotate(trac%x,trac%y,trac%z,0d0,1d0,0d0,inc*pi/180d0)
	ct=abs(trac%z)/R_sphere(nR+1)
	trac%i=nR
	do k=0,nTheta
		if(ct.lt.Theta(trac%i,k).and.ct.ge.Theta(trac%i,k+1)) then
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

	trac_count=trac
	call tracepath(trac_count,PP,.true.)
	call tracepath(trac,PP,.false.)

	deallocate(startphi)

	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine tracepath(trac,PP,onlycount)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Tracer) trac
	real*8 x,y,z,phi,d,vtot,taumin,v1,v2,v
	integer inext,jnext,ntrace,i,j,k,ipop,imol,kk,nk
	logical hitstar,onlycount
	type(Path) PP

	PP%vx=trac%vx
	PP%vy=trac%vy
	PP%vz=trac%vz

	if(.not.onlycount) then
		allocate(PP%v(PP%n))
		allocate(PP%v1(PP%n))
		allocate(PP%v2(PP%n))
		allocate(PP%d(PP%n))
		allocate(PP%i(PP%n))
		allocate(PP%j(PP%n))
	endif
	PP%n=1

	taumin=0d0

1	continue

	call Trace2edge(trac,d,inext,jnext)

	k=PP%n
	
	x=trac%x
	y=trac%y
	z=trac%z

	v1=C(trac%i,trac%j)%v*(trac%vx*y-trac%vy*x)/sqrt(x**2+y**2)
	if(.not.onlycount) then
		PP%v1(k)=v1
	endif

	x=trac%x+trac%vx*d
	y=trac%y+trac%vy*d
	z=trac%z+trac%vz*d

	v2=C(trac%i,trac%j)%v*(trac%vx*y-trac%vy*x)/sqrt(x**2+y**2)
	if(.not.onlycount) then
		PP%v2(k)=v2
	endif
	
	nk=1
	do imol=1,nmol
		i=abs(v1-v2)*5d0/C(trac%i,trac%j)%line_width(imol)
		if(i.gt.nk) nk=i
	enddo
	i=abs(v1-v2)/(vres_profile)+1
	if(nk.gt.i) nk=i

	x=trac%x
	y=trac%y
	z=trac%z

	do kk=1,nk
		k=PP%n

		if(.not.onlycount) then
			PP%d(k)=d/real(nk)
			PP%i(k)=trac%i
			PP%j(k)=trac%j
		endif

		x=trac%x+(real(kk)-0.5)*trac%vx*d/real(nk)
		y=trac%y+(real(kk)-0.5)*trac%vy*d/real(nk)
		z=trac%z+(real(kk)-0.5)*trac%vz*d/real(nk)

		call checkcell(x,y,z,trac%i,trac%j)

		if(.not.onlycount) then
			v=sqrt(G*Mstar*Msun*(x*x+y*y)/((x**2+y**2+z**2)**(3d0/2d0)))
c			PP%v(k)=C(trac%i,trac%j)%v*(trac%vy*x-trac%vx*y)/sqrt(x**2+y**2)
			PP%v(k)=v*(trac%vy*x-trac%vx*y)/sqrt(x**2+y**2)

			do imol=1,nmol
				if(C(trac%i,trac%j)%N(imol).gt.1d-50
     &		.and.trac%j.gt.0.and.trac%i.lt.nR.and.trac%i.gt.0) then
					vtot=abs(PP%v(k))+3d0*C(trac%i,trac%j)%line_width(imol)
					if(vtot.gt.PP%vmax(imol).and.trac%i.gt.0.and.taumin.lt.tau_max) PP%vmax(imol)=vtot
					vtot=abs(PP%v(k))-3d0*C(trac%i,trac%j)%line_width(imol)
					if(vtot.lt.PP%vmin(imol).and.trac%i.gt.0.and.taumin.lt.tau_max) PP%vmin(imol)=vtot
					ipop=C(trac%i,trac%j)%npopmax(imol)
					if(ipop.gt.PP%npopmax(imol)) PP%npopmax(imol)=ipop
				endif
			enddo
		endif

		if(C(trac%i,trac%j)%dens.gt.1d-50
     &		.and.trac%i.gt.0.and.trac%i.lt.nR.and.trac%j.gt.0) then
			PP%n=PP%n+1
		endif
	enddo

	trac%x=trac%x+trac%vx*d
	trac%y=trac%y+trac%vy*d
	trac%z=trac%z+trac%vz*d

	if(trac%i.gt.0) taumin=taumin+minval(C(trac%i,trac%j)%kext(ilam1:ilam2))*d

	if(inext.gt.nR) then
		return
	endif
	if(inext.lt.0) then
		return
	endif

	trac%i=inext
	trac%j=jnext

	goto 1

	PP%n=PP%n-1

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
	real*8 b,v,rr,R1,R2,T1,T2,vR1,vR2,vT1,vT2,rcyl,bcyl,R1s,R2s,ct
	integer inext,jnext,k
	logical hitR1,hitR2,hitR,hitT1,hitT2,hitT,hitTsame

	rr=trac%x**2+trac%y**2+trac%z**2
	rcyl=trac%x**2+trac%y**2

	if(trac%i.gt.0) then
		R1=R(trac%i)**2
		R1s=R_sphere(trac%i)**2
	else
		R1=(Rstar*Rsun)**2
		R1s=(Rstar*Rsun)**2
	endif
	R2=R(trac%i+1)**2
	R2s=R_sphere(trac%i+1)**2
	T1=Theta(trac%i,trac%j)**2
	T2=Theta(trac%i,trac%j+1)**2

	b=2d0*(trac%x*trac%vx+trac%y*trac%vy+trac%z*trac%vz)
	bcyl=2d0*(trac%x*trac%vx+trac%y*trac%vy)/sqrt(trac%vx**2+trac%vy**2)

	if(.not.trac%onEdge) then
		if(cylindrical.and.trac%i.gt.0.and.trac%j.gt.0) then
			hitR1=hitR(trac,R1,rcyl,bcyl,vR1)
			vR1=vR1/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR1=hitR(trac,R1s,rr,b,vR1)
		endif
		if(cylindrical.and.trac%i.lt.nR.and.trac%j.gt.0) then
			hitR2=hitR(trac,R2,rcyl,bcyl,vR2)
			vR2=vR2/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR2=hitR(trac,R2s,rr,b,vR2)
		endif
		hitT1=hitT(trac,T1,rr,b,vT1)
		hitT2=hitT(trac,T2,rr,b,vT2)
	else
	if(trac%edgeNr.eq.1) then
		hitR1=.false.
		vR1=1d200
		if(cylindrical.and.trac%i.lt.nR.and.trac%j.gt.0) then
			hitR2=hitR(trac,R2,rcyl,bcyl,vR2)
			vR2=vR2/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR2=hitR(trac,R2s,rr,b,vR2)
		endif
		hitT1=hitT(trac,T1,rr,b,vT1)
		hitT2=hitT(trac,T2,rr,b,vT2)
	else if(trac%edgeNr.eq.2) then
		if(cylindrical.and.trac%i.gt.0.and.trac%j.gt.0) then
			hitR1=hitR(trac,R1,rcyl,bcyl,vR1)
			vR1=vR1/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR1=hitR(trac,R1s,rr,b,vR1)
		endif
		hitR2=.true.
		if(cylindrical.and.trac%i.lt.nR.and.trac%j.gt.0) then
			vR2=-bcyl
			vR2=vR2/sqrt(trac%vx**2+trac%vy**2)
		else
			vR2=-b
		endif
		hitT1=hitT(trac,T1,rr,b,vT1)
		hitT2=hitT(trac,T2,rr,b,vT2)
	else if(trac%edgeNr.eq.3) then
		if(cylindrical.and.trac%i.gt.0.and.trac%j.gt.0) then
			hitR1=hitR(trac,R1,rcyl,bcyl,vR1)
			vR1=vR1/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR1=hitR(trac,R1s,rr,b,vR1)
		endif
		if(cylindrical.and.trac%i.lt.nR.and.trac%j.gt.0) then
			hitR2=hitR(trac,R2,rcyl,bcyl,vR2)
			vR2=vR2/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR2=hitR(trac,R2s,rr,b,vR2)
		endif
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
		if(cylindrical.and.trac%i.gt.0.and.trac%j.gt.0) then
			hitR1=hitR(trac,R1,rcyl,bcyl,vR1)
			vR1=vR1/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR1=hitR(trac,R1s,rr,b,vR1)
		endif
		if(cylindrical.and.trac%i.lt.nR.and.trac%j.gt.0) then
			hitR2=hitR(trac,R2,rcyl,bcyl,vR2)
			vR2=vR2/sqrt(trac%vx**2+trac%vy**2)
		else
			hitR2=hitR(trac,R2s,rr,b,vR2)
		endif
		if(trac%j.ne.0) then
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
		if(inext.ge.0) then
			ct=abs(trac%z+v*trac%vz)/sqrt((trac%x+v*trac%vx)**2+(trac%y+v*trac%vy)**2+(trac%z+v*trac%vz)**2)
			do k=0,nTheta
				if(ct.lt.Theta(inext,k).and.ct.ge.Theta(inext,k+1)) then
					jnext=k
				endif
			enddo
		endif
		trac%edgeNr=2
	endif
	if(hitR2.and.vR2.lt.v) then
		v=vR2
		inext=trac%i+1
		jnext=trac%j
		if(inext.le.nR) then
			ct=abs(trac%z+v*trac%vz)/sqrt((trac%x+v*trac%vx)**2+(trac%y+v*trac%vy)**2+(trac%z+v*trac%vz)**2)
			do k=0,nTheta
				if(ct.lt.Theta(inext,k).and.ct.ge.Theta(inext,k+1)) then
					jnext=k
				endif
			enddo
		endif
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
		if(trac%edgeNr.eq.3) trac%edgeNr=4
	endif
	if(jnext.lt.0) then
		jnext=0
		if(trac%edgeNr.eq.4) trac%edgeNr=3
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

	subroutine trace2tau1(j,rr)
	use GlobalSetup
	IMPLICIT NONE
	real*8 rr,tau,tau0
	integer j,i
	
	tau=0d0
	do i=1,nR
		if(cylindrical.and.j.ne.0) then
			tau0=maxval(C(i,j)%kext)*(R(i+1)-R(i))/sin(theta_av(i,j))
		else
			tau0=maxval(C(i,j)%kext)*(R_sphere(i+1)-R_sphere(i))
		endif
		if(tau+tau0.gt.1d0) then
			if(cylindrical) then
				rr=R(i)+(R(i+1)-R(i))*(1d0-tau)/tau0
				rr=rr/sin(theta_av(i,j))
				return
			else
				rr=R_sphere(i)+(R_sphere(i+1)-R_sphere(i))*(1d0-tau)/tau0
				return
			endif
		endif
		tau=tau+tau0
	enddo
	rr=0d0

	return
	end	
	
!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

	real*8 function ComputeIncFact(R0)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 R0

	ComputeIncFact=(cos(inc*pi/180d0)+(1d0-cos(inc*pi/180d0))
     &		*((R0-R_sphere(1))/(R_sphere(nR+1)-R_sphere(1)))**2)
	if(ComputeIncFact.lt.cos(inc*pi/180d0)) ComputeIncFact=cos(inc*pi/180d0)
		
	return
	end
	
	
	subroutine checkcell(x,y,z,i,j)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 rr,rho,thet,x,y,z
	integer i,j
	
	return
	
	rr=sqrt(x*x+y*y+z*z)
	rho=sqrt(x*x+y*y)
	thet=abs(z/rr)

	if(thet.gt.Theta(i,j)) then
		print*,'theta too big',i,j
		stop
	endif
	
	if(thet.lt.Theta(i,j+1)) then
		print*,'theta too small',i,j
		print*,thet,Theta(i,j+1)
		stop
	endif
	
	if(i.eq.0.or.j.eq.0) then
		if(rr.lt.R_sphere(i)) then
			print*,'R too small',i,j
			print*,rr/AU,R_sphere(i+1)/AU
			stop
		endif
	else
		if(rho.lt.R(i)) then
			print*,'Rho too small',i,j
			print*,rho/AU,R(i+1)/AU
			stop
		endif
	endif
		
	if(i.eq.nr.or.j.eq.0) then
		if(rr.gt.R_sphere(i+1)) then
			print*,'R too big',i,j
			print*,rr/AU,R_sphere(i+1)/AU
			stop
		endif
	else
		if(rho.gt.R(i+1)) then
			print*,'Rho too big',i,j
			print*,rho/AU,R(i+1)/AU
			stop
		endif
	endif
		
	return
	end
	
