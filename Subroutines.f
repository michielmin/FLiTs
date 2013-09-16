c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function Planck(T,lam)
	IMPLICIT NONE
	real*8 T,k,c,h,nu,lam,x
	k=1.3807d-16
	c=2.9979d10
	h=6.6261d-27
	nu=c/(lam*1d-4)
	x=h*nu/(k*T)
	Planck=(2d0*h*nu**3/c**2)/(exp(x)-1d0)
c	Planck=Planck*1e23

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotate(x,y,z,u,v,w,theta)
	IMPLICIT NONE
	real*8 x,y,z,u,v,w,yy(3),theta,inp
	real*8 cost,sint,u2,v2,w2
	cost=cos(theta)
	sint=sin(theta)
	u2=u*u
	v2=v*v
	w2=w*w

	inp=x*u+y*v+z*w
	yy(1)=u*inp
     & +(x*(v2+w2)-u*(v*y+w*z))*cost
     & +(v*z-w*y)*sint
	yy(2)=v*inp
     & +(y*(u2+w2)-v*(u*x+w*z))*cost
     & +(w*x-u*z)*sint
	yy(3)=w*inp
     & +(z*(u2+v2)-w*(u*x+v*y))*cost
     & +(u*y-v*x)*sint
	x=yy(1)
	y=yy(2)
	z=yy(3)
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine sort(x,n)
	IMPLICIT NONE
	integer n,i,j,imin
	real*8 x(n),min
	
	do j=1,n-1
	min=x(j)
	imin=j
	do i=j,n
		if(x(i).lt.min) then
			min=x(i)
			imin=i
		endif
	enddo
	min=x(j)
	x(j)=x(imin)
	x(imin)=min
	enddo
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje(i,n)
	IMPLICIT NONE
	integer i,n,f
	
	if(i.eq.1) call output("....................")
	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*real(i-1)/real(n).lt.real(f)
     &   .and.20d0*real(i+1)/real(n).gt.real(f)) then
		call outputform(".",'(a1,$)')
	endif
	
	if(i.eq.n) call output("")

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje_time(i,n,ii,nn,starttime)
	IMPLICIT NONE
	integer i,n,f
	integer ii,nn
	real*8 starttime,stoptime,xx
	character*20 dbl2string
	
	if(i.eq.1) then
		call cpu_time(stoptime)
		xx=100d0*real(ii)/real(nn)
		call output(trim(dbl2string((stoptime-starttime)/real(ii),'(f8.2)'))
     &			//" s per line. Approx " // 
     &			trim(dbl2string((stoptime-starttime)*(nn-ii)/real(ii),'(f8.2)'))
     &			//" s left. (" //
     &			trim(dbl2string(xx,'(f5.1)')) // " %)")
	endif
	
	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*real(i-1)/real(n).lt.real(f)
     &   .and.20d0*real(i+1)/real(n).gt.real(f)) then
		call cpu_time(stoptime)
		xx=100d0*real(ii)/real(nn)
		call output(trim(dbl2string((stoptime-starttime)/real(ii),'(f8.2)'))
     &			//" s per line. Approx " // 
     &			trim(dbl2string((stoptime-starttime)*(nn-ii)/real(ii),'(f8.2)'))
     &			//" s left. (" //
     &			trim(dbl2string(xx,'(f5.1)')) // " %)")
	endif

	return
	end

