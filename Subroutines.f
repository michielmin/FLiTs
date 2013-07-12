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

