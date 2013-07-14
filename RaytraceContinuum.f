	subroutine RaytraceContinuum()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam
	real*8 tot,TraceFlux,lam0
	
	open(unit=20,file='out.dat')
	lam0=lmin
	ilam=1
	do while(lam0.lt.lmax)
		tot=0d0
		do i=1,nImR
			do j=1,nImPhi
				tot=tot+TraceFlux(i,j,lam0,ilam)*P(i,j)%A
			enddo
		enddo
		write(20,*) lam0,tot
		lam0=lam0+lam0/rlines
		do while(lam0.gt.lam_cont(ilam+1))
			ilam=ilam+1
		enddo
	enddo
	close(unit=20)
	
	return
	end
	
	
	real*8 function TraceFlux(i,j,lam0,ilam)
	use GlobalSetup
	use Constants
	integer i,j,k,ilam
	real*8 lam0,tau,tot
	type(PathElement),pointer :: current
	
	current => P(i,j)%start
	tau=0d0
	tot=0d0
	do k=1,P(i,j)%n
		tot=tot+exp(-lam0)*exp(-tau)
		tau=tau+current%C%kext(ilam)*current%d
		if(tau.gt.10d0) exit
	enddo
	
	TraceFlux=tot
	
	return
	end
	