	subroutine RaytraceLines()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,ilam,k,ilines,vmult,iv,nv
	real*8 lam,T,Planck,wl1,wl2,v,flux0,starttime,stoptime
	real*8,allocatable :: flux(:)
	type(Path),pointer :: PP
	
	call output("==================================================================")
	call output("Tracing the lines")

	call cpu_time(starttime)

	open(unit=20,file='out.dat')

	nv=int(vmax/vresolution)

	allocate(flux(-nv:nv))

	lam=lmin
	ilam=1
	do while(lam.gt.lam_cont(ilam+1))
		ilam=ilam+1
	enddo
	do i=0,nR
		do j=1,nTheta
			allocate(C(i,j)%profile(-nv:nv))
		enddo
	enddo

	do i=0,nR
		do j=1,nTheta
			do k=-nv,nv
				C(i,j)%profile(k)=((real(k))/5d0)**2
			enddo
			C(i,j)%profile=exp(-C(i,j)%profile)
		enddo
	enddo
	do ilines=1,100
		call tellertje(ilines,100)
		flux=0d0

		do i=0,nR
			do j=1,nTheta
				C(i,j)%kext_l=C(i,j)%kext(ilam)
				C(i,j)%albedo_l=C(i,j)%albedo(ilam)
				C(i,j)%BB_l=BB(ilam,C(i,j)%iT)
			enddo
		enddo

		do i=1,nImR
!$OMP PARALLEL IF(.false.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(flux0,iv,j,PP,vmult)
!$OMP& SHARED(i,ilam,P,flux,nImR,nImPhi,nv)
!$OMP DO
			do j=1,nImPhi
				PP => P(i,j)
				do vmult=-1,1,2
					do iv=-nv,nv
						if(real(iv*vmult).gt.(PP%vmax+15d0).or.real(iv*vmult).lt.(PP%vmin-15d0)) then
							flux0=PP%flux_cont(ilam)
						else
							call TraceFluxLines(PP,flux0,iv,nv,vmult)
						endif
						flux(iv)=flux(iv)+flux0*PP%A/2d0
					enddo
				enddo
			enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		enddo
		do i=-nv,nv
			write(20,*) i,flux(i)
		enddo
	enddo
	close(unit=20)
	
	call cpu_time(stoptime)

	call output("Time used for the lines: "//trim(dbl2string(stoptime-starttime,'(f7.1)'))
     &			//" s")


	return
	end
	
	
	subroutine TraceFluxLines(p0,flux,ii,nv,vmult)
	use GlobalSetup
	use Constants
	integer i,j,k,vmult,iv,ii,nv
	real*8 tau,exptau,flux,fact,profile,xx
	type(Path) p0
	type(Cell),pointer :: CC

	fact=1d0
	flux=0d0

	do k=1,p0%n
		i=p0%i(k)
		if(i.ne.0) then
			j=p0%j(k)
			CC => C(i,j)
			iv=int(p0%v(k))*vmult
			jj=ii+iv
			if(jj.lt.-nv) jj=-nv
			if(jj.gt.nv) jj=nv
			profile=CC%profile(jj)
			tau=CC%kext_l*(1d0+profile)*p0%d(k)
			exptau=exp(-tau)
			xx=CC%BB_l*(1d0-CC%albedo_l)
			flux=flux+xx*(1d0-exptau)*fact
			fact=fact*exptau
			if(fact.lt.1d-6) exit
		endif
	enddo
	
	return
	end
	
