subroutine RaytraceContinuum()
  use GlobalSetup
  use Constants
  use InOut
  implicit none
  integer          :: i, j, ilam, k, ilamstart, ilamend
  double precision :: flux, lam0, T, Planck, wl1, wl2
  logical          :: doit
	double precision :: time, utime

  call output("==================================================================")
  call clock(time,utime)
	!call output("Raytracing the continuum")

! BB is not used and this produces a floating overflow, already for T=1
!  allocate(BB(nlam,MAXT))
!  do ilam=1,nlam
!    do i=1,MAXT
!      T=real(i)
!      BB(ilam,i)=Planck(T,lam_cont(ilam))
!    end do
!  end do

  do i = 1, ngrids
    do j = 1, npoints(i)
      allocate(P(i,j)%flux_cont(nlam))
    end do
  end do
  allocate(path2star%flux_cont(nlam))

  do ilamstart = 1, nlam-1
    if (lam_cont(ilamstart+1) > lmin) exit
  end do

  do ilamend = ilamstart, nlam
    if (lam_cont(ilamend) > lmax) exit
  end do
  ilamend = min(ilamend, nlam) ! case of last point

  call output("Raytracing the continuum for " // trim(int2string(ilamend-ilamstart,'(i5)')) // " wavelengths")	
  do ilam = ilamstart, ilamend
    call tellertje(ilam-ilamstart+1, ilamend-ilamstart)
    flux = 0d0
    do i = 1, ngrids
!$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(+:flux) &
!$OMP& PRIVATE(j) &
!$OMP& SHARED(i, ilam, npoints, P, ngrids)
      do j = 1, npoints(i)
        call TraceFluxCont(P(i,j), ilam, P(i,j)%flux_cont(ilam))
        flux = flux + P(i,j)%flux_cont(ilam)*P(i,j)%A/real(ngrids)
      end do
!$OMP END PARALLEL DO
    end do
    call Trace2StarCont(path2star, ilam, path2star%flux_cont(ilam))
    flux = flux + path2star%flux_cont(ilam)*path2star%A
  end do
	call clock_write(time,utime,"RaytraceContinuum:")
  

end subroutine RaytraceContinuum

!-----------------------------------------------------------------------

subroutine TraceFluxCont(p0, ilam, flux)
  use GlobalSetup
  use Constants
  implicit none
  integer          :: i, j, k, ilam
  double precision :: tau, exptau, flux, fact, S
  type(Path)       :: p0

  fact = 1d0
  flux = 0d0
  do k = 1, p0%n
    i   = p0%i(k)
    j   = p0%j(k)
    tau = C(i,j)%kext(ilam)*p0%d(k)
!   S=BB(ilam,C(i,j)%iT)*(1d0-C(i,j)%albedo(ilam))
!   S=S+C(i,j)%LRF(ilam)*C(i,j)%albedo(ilam)
!   Using the source function for now.
    S = C(i,j)%S(ilam)
    if (tau > 1d-6) then
      exptau = exp(-tau)
    else
      exptau = 1d0 - tau
    end if
    flux = flux + S*(1d0-exptau)*fact
    fact = fact*exptau
    if (fact < 1d-6) exit
  end do

end subroutine TraceFluxCont

!-----------------------------------------------------------------------

subroutine Trace2StarCont(p0, ilam, flux)
  use GlobalSetup
  use Constants
  implicit none
  integer          :: i, j, k, ilam
  double precision :: tau, flux
  type(Path)       :: p0

  tau = 0d0
  do k = 1, p0%n
    i = p0%i(k)
    if (i == 0) exit
    j   = p0%j(k)
    tau = tau + C(i,j)%kext(ilam)*p0%d(k)
  end do

  flux = Fstar(ilam)*exp(-tau)

end subroutine Trace2StarCont
