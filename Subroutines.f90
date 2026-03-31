!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

real(kind=8) function Planck(T, lam)
  implicit none
  real(kind=8) :: T, k, c, h, nu, lam, x
  k = 1.3807d-16
  c = 2.9979d10
  h = 6.6261d-27
  nu = c/(lam*1d-4)
  x = h*nu/(k*T)
  if (x > 40d0) then
    Planck = (2d0*h*nu**3/c**2)*exp(-x)
  else if (x < 0.1) then
    Planck = (2d0*h*nu**3/c**2)*(-0.5d0+1.d0/x+x/12.d0-x**3/720.d0)
  else
    Planck = (2d0*h*nu**3/c**2)/(exp(x)-1d0)
  endif
!  Planck=Planck*1e23

  return
end function Planck

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine rotate(x, y, z, u, v, w, theta)
  implicit none
  real(kind=8) :: x, y, z, u, v, w, yy(3), theta, inp
  real(kind=8) :: cost, sint, u2, v2, w2
  cost = cos(theta)
  sint = sin(theta)
  u2 = u*u
  v2 = v*v
  w2 = w*w

  inp = x*u+y*v+z*w
  yy(1) = u*inp &
       + (x*(v2+w2)-u*(v*y+w*z))*cost &
       + (v*z-w*y)*sint
  yy(2) = v*inp &
       + (y*(u2+w2)-v*(u*x+w*z))*cost &
       + (w*x-u*z)*sint
  yy(3) = w*inp &
       + (z*(u2+v2)-w*(u*x+v*y))*cost &
       + (u*y-v*x)*sint
  x = yy(1)
  y = yy(2)
  z = yy(3)
  return
end subroutine rotate

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine sort(x, n)
  implicit none
  integer :: n, i, j, imin
  real(kind=8) :: x(n), min

  do j = 1, n-1
    min = x(j)
    imin = j
    do i = j, n
      if (x(i) < min) then
        min = x(i)
        imin = i
      endif
    enddo
    min = x(j)
    x(j) = x(imin)
    x(imin) = min
  enddo

  return
end subroutine sort

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

real(kind=8) function ran1(idum)
  implicit none
  integer :: idum, IA, IM, IQ, IR, NTAB, NDIV
  real(kind=4) :: AM, EPS, RNMX
  parameter (IA=16807, IM=2147483647, AM=1./IM, IQ=127773, IR=2836, &
             NTAB=32, NDIV=1+(IM-1)/NTAB, EPS=1.2e-7, RNMX=1.-EPS)
  integer :: j, k
  integer, save :: iv(NTAB), iy
  data iv/NTAB*0/, iy/0/
!$omp threadprivate(iv,iy)
  if (idum <= 0 .or. iy == 0) then
    idum = max(-idum, 1)
    do j = NTAB+8, 1, -1
      k = idum/IQ
      idum = IA*(idum-k*IQ)-IR*k
      if (idum < 0) idum = idum+IM
      if (j <= NTAB) iv(j) = idum
    enddo
    iy = iv(1)
  endif
  k = idum/IQ
  idum = IA*(idum-k*IQ)-IR*k
  if (idum < 0) idum = idum+IM
  j = 1+iy/NDIV
  iy = iv(j)
  iv(j) = idum
  ran1 = min(AM*iy, RNMX)
  return
end function ran1


subroutine tellertje(i, n)
  use InOut
  implicit none
  integer :: i, n, f

  if (i == 1) call output("....................")
  f = int(20d0*dble(i)/dble(n))

  if (20d0*real(i-1)/real(n) < real(f) &
      .and. 20d0*real(i+1)/real(n) > real(f)) then
    call outputform(".", '(a1,$)')
    flush(6)
  endif

  if (i == n) call output("")

  return
end subroutine tellertje

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine tellertje_time(i, n, ii, nn, starttime)
  use InOut
  implicit none
  integer :: i, n, f
  integer :: ii, nn
  integer(kind=4) :: counts, count_rate, count_max
  real(kind=8) :: starttime, stoptime, xx

  if (i == 1) then
    !call cpu_time(stoptime)
    call SYSTEM_CLOCK(counts, count_rate, count_max)
    stoptime = dble(counts)/dble(count_rate)
    xx = 100d0*real(ii)/real(nn)
    call output(trim(dbl2string((stoptime-starttime)/real(i), '(f8.3)')) &
         //" s per line. Approx " // &
         trim(dbl2string((stoptime-starttime)*(nn-ii)/real(ii), '(f10.2)')) &
         //" s left. (" // &
         trim(dbl2string(xx, '(f5.1)')) // " %)")
  endif

  f = int(20d0*dble(i)/dble(n))

  if (20d0*real(i-1)/real(n) < real(f) &
      .and. 20d0*real(i+1)/real(n) > real(f)) then
    !call cpu_time(stoptime)
    call SYSTEM_CLOCK(counts, count_rate, count_max)
    stoptime = dble(counts)/dble(count_rate)
    xx = 100d0*real(ii)/real(nn)
    call output(trim(dbl2string((stoptime-starttime)/real(i), '(f8.3)')) &
         //" s per line. Approx " // &
         trim(dbl2string((stoptime-starttime)*(nn-ii)/real(ii), '(f10.2)')) &
         //" s left. (" // &
         trim(dbl2string(xx, '(f5.1)')) // " %)")
    flush(6)
  endif

  return
end subroutine tellertje_time

