subroutine ComputeLTE()
  use GlobalSetup
  use Constants
  use InOut
  implicit none
  integer :: i, j, k, imol, nLTE
  real(kind=8) :: UT, T

  nLTE = 0
  do imol = 1, nmol
    if (Mol(imol)%LTE) then
      call output("Computing LTE level populations for "//trim(Mol(imol)%name))
      nLTE = nLTE + 1
    end if
  end do

  if (nLTE > 0) then
    do i = 0, nR
      call tellertje(i + 1, nR + 1)
      do j = 1, nTheta
        T = C(i, j)%Tgas
        if (T .lt. 3d0) T = 3d0
        do imol = 1, nmol
          if (Mol(imol)%LTE) then
            UT = 0d0
            do k = 1, Mol(imol)%nlevels
              UT = UT + Mol(imol)%g(k)*exp(-Mol(imol)%E(k)/T)
            end do

            do k = 1, Mol(imol)%nlevels
              C(i, j)%npop(imol)%N(k) = Mol(imol)%g(k)*exp(-Mol(imol)%E(k)/T)/UT
            end do
          end if
        end do
      end do
    end do
  end if
  return
end

