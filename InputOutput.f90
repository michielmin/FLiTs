!-----------------------------------------------------------------------
MODULE InOut
contains
  subroutine outputform(string, form)
    implicit none
    character(len=*), intent(in)           :: string
    character(len=*), intent(in), optional :: form

    if (present(form)) then
      write (*, form) trim(string)
    else
      write (*, '(a)') trim(string)
    end if

  end subroutine outputform

  !-----------------------------------------------------------------------

  subroutine output(string)
    implicit none
    character(len=*), intent(in) :: string

    write (*, '(a)') trim(string)

  end subroutine output

  !-----------------------------------------------------------------------

  subroutine ignorestar(un)
    implicit none
    integer, intent(in) :: un
    character           :: c
    integer             :: ios

    do
      read (un, '(a1)', iostat=ios) c
      if (ios /= 0) exit
      if (c /= '*') then
        backspace (unit=un)
        exit
      end if
    end do

  end subroutine ignorestar

  !-----------------------------------------------------------------------

  function int2string(i, form)
    implicit none
    integer, intent(in) :: i
    character(len=*), intent(in), optional :: form
    character(len=:), allocatable :: int2string
    character(len=range(i) + 2) :: tmp

    if (present(form)) then
      write (tmp, form) i
    else
      write (tmp, *) i
    end if

    allocate (character(len=len(trim(tmp))) :: int2string)
    int2string = trim(tmp)

  end function int2string

  !-----------------------------------------------------------------------

  function dbl2string(x, form)
    implicit none
    real(kind=8), intent(in)  :: x
    character(len=*), intent(in), optional :: form
    character(len=:), allocatable  :: dbl2string
    character(len=range(x) + 2)     :: tmp

    if (present(form)) then
      write (tmp, form) x
    else
      write (tmp, *) x
    end if

    allocate (character(len=len(trim(tmp))) :: dbl2string)
    dbl2string = trim(tmp)

  end function dbl2string

  !-----------------------------------------------------------------------
  subroutine clock(t, ut)
    ! clock both the CPU-time (t) and user time (ut)
    ! both values will be returned
    implicit none
    real(kind=8), intent(out) :: t, ut
    integer :: count, count_rate, count_max

    call CPU_TIME(t)
    call SYSTEM_CLOCK(count, count_rate, count_max)
    ut = real(count)/real(count_rate)
  end subroutine clock

  !-----------------------------------------------------------------------
  subroutine clock_write(t0, ut0, info)
    ! write out some time information (CPU-time and user time)
    ! t0, ut0 have to be initialised via CLOCK (see above)
    ! info string is optional and is written just before the time information
    implicit none
    real(kind=8), intent(in)           :: t0, ut0
    character(len=*), intent(in), optional :: info
    character(len=200) :: outstr
    real(kind=8) :: t1, ut1

    call clock(t1, ut1)

    write (outstr, '("time =",F10.2," s CPU-time=",F10.2," s (",F6.2," x speedup)")') &
      ut1 - ut0, t1 - t0, (t1 - t0)/(ut1 - ut0)

    if (present(info)) then
      outstr = trim(info)//" "//outstr
    end if
    call output(outstr)
    call output("")
    flush(6)
  end subroutine clock_write
end MODULE InOut
