program shapefuncs_test
!! Validate shape functions

  use iso_fortran_env, only: ik => int32, rk => real64, stdout => output_unit
  use cubatures, only: cubature, elmtypes
  use shapefuncs, only: shapefunc, shaporders

  implicit none

  integer :: g, i, j
  integer :: ord
  integer :: maxo

  ! real(rk) ::

  character(1) :: order

  character(*), parameter :: fmt1 = "(A15,': ',$)"

  type(cubature) :: scheme
  type(shapefunc) :: shape

  write(*,"(A)") "1D elements"
  write(*,fmt1) "Line"

  j = maxval(maxloc(index(elmtypes, "LIN")))
  maxo = shaporders(j)
  do i = 1,maxo
    write(order,"(I1)") i

    call scheme%set("LIN", [i])

    ! Standard domain
    call shape%set(scheme, i)

    do g = 1,shape%points
      if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
        write(*,*)
        write(*,*) sum(shape%func)
        error stop "order "//order//" standard: partition of unity failed"
      end if

      if (.not.is_zero(sum(shape%derv(1,:,g)))) then
        write(*,*)
        write(*,*) sum(shape%derv(1,:,g))
        error stop "order "//order//" standard: derivatives failed"
      end if
    end do ! g

    ! Infinite domain
    call shape%set(scheme, i, infin=["INF"])

    do g = 1,shape%points
      if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
        write(*,*)
        write(*,*) sum(shape%func)
        error stop "order "//order//" infinite: partition of unity failed"
      end if

      if (.not.is_zero(sum(shape%derv(1,:,g)))) then
        write(*,*)
        write(*,*) sum(shape%derv(1,:,g))
        error stop "order "//order//" infinite: derivatives failed"
      end if
    end do ! g

  end do ! i
  write(*,"(A)") "passed"

contains

pure elemental function is_zero(r) result(z)
  real(rk), intent(in) :: r
  logical :: z
  z = abs(r) < 10.0_rk**(-12)
end function is_zero

end program shapefuncs_test