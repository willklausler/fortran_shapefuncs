program shapefuncs_test
!! Validate shape functions

  use iso_fortran_env, only: rk => real64, stdout => output_unit
  use cubatures, only: cubature, elmtypes
  use shapefuncs, only: shapefunc, shaporders

  implicit none

  integer :: g, i, j, k, l, m, n
  integer :: ord
  integer :: maxo

  character(1) :: order

  character(*), parameter :: fmt1 = "(A15,': ',$)"
  character(3), parameter :: infi(2) = ["FIN", "INF"]

  type(cubature) :: scheme
  type(shapefunc) :: shape

  write(*,"(A)") "1D elements"
  write(*,fmt1) "Line"

  j = maxval(maxloc(index(elmtypes, "LIN")))
  maxo = shaporders(j)
  do i = 1,maxo
    write(order,"(I1)") i

    call scheme%set("LIN", [i])

    do k = 1,2
      call shape%set(scheme, i, infin=[infi(k)])

      do g = 1,shape%points
        if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
          write(*,*)
          write(*,*) sum(shape%func)
          error stop "order "//order//"-"//infi(k)//": partition of unity failed"
        end if

        if (.not.is_zero(sum(shape%derv(1,:,g)))) then
          write(*,*)
          write(*,*) sum(shape%derv(1,:,g))
          error stop "order "//order//"-"//infi(k)//": derivatives failed"
        end if

        if (.not.is_zero(sum(shape%curv(1,1,:,g)))) then
          write(*,*)
          write(*,*) sum(shape%curv(1,1,:,g))
          error stop "order "//order//"-"//infi(k)//": second derivatives failed"
        end if
      end do ! g
    end do ! k
  end do ! i
  write(*,"(A)") "passed"

  write(*,"(/,A)") "2D elements"

  write(*,fmt1) "Triangle"
  j = maxval(maxloc(index(elmtypes, "TRI")))
  maxo = shaporders(j)
  do i = 1,maxo
    write(order,"(I1)") i

    call scheme%set("TRI", [i])

    ! Standard domain
    call shape%set(scheme, i)

    do g = 1,shape%points
      if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
        write(*,*)
        write(*,*) sum(shape%func)
        error stop "order "//order//": partition of unity failed"
      end if

      do j = 1,2
        if (.not.is_zero(sum(shape%derv(j,:,g)))) then
          write(*,*)
          write(*,*) sum(shape%derv(j,:,g))
          error stop "order "//order//": derivatives failed"
        end if
      end do ! j

      do j = 1,2
        do k = 1,2
          if (.not.is_zero(sum(shape%curv(j,1,:,g)))) then
            write(*,*)
            write(*,*) sum(shape%curv(j,1,:,g))
            error stop "order "//order//": second derivatives failed"
          end if
        end do ! k
      end do ! j
    end do ! g
  end do ! i
  write(*,"(A)") "passed"

  write(*,fmt1) "Quadrilateral"
  j = maxval(maxloc(index(elmtypes, "QUA")))
  maxo = shaporders(j)
  do i = 1,maxo
    write(order,"(I1)") i

    call scheme%set("QUA", [i])

    do k = 1,2
      do l = 1,2

        ! Standard domain
        call shape%set(scheme, i, infin=[infi(k), infi(l)])

        do g = 1,shape%points
          if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
            write(*,*)
            write(*,*) sum(shape%func)
            error stop "order "//order//"-"//infi(k)//"-"//infi(l)// &
                       ": partition of unity failed"
          end if

          do j = 1,2
            if (.not.is_zero(sum(shape%derv(j,:,g)))) then
              write(*,*)
              write(*,*) sum(shape%derv(j,:,g))
              error stop "order "//order//"-"//infi(k)//"-"//infi(l)// &
                         ": derivatives failed"
            end if
          end do ! j

          do j = 1,2
            do m = 1,2
              if (.not.is_zero(sum(shape%curv(j,m,:,g)))) then
                write(*,*)
                write(*,*) sum(shape%curv(j,m,:,g))
                error stop "order "//order//"-"//infi(k)//"-"//infi(l)// &
                           ": second derivatives failed"
              end if
            end do ! m
          end do ! j
        end do ! g
      end do ! l
    end do ! k
  end do ! i
  write(*,"(A)") "passed"

  write(*,"(/,A)") "3D elements"

  write(*,fmt1) "Tetrahedron"
  j = maxval(maxloc(index(elmtypes, "TET")))
  maxo = shaporders(j)
  do i = 1,maxo
    write(order,"(I1)") i

    call scheme%set("TET", [i])

    ! Standard domain
    call shape%set(scheme, i)

    do g = 1,shape%points
      if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
        write(*,*)
        write(*,*) sum(shape%func)
        error stop "order "//order//": partition of unity failed"
      end if

      do j = 1,3
        if (.not.is_zero(sum(shape%derv(j,:,g)))) then
          write(*,*)
          write(*,*) sum(shape%derv(j,:,g))
          error stop "order "//order//": derivatives failed"
        end if
      end do ! j

      do j = 1,3
        do k = 1,3
          if (.not.is_zero(sum(shape%curv(j,k,:,g)))) then
            write(*,*)
            write(*,*) sum(shape%curv(j,k,:,g))
            error stop "order "//order//": second derivatives failed"
          end if
        end do ! k
      end do ! j
    end do ! g
  end do ! i
  write(*,"(A)") "passed"

  write(*,fmt1) "Hexahedron"
  j = maxval(maxloc(index(elmtypes, "HEX")))
  maxo = shaporders(j)
  do i = 1,maxo
    write(order,"(I1)") i

    call scheme%set("HEX", [i])

    do k = 1,2
      do l = 1,2
        do m = 1,2

          ! Standard domain
          call shape%set(scheme, i, infin=[infi(k), infi(l), infi(m)])

          do g = 1,shape%points
            if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
              write(*,*)
              write(*,*) sum(shape%func)
              error stop "order "//order//"-"//infi(k)//"-"//infi(l)//"-"//infi(m)// &
                         ": partition of unity failed"
            end if

            do j = 1,3
              if (.not.is_zero(sum(shape%derv(j,:,g)))) then
                write(*,*)
                write(*,*) sum(shape%derv(j,:,g))
                error stop "order "//order//"-"//infi(k)//"-"//infi(l)//"-"//infi(m)// &
                           ": derivatives failed"
              end if
            end do ! j

            do j = 1,3
              do n = 1,3
                if (.not.is_zero(sum(shape%curv(j,n,:,g)))) then
                  write(*,*)
                  write(*,*) sum(shape%curv(j,n,:,g))
                  error stop "order "//order//"-"//infi(k)//"-"//infi(l)//"-"//infi(m)// &
                            ": second derivatives failed"
                end if
              end do ! n
            end do ! j
          end do ! g
        end do ! m
      end do ! l
    end do ! k
  end do ! i
  write(*,"(A)") "passed"

  write(*,fmt1) "Prism"
  j = maxval(maxloc(index(elmtypes, "WEJ")))
  maxo = shaporders(j)
  do i = 1,maxo
    write(order,"(I1)") i

    call scheme%set("WEJ", [i])

    ! Standard domain
    call shape%set(scheme, i)

    do g = 1,shape%points
      if (.not.is_zero(sum(shape%func(:,g)) - 1.0_rk)) then
        write(*,*)
        write(*,*) sum(shape%func)
        error stop "order "//order//": partition of unity failed"
      end if

      do j = 1,3
        if (.not.is_zero(sum(shape%derv(j,:,g)))) then
          write(*,*)
          write(*,*) sum(shape%derv(j,:,g))
          error stop "order "//order//": derivatives failed"
        end if
      end do ! j

      do j = 1,3
        do k = 1,3
          if (.not.is_zero(sum(shape%curv(j,k,:,g)))) then
            write(*,*)
            write(*,*) sum(shape%curv(j,k,:,g))
            error stop "order "//order//": second derivatives failed"
          end if
        end do ! k
      end do ! j
    end do ! g
  end do ! i
  write(*,"(A)") "passed"

contains

pure elemental function is_zero(r) result(z)
  real(rk), intent(in) :: r
  logical :: z
  z = abs(r) < 10.0_rk**(-8)
end function is_zero

end program shapefuncs_test