program shapefuncs_example
!! Use [[shapefuncs]] the way a finite element code does:
!!
!! 1. inspect the node numbering of an element,
!! 2. compute the area of a curved quadrilateral through an isoparametric
!!    map,
!! 3. integrate a decaying field over an unbounded domain with a mapped
!!    infinite element,
!! 4. keep one set of shape functions per element type in a table indexed
!!    by `CUB_*`.

  use cubatures, only: cubature, rk, CUB_LIN, CUB_QUA, CUB_TET, CUB_WED
  use shapefuncs, only: shapefunc, SHP_INF, SHP_CHP

  implicit none

  call inspect_numbering()
  call curved_area()
  call infinite_element()
  call function_table()

contains

!***********************************************************************

subroutine inspect_numbering()
!! Print the nodes of the 10-node tetrahedron

  type(shapefunc) :: s

  write(*,"(/,A)") "1. Nodes of the quadratic tetrahedron"
  s = shapefunc(cubature(CUB_TET, 2), 2)
  call s%numbering()

end subroutine inspect_numbering

!***********************************************************************

subroutine curved_area()
!! Area of a quarter annulus \(1 \le r \le 2\), \(0 \le \theta \le \pi/2\),
!! meshed by one 9-node quadrilateral whose nodes lie on the exact
!! geometry, compared to the exact area \(3\pi/4\)

  real(rk), parameter :: pi = acos(-1.0_rk)
  type(cubature) :: q
  type(shapefunc) :: s
  real(rk) :: xn(2,9), jac(2,2), area, r, t
  integer :: g, n

  write(*,"(/,A)") "2. Area of a quarter annulus, one 9-node quadrilateral"

  q = cubature(CUB_QUA, 6)
  s = shapefunc(q, 2)

  ! Place each node on the annulus: xi_1 -> radius, xi_2 -> angle
  do n = 1, s%nnodes
    r = 1.5_rk + 0.5_rk*s%coords(1,n)
    t = pi/4*(1 + s%coords(2,n))
    xn(:,n) = r*[cos(t), sin(t)]
  end do

  ! Area = sum_g w_g det(dx/dxi); dx/dxi = sum_n x_n dN_n/dxi
  area = 0
  do g = 1, q%npoints
    jac = matmul(xn, transpose(s%derv(:,:,g)))
    area = area + (jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1))*q%weights(g)
  end do

  write(*,"(A,F12.8)") "Isoparametric: ", area
  write(*,"(A,F12.8)") "Exact:         ", 3*pi/4

end subroutine curved_area

!***********************************************************************

subroutine infinite_element()
!! \(\int_1^\infty u^2\,dx\) for \(u = 1/x\), whose exact value is 1, on
!! one infinite line element with nodes at \(x = 1\) and \(x = 2\).
!!
!! Linear mapping functions give \(x(\xi) = 2/(1-\xi)\), so \(\xi \to 1\)
!! is \(x \to \infty\). Linear chopped functions interpolate \(u\) from
!! its nodal values; they represent \(1/x\) exactly.

  type(cubature) :: q
  type(shapefunc) :: map, field
  real(rk), parameter :: xn(2) = [1.0_rk, 2.0_rk]   ! Nodal coordinates
  real(rk) :: un(2), u, dxdxi, total
  integer :: g

  write(*,"(/,A)") "3. Integral of 1/x**2 from 1 to infinity, one infinite element"

  q = cubature(CUB_LIN, 2)
  map = shapefunc(q, 1, [SHP_INF])
  field = shapefunc(q, 1, [SHP_CHP])
  un = 1/xn

  total = 0
  do g = 1, q%npoints
    dxdxi = dot_product(map%derv(1,:,g), xn)
    u = dot_product(field%func(:,g), un)
    total = total + u**2*dxdxi*q%weights(g)
  end do

  write(*,"(A,F12.8)") "Infinite element: ", total
  write(*,"(A,F12.8)") "Exact:            ", 1.0_rk

end subroutine infinite_element

!***********************************************************************

subroutine function_table()
!! One set of quadratic shape functions per element type, with a matching
!! cubature, indexed by the element constants

  type(shapefunc) :: table(6)
  integer :: elm

  write(*,"(/,A)") "4. Quadratic shape functions for every element type"

  do elm = CUB_LIN, CUB_WED
    table(elm) = shapefunc(cubature(elm, 4), 2)
  end do

  associate (s => table(CUB_WED))
    call s%summary()
  end associate

  write(*,"(A,6I4)") "Nodes per element:", table%nnodes
  write(*,"(A,6I4)") "Points per element:", table%npoints

end subroutine function_table

!***********************************************************************

end program shapefuncs_example
