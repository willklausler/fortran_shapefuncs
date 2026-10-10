program demo_infinite
  !! Integral of u**2 from 1 to infinity for u = 1/x, which is 1, on one
  !! infinite line element with nodes at x = 1 and x = 2. Linear mapping
  !! functions give x(xi) = 2/(1 - xi), so xi -> 1 is x -> infinity; linear
  !! chopped functions interpolate u and represent 1/x exactly.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  use shapefuncs
  implicit none
  real(dp), parameter :: xn(2) = [1._dp, 2._dp]   ! nodal coordinates
  type(cubature) :: q
  type(shapefunc) :: map, field
  real(dp) :: un(2), u, dxdxi, total
  integer :: g

  q = cubature(CUB_LIN, 2)
  map = shapefunc(q, 1, [SHP_INF])                ! geometry
  field = shapefunc(q, 1, [SHP_CHP])              ! solution, vanishes at infinity
  un = 1/xn

  total = 0
  do g = 1, q%npoints
    dxdxi = dot_product(map%derv(1,:,g), xn)
    u = dot_product(field%func(:,g), un)
    total = total + u**2*dxdxi*q%weights(g)
  end do

  print '(a,f12.8)', 'infinite element:', total   ! 1.00000000
end program demo_infinite
