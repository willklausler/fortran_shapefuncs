program demo_isoparametric
  !! Area of the quarter annulus 1 <= r <= 2, 0 <= theta <= pi/2, meshed by
  !! one 9-node quadrilateral whose nodes lie on the exact geometry
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  use shapefuncs
  implicit none
  real(dp), parameter :: pi = acos(-1._dp)
  type(cubature) :: q
  type(shapefunc) :: s
  real(dp) :: xn(2,9), jac(2,2), area, r, t
  integer :: g, n

  q = cubature(CUB_QUA, 6)
  s = shapefunc(q, 2)

  do n = 1, s%nnodes                           ! xi_1 -> radius, xi_2 -> angle
    r = 1.5_dp + 0.5_dp*s%coords(1,n)
    t = pi/4*(1 + s%coords(2,n))
    xn(:,n) = r*[cos(t), sin(t)]
  end do

  area = 0                                     ! sum_g w_g det(dx/dxi)
  do g = 1, q%npoints
    jac = matmul(xn, transpose(s%derv(:,:,g))) ! dx/dxi = sum_n x_n dN_n/dxi
    area = area + (jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1))*q%weights(g)
  end do

  print '(a,f12.8)', 'isoparametric:', area    ! 2.32842712
  print '(a,f12.8)', 'exact 3*pi/4: ', 3*pi/4  ! 2.35619449
end program demo_isoparametric
