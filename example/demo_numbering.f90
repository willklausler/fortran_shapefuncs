program demo_numbering
  !! Print the node numbering of the quadratic tetrahedron and check that
  !! each shape function is 1 at its own node and 0 at the others
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  use shapefuncs
  implicit none
  type(shapefunc) :: s
  real(dp) :: func(10), derv(3,10), curv(3,3,10), err
  integer :: i

  s = shapefunc(cubature(CUB_TET, 2), 2)       ! 10-node tetrahedron
  call s%numbering()                           ! vertices 1-4, then edge midpoints 5-10

  err = 0
  do i = 1, s%nnodes
    call s%eval(s%coords(:,i), func, derv, curv)
    func(i) = func(i) - 1
    err = max(err, maxval(abs(func)))
  end do
  print '(a,es8.1)', 'max |N_j(x_i) - delta_ij| =', err  ! 0.0E+00
end program demo_numbering
