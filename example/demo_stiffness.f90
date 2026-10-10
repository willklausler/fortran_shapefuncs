program demo_stiffness
  !! The element loop of a finite element code: the stiffness matrix of the
  !! Laplacian on a bilinear quadrilateral, from the isoparametric map and
  !! the physical gradients dN/dx = J**(-T) dN/dxi
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  use shapefuncs
  implicit none
  real(dp), parameter :: xn(2,4) = reshape([0._dp, 0._dp, 1._dp, 0._dp, &
                                            1._dp, 1._dp, 0._dp, 1._dp], [2, 4])  ! unit square
  type(cubature) :: q
  type(shapefunc) :: s
  real(dp) :: k(4,4), jac(2,2), jinv(2,2), detj, dndx(2,4)
  integer :: g

  q = cubature(CUB_QUA, 2)                     ! degree 2p for order p = 1
  s = shapefunc(q, 1)
  k = 0
  do g = 1, q%npoints
    jac  = matmul(xn, transpose(s%derv(:,:,g)))
    detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
    jinv = reshape([jac(2,2), -jac(2,1), -jac(1,2), jac(1,1)], [2, 2])/detj
    dndx = matmul(transpose(jinv), s%derv(:,:,g))
    k = k + matmul(transpose(dndx), dndx)*detj*q%weights(g)
  end do

  print '(4f9.5)', k(1,:)                      ! 0.66667 -0.16667 -0.33333 -0.16667
  print '(a,es8.1)', 'max row sum:', maxval(abs(sum(k, dim=2)))  ! 0 up to rounding
end program demo_stiffness
