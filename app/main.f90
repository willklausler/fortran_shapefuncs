program interpolation_convergence
  !! Interpolate exp((x + y + z)/2) on every reference element with shape
  !! functions of increasing order and print the maximum error at the points
  !! of a fine cubature, with the number of nodes: the table a finite element
  !! user needs to choose an element order for a smooth field.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures, only: cubature, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
  use shapefuncs, only: shapefunc
  implicit none

  integer, parameter :: elements(6) = [CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED]
  character(*), parameter :: names(6) = ['line ', 'tri  ', 'quad ', 'tet  ', 'hex  ', 'prism']
  integer, parameter :: pmax = 6
  integer, parameter :: degree = 12            ! evaluation points, independent of the order
  type(cubature) :: q
  type(shapefunc) :: s
  integer :: i, p

  write(*, '(a)') 'Maximum interpolation error (nodes) of exp((x + y + z)/2)'
  write(*, '(a5, 6a16)') 'order', names
  do p = 1, pmax
    write(*, '(i5)', advance='no') p
    do i = 1, size(elements)
      q = cubature(elements(i), degree)
      s = shapefunc(q, p)
      write(*, '(es11.2, " (", i3, ")")', advance='no') &
        maxval(abs(matmul(f(s%coords), s%func) - f(q%abscissae))), s%nnodes
    end do
    write(*, '(a)') ''
  end do

contains

pure function f(x) result(y)
  !! The interpolated field at the points `x(dim, n)`
  real(dp), intent(in) :: x(:,:)
  real(dp) :: y(size(x, 2))
  y = exp(sum(x, dim=1)/2)
end function f

end program interpolation_convergence
