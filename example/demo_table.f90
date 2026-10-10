program demo_table
  !! Keep one set of quadratic shape functions per element type, with a
  !! matching cubature, in a table indexed by the element constants
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  use shapefuncs
  implicit none
  type(shapefunc) :: table(6)
  integer :: elm

  do elm = CUB_LIN, CUB_WED                    ! CUB_* are 1 to 6
    table(elm) = shapefunc(cubature(elm, 4), 2)
  end do

  call table(CUB_WED)%summary()                ! WED, order 2, 18 nodes, 18 points
  print '(a,6i4)', 'nodes: ', table%nnodes     ! 3   6   9  10  27  18
  print '(a,6i4)', 'points:', table%npoints    ! 3   6   9  27  27  18
end program demo_table
