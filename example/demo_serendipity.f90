program demo_serendipity
  !! Use the 20-node serendipity hexahedron in the node ordering of a mesh
  !! format that lists the edge nodes by their first vertex, as Gmsh does.
  !! The permutation for `nodes` comes from matching reference coordinates.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  use shapefuncs
  implicit none
  integer, parameter :: edges(2,12) = reshape([1,2, 1,4, 1,5, 2,3, 2,6, 3,4, &
                                               3,7, 4,8, 5,6, 5,8, 6,7, 7,8], [2, 12])
  type(cubature) :: q
  type(shapefunc) :: s
  real(dp) :: xmesh(3,20)
  integer :: perm(20), k

  q = cubature(CUB_HEX, 4)
  s = shapefunc(q, 2, family=SHP_SERENDIPITY)

  xmesh(:,1:8) = s%coords(:,1:8)               ! same vertices, then the mesh's edges
  do k = 1, 12
    xmesh(:,8+k) = (s%coords(:,edges(1,k)) + s%coords(:,edges(2,k)))/2
  end do

  do k = 1, 20                                 ! mesh node k is node perm(k) of s
    perm(k) = minloc(sum(abs(s%coords - spread(xmesh(:,k), 2, 20)), dim=1), 1)
  end do
  print '(a,20i3)', 'nodes =', perm            ! 1 ... 8  9 12 17 10 18 11 19 20 13 16 14 15

  s = shapefunc(q, 2, family=SHP_SERENDIPITY, nodes=perm)
  print '(a,l2)', 'same order:', all(abs(s%coords - xmesh) < 1e-12_dp)  ! T
end program demo_serendipity
