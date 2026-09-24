program shapefuncs_test
!! Unit tests for [[shapefuncs]].
!!
!! For every element type, a sweep of orders and every combination of
!! finite, infinite and chopped directions, checks that the shape functions
!!
!! - are consistent (`is_valid`) and have the expected number of nodes,
!! - are 1 at their own node and 0 at the others,
!! - reproduce every function of their space exactly, with its first and
!!   second derivatives,
!! - have first and second derivatives that agree with finite differences
!!   and symmetric second derivatives,
!! - sum to 1 (finite and mapping functions) or vanish at infinity
!!   (chopped functions).
!!
!! The same checks run on the serendipity quadrilaterals and hexahedra.
!! Also checks the classical node numbering of orders 1 to 3, closed forms
!! from the literature, integration with a cubature, caller-defined node
!! ordering, the constructor and `set` interfaces, reuse, `destroy`, and
!! output.

  use cubatures, only: cubature, rk, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
  use shapefuncs

  implicit none

  real(rk), parameter :: tol = 1.0e-10_rk   !! Relative tolerance on exact identities
  real(rk), parameter :: fdtol = 1.0e-6_rk  !! Relative tolerance on finite differences
  integer, parameter :: fam(3) = [SHP_FIN, SHP_INF, SHP_CHP]
  integer :: nfail = 0                      !! Number of failed checks
  integer :: p, i, j, k

  call section("Line")
  do p = 0, 8
    do i = 1, 3
      if (i > 1 .and. p == 0) cycle
      call check_element(CUB_LIN, p, [fam(i)])
    end do
  end do

  call section("Quadrilateral")
  do p = 0, 6
    do j = 1, 3
      do i = 1, 3
        if ((i > 1 .or. j > 1) .and. p == 0) cycle
        call check_element(CUB_QUA, p, [fam(i), fam(j)])
      end do
    end do
  end do

  call section("Hexahedron")
  do p = 0, 4
    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          if ((i > 1 .or. j > 1 .or. k > 1) .and. (p == 0 .or. p > 2)) cycle
          call check_element(CUB_HEX, p, [fam(i), fam(j), fam(k)])
        end do
      end do
    end do
  end do

  call section("Triangle")
  do p = 1, 8
    call check_element(CUB_TRI, p, [SHP_FIN])
  end do

  call section("Tetrahedron")
  do p = 1, 6
    call check_element(CUB_TET, p, [SHP_FIN])
  end do

  call section("Wedge")
  do p = 1, 5
    do i = 1, 3
      call check_element(CUB_WED, p, [SHP_FIN, SHP_FIN, fam(i)])
    end do
  end do

  call section("Serendipity")
  do p = 1, 3
    call check_element(CUB_QUA, p, [SHP_FIN], SHP_SERENDIPITY)
    call check_element(CUB_HEX, p, [SHP_FIN], SHP_SERENDIPITY)
  end do

  call section("Node counts")
  call expect(nodes(CUB_LIN, 4) == 5,    "LIN 4")
  call expect(nodes(CUB_QUA, 4) == 25,   "QUA 4")
  call expect(nodes(CUB_HEX, 3) == 64,   "HEX 3")
  call expect(nodes(CUB_TRI, 4) == 15,   "TRI 4")
  call expect(nodes(CUB_TET, 4) == 35,   "TET 4")
  call expect(nodes(CUB_WED, 3) == 40,   "WED 3")
  call expect(nodes(CUB_WED, 4) == 75,   "WED 4")
  call expect(nodes(CUB_QUA, 2, SHP_SERENDIPITY) == 8,  "QUA 2 serendipity")
  call expect(nodes(CUB_QUA, 3, SHP_SERENDIPITY) == 12, "QUA 3 serendipity")
  call expect(nodes(CUB_HEX, 2, SHP_SERENDIPITY) == 20, "HEX 2 serendipity")
  call expect(nodes(CUB_HEX, 3, SHP_SERENDIPITY) == 32, "HEX 3 serendipity")

  call section("Classical numbering")
  call test_numbering()

  call section("Closed forms")
  call test_closed_forms()

  call section("Node ordering")
  call test_ordering()

  call section("Integration")
  call test_integration()

  call section("Interface")
  call test_interface()

  call section("Output")
  call test_output()

  write(*,*)
  if (nfail > 0) then
    write(*,"(I0,A)") nfail, " check(s) failed"
    error stop 1
  end if
  write(*,"(A)") "All tests passed"

contains

!***********************************************************************

subroutine section(name)
!! Start a group of checks
  character(*), intent(in) :: name
  write(*,"(A)") name
end subroutine section

!***********************************************************************

subroutine expect(cond, msg)
!! Record a failed check
  logical, intent(in) :: cond
  character(*), intent(in) :: msg
  if (cond) return
  nfail = nfail + 1
  write(*,"(2X,'FAIL: ',A)") msg
end subroutine expect

!***********************************************************************

logical function close(a, b, rtol)
!! True if the arrays agree to `rtol` relative to the largest magnitude, at least 1
  real(rk), intent(in) :: a(:), b(:), rtol
  close = maxval(abs(a - b)) <= rtol*max(1.0_rk, maxval(abs(a)), maxval(abs(b)))
end function close

!***********************************************************************

integer function nodes(elm, order, family)
!! Number of nodes of an element
  integer, intent(in) :: elm, order
  integer, intent(in), optional :: family
  type(shapefunc) :: s
  s = shapefunc(cubature(elm, 1), order, family=family)
  nodes = s%nnodes
end function nodes

!***********************************************************************

integer function expected_nodes(elm, p, family)
!! Closed-form number of nodes
  integer, intent(in) :: elm, p, family
  if (family == SHP_SERENDIPITY) then
    expected_nodes = merge(4*p, 8 + 12*(p - 1), elm == CUB_QUA)
    return
  end if
  select case (elm)
  case (CUB_LIN); expected_nodes = p + 1
  case (CUB_QUA); expected_nodes = (p + 1)**2
  case (CUB_HEX); expected_nodes = (p + 1)**3
  case (CUB_TRI); expected_nodes = (p + 1)*(p + 2)/2
  case (CUB_TET); expected_nodes = (p + 1)*(p + 2)*(p + 3)/6
  case (CUB_WED); expected_nodes = (p + 1)**2*(p + 2)/2
  case default;   expected_nodes = -1
  end select
end function expected_nodes

!***********************************************************************

subroutine check_element(elm, p, infin, family)
!! Run every per-element check on one element, order, infinitude and family

  integer, intent(in) :: elm, p, infin(:)
  integer, intent(in), optional :: family

  type(cubature) :: q
  type(shapefunc) :: s
  character(40) :: label

  ! Degree 2p integrates the mass matrix of affine elements
  q = cubature(elm, 2*p)
  s = shapefunc(q, p, infin, family)

  write(label,"(I0,' p=',I0,' fam=',I0,' inf=',*(I0,:,','))") elm, p, s%family, infin

  call expect(s%is_valid(), trim(label)//" is_valid")
  call expect(s%nnodes == expected_nodes(elm, p, s%family), trim(label)//" node count")
  call expect(s%npoints == q%npoints, trim(label)//" point count")
  call check_kronecker(s, label)
  call check_reproduction(s, q, label)
  call check_finite_differences(s, q, label)
  call check_sums(s, q, label)

end subroutine check_element

!***********************************************************************

subroutine check_kronecker(s, label)
!! N_i(x_j) = delta_ij

  type(shapefunc), intent(in) :: s
  character(*), intent(in) :: label

  integer :: j
  real(rk) :: f(s%nnodes), d(s%dim,s%nnodes), c(s%dim,s%dim,s%nnodes), e(s%nnodes)
  logical :: ok

  ok = .true.
  do j = 1, s%nnodes
    call s%eval(s%coords(:,j), f, d, c)
    e = 0
    e(j) = 1
    ok = ok .and. close(f, e, tol)
  end do
  call expect(ok, trim(label)//" Kronecker delta at nodes")

end subroutine check_kronecker

!***********************************************************************

subroutine check_reproduction(s, q, label)
!! Interpolating any function of the element's space reproduces it and
!! its first and second derivatives at every cubature point.
!!
!! The space is spanned by products of one function per direction: the
!! monomials \(\xi^e\) along finite directions, \(\xi^e/(1-\xi)\) along
!! mapped infinite directions and \((1-\xi)\xi^e\) along chopped
!! directions, with \(e \le p\) per direction and total degree \(\le p\)
!! over the simplex directions. Serendipity elements span the monomials of
!! superlinear degree \(\le p\).

  type(shapefunc), intent(in) :: s
  type(cubature), intent(in) :: q
  character(*), intent(in) :: label

  integer :: e(3), g, i, j, k, p, ntri
  real(rk) :: tv, tg(3), th(3,3), val, grad(3), hess(3,3)
  logical :: ok

  p = s%order
  ntri = 0
  if (s%elm == CUB_TRI .or. s%elm == CUB_TET .or. s%elm == CUB_WED) ntri = min(s%dim, 2)
  if (s%elm == CUB_TET) ntri = 3

  ok = .true.
  e = 0
  do k = 0, merge(p, 0, s%dim > 2)
    do j = 0, merge(p, 0, s%dim > 1)
      do i = 0, p
        e = [i, j, k]
        if (ntri > 0) then
          if (sum(e(1:ntri)) > p) cycle
        end if
        if (s%family == SHP_SERENDIPITY) then
          if (sum(e, mask=e >= 2) > p) cycle
        end if
        do g = 1, s%npoints
          ! Interpolant
          val = 0
          grad = 0
          hess = 0
          block
            integer :: n
            real(rk) :: gdum(3), hdum(3,3)
            do n = 1, s%nnodes
              call test_function(s, e, s%coords(:,n), tv, gdum, hdum)
              val = val + tv*s%func(n,g)
              grad(1:s%dim) = grad(1:s%dim) + tv*s%derv(:,n,g)
              hess(1:s%dim,1:s%dim) = hess(1:s%dim,1:s%dim) + tv*s%curv(:,:,n,g)
            end do
          end block
          ! Exact
          call test_function(s, e, q%abscissae(:,g), tv, tg, th)
          ok = ok .and. close([val], [tv], tol) .and. close(grad(1:s%dim), tg(1:s%dim), tol) &
                  .and. close(reshape(hess(1:s%dim,1:s%dim), [s%dim**2]), &
                              reshape(th(1:s%dim,1:s%dim), [s%dim**2]), tol)
        end do
      end do
    end do
  end do
  call expect(ok, trim(label)//" polynomial reproduction")

end subroutine check_reproduction

!***********************************************************************

subroutine test_function(s, e, x, v, g, h)
!! Product test function of exponents `e` at `x`, with gradient and Hessian

  type(shapefunc), intent(in) :: s
  integer, intent(in) :: e(3)
  real(rk), intent(in) :: x(:)
  real(rk), intent(out) :: v, g(3), h(3,3)

  real(rk) :: u(3), du(3), cu(3)
  integer :: j, k, l, dir

  u = 1
  du = 0
  cu = 0
  do j = 1, s%dim
    dir = s%infin(j)
    call factor(dir, e(j), x(j), u(j), du(j), cu(j))
  end do

  v = product(u)
  g = 0
  h = 0
  do k = 1, s%dim
    g(k) = du(k)*product(u, mask=[(l /= k, l = 1, 3)])
    h(k,k) = cu(k)*product(u, mask=[(l /= k, l = 1, 3)])
    do l = 1, s%dim
      if (l /= k) h(k,l) = du(k)*du(l)*product(u, mask=[(j /= k .and. j /= l, j = 1, 3)])
    end do
  end do

end subroutine test_function

!***********************************************************************

pure subroutine factor(family, e, x, u, du, cu)
!! One-dimensional test function of a family, with derivatives

  integer, intent(in) :: family, e
  real(rk), intent(in) :: x
  real(rk), intent(out) :: u, du, cu

  real(rk) :: m, dm, cm   ! Monomial x**e

  m = x**e
  dm = 0
  cm = 0
  if (e >= 1) dm = e*x**(e - 1)
  if (e >= 2) cm = e*(e - 1)*x**(e - 2)

  select case (family)
  case (SHP_INF)    ! m/(1-x)
    u  = m/(1 - x)
    du = dm/(1 - x) + m/(1 - x)**2
    cu = cm/(1 - x) + 2*dm/(1 - x)**2 + 2*m/(1 - x)**3
  case (SHP_CHP)    ! (1-x) m
    u  = (1 - x)*m
    du = (1 - x)*dm - m
    cu = (1 - x)*cm - 2*dm
  case default
    u  = m
    du = dm
    cu = cm
  end select

end subroutine factor

!***********************************************************************

subroutine check_finite_differences(s, q, label)
!! Derivatives against fourth-order central differences, with a step
!! that shrinks near the singularity at \(\xi=1\), and symmetry of the
!! second derivatives

  type(shapefunc), intent(in) :: s
  type(cubature), intent(in) :: q
  character(*), intent(in) :: label

  real(rk), parameter :: w(4) = [1, -8, 8, -1]/12.0_rk
  real(rk), parameter :: off(4) = [-2, -1, 1, 2]
  integer :: g, j, m
  real(rk) :: h
  real(rk) :: x(s%dim), f(s%nnodes), d(s%dim,s%nnodes), c(s%dim,s%dim,s%nnodes)
  real(rk) :: fd(s%nnodes), dd(s%dim,s%nnodes)
  logical :: ok, sym

  ok = .true.
  sym = .true.
  do g = 1, s%npoints
    do j = 1, s%dim
      fd = 0
      dd = 0
      h = 1.0e-3_rk*min(1.0_rk, 1 - q%abscissae(j,g))
      do m = 1, 4
        x = q%abscissae(:,g)
        x(j) = x(j) + off(m)*h
        call s%eval(x, f, d, c)
        fd = fd + w(m)/h*f
        dd = dd + w(m)/h*d
      end do
      ok = ok .and. close(s%derv(j,:,g), fd, fdtol) &
              .and. close(reshape(s%curv(:,j,:,g), [s%dim*s%nnodes]), &
                          reshape(dd, [s%dim*s%nnodes]), fdtol)
      do m = 1, s%dim
        sym = sym .and. all(s%curv(j,m,:,g) == s%curv(m,j,:,g))
      end do
    end do
  end do
  call expect(ok, trim(label)//" finite differences")
  call expect(sym, trim(label)//" symmetric second derivatives")

end subroutine check_finite_differences

!***********************************************************************

subroutine check_sums(s, q, label)
!! Without chopped directions the functions sum to 1 and their derivatives
!! to 0; chopped functions vanish at infinity

  type(shapefunc), intent(in) :: s
  type(cubature), intent(in) :: q
  character(*), intent(in) :: label

  integer :: g
  real(rk) :: x(s%dim), f(s%nnodes), d(s%dim,s%nnodes), c(s%dim,s%dim,s%nnodes)
  logical :: chopped(s%dim), ok

  chopped = s%infin(1:s%dim) == SHP_CHP
  ok = .true.
  do g = 1, s%npoints
    if (any(chopped)) then
      x = merge(1.0_rk, q%abscissae(:,g), chopped)
      call s%eval(x, f, d, c)
      ok = ok .and. close(f, 0*f, tol)
    else
      ! Tolerance relative to the terms, which grow near infinity
      ok = ok .and. abs(sum(s%func(:,g)) - 1) <= tol*max(1.0_rk, sum(abs(s%func(:,g)))) &
              .and. all(abs(sum(s%derv(:,:,g), dim=2)) <= tol*max(1.0_rk, sum(abs(s%derv(:,:,g))))) &
              .and. all(abs(sum(s%curv(:,:,:,g), dim=3)) <= tol*max(1.0_rk, sum(abs(s%curv(:,:,:,g)))))
    end if
  end do
  if (any(chopped)) then
    call expect(ok, trim(label)//" chopped functions vanish at infinity")
  else
    call expect(ok, trim(label)//" partition of unity")
  end if

end subroutine check_sums

!***********************************************************************

subroutine test_numbering()
!! Orders 1 to 3 reproduce the classical node numbering

  ! Node number at each tensor-product position, first coordinate fastest
  call tensor(CUB_LIN, 1, [1, 2])
  call tensor(CUB_LIN, 2, [1, 3, 2])
  call tensor(CUB_LIN, 3, [1, 3, 4, 2])
  call tensor(CUB_LIN, 4, [1, 3, 4, 5, 2])
  call tensor(CUB_QUA, 1, [1, 2, 4, 3])
  call tensor(CUB_QUA, 2, [1, 5, 2, 8, 9, 6, 4, 7, 3])
  call tensor(CUB_QUA, 3, [1, 5, 6, 2, 12, 13, 14, 7, 11, 16, 15, 8, 4, 10, 9, 3])
  call tensor(CUB_HEX, 1, [1, 2, 4, 3, 5, 6, 8, 7])
  call tensor(CUB_HEX, 2, [ 1,  9,  2, 12, 21, 10,  4, 11,  3, &
                           17, 25, 18, 23, 27, 24, 20, 26, 19, &
                            5, 13,  6, 16, 22, 14,  8, 15,  7])

  ! Serendipity nodes are the vertex and edge nodes of the Lagrange element
  call serendipity_nodes(CUB_QUA, 2)
  call serendipity_nodes(CUB_QUA, 3)
  call serendipity_nodes(CUB_HEX, 2)
  call serendipity_nodes(CUB_HEX, 3)

  ! Lattice indices of each node
  call lattice(CUB_TRI, 1, [1,0,0, 0,1,0, 0,0,1])
  call lattice(CUB_TRI, 2, [2,0,0, 0,2,0, 0,0,2, 1,1,0, 0,1,1, 1,0,1])
  call lattice(CUB_TRI, 3, [3,0,0, 0,3,0, 0,0,3, 2,1,0, 1,2,0, 0,2,1, 0,1,2, 1,0,2, 2,0,1, 1,1,1])
  call lattice(CUB_TET, 1, [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1])
  call lattice(CUB_TET, 2, [2,0,0,0, 0,2,0,0, 0,0,2,0, 0,0,0,2, &
                            1,1,0,0, 0,1,1,0, 1,0,1,0, 1,0,0,1, 0,1,0,1, 0,0,1,1])
  call lattice(CUB_TET, 3, [3,0,0,0, 0,3,0,0, 0,0,3,0, 0,0,0,3, &
                            2,1,0,0, 1,2,0,0, 0,2,1,0, 0,1,2,0, 1,0,2,0, 2,0,1,0, &
                            2,0,0,1, 1,0,0,2, 0,2,0,1, 0,1,0,2, 0,0,2,1, 0,0,1,2, &
                            1,1,1,0, 1,1,0,1, 0,1,1,1, 1,0,1,1])
  call lattice(CUB_WED, 1, [1,0,0,0, 0,1,0,0, 0,0,1,0, 1,0,0,1, 0,1,0,1, 0,0,1,1])
  call lattice(CUB_WED, 2, [2,0,0,0, 0,2,0,0, 0,0,2,0, 2,0,0,2, 0,2,0,2, 0,0,2,2, &
                            1,1,0,0, 0,1,1,0, 1,0,1,0, 1,1,0,2, 0,1,1,2, 1,0,1,2, &
                            2,0,0,1, 0,2,0,1, 0,0,2,1, 1,1,0,1, 0,1,1,1, 1,0,1,1])

end subroutine test_numbering

!***********************************************************************

subroutine tensor(elm, p, seq)
!! Check the node number at each tensor-product position
  integer, intent(in) :: elm, p, seq(:)
  type(shapefunc) :: s
  integer :: n, a(3)
  logical :: ok
  character(20) :: label
  s = shapefunc(cubature(elm, 1), p)
  ok = size(seq) == s%nnodes
  do n = 1, size(seq)
    a = [mod(n - 1, p + 1), mod((n - 1)/(p + 1), p + 1), (n - 1)/(p + 1)**2]
    ok = ok .and. all(s%lattice(:,seq(n)) == a(1:s%dim))
  end do
  write(label,"(I0,' p=',I0)") elm, p
  call expect(ok, trim(label)//" numbering")
end subroutine tensor

!***********************************************************************

subroutine serendipity_nodes(elm, p)
!! Serendipity nodes are the leading Lagrange nodes, on vertices and edges
  integer, intent(in) :: elm, p
  type(shapefunc) :: s, l
  character(20) :: label
  s = shapefunc(cubature(elm, 1), p, family=SHP_SERENDIPITY)
  l = shapefunc(cubature(elm, 1), p)
  write(label,"(I0,' p=',I0)") elm, p
  call expect(all(s%lattice == l%lattice(:,1:s%nnodes)) &
              .and. all(count(s%lattice > 0 .and. s%lattice < p, dim=1) <= 1), &
              trim(label)//" serendipity numbering")
end subroutine serendipity_nodes

!***********************************************************************

subroutine lattice(elm, p, a)
!! Check the lattice indices of every node
  integer, intent(in) :: elm, p, a(:)
  type(shapefunc) :: s
  character(20) :: label
  s = shapefunc(cubature(elm, 1), p)
  write(label,"(I0,' p=',I0)") elm, p
  call expect(size(a) == size(s%lattice), trim(label)//" numbering size")
  if (size(a) == size(s%lattice)) &
    call expect(all(reshape(a, shape(s%lattice)) == s%lattice), trim(label)//" numbering")
end subroutine lattice

!***********************************************************************

subroutine test_closed_forms()
!! Compare with closed forms: linear mapping and chopped functions of
!! Zienkiewicz, Emson and Bettess, quadratic Lagrange, linear tetrahedron

  type(shapefunc) :: s
  real(rk) :: x, xi(3)
  real(rk), allocatable :: f(:), d(:,:), c(:,:,:)

  x = 0.3_rk

  ! Mapping functions M = [-2x, 1+x]/(1-x), nodes at -1 and 0
  s = shapefunc(cubature(CUB_LIN, 1), 1, [SHP_INF])
  call evaluate(s, [x], f, d, c)
  call expect(close(f,        [-2*x, 1 + x]/(1 - x), tol),         "INF 1 values")
  call expect(close(d(1,:),   [-2.0_rk, 2.0_rk]/(1 - x)**2, tol),  "INF 1 first derivatives")
  call expect(close(c(1,1,:), [-4.0_rk, 4.0_rk]/(1 - x)**3, tol),  "INF 1 second derivatives")
  call expect(close(s%coords(1,:), [-1.0_rk, 0.0_rk], tol),          "INF 1 nodes")

  ! Chopped linear: Lagrange on -1, 0, 1 without the node at 1
  s = shapefunc(cubature(CUB_LIN, 1), 1, [SHP_CHP])
  call evaluate(s, [x], f, d, c)
  call expect(close(f,        [x*(x - 1)/2, 1 - x**2], tol),       "CHP 1 values")
  call expect(close(d(1,:),   [x - 0.5_rk, -2*x], tol),            "CHP 1 first derivatives")
  call expect(close(c(1,1,:), [1.0_rk, -2.0_rk], tol),             "CHP 1 second derivatives")

  ! Quadratic Lagrange on -1, 1, 0 (vertices first)
  s = shapefunc(cubature(CUB_LIN, 1), 2)
  call evaluate(s, [x], f, d, c)
  call expect(close(f,        [x*(x - 1)/2, x*(x + 1)/2, 1 - x**2], tol), "LIN 2 values")
  call expect(close(d(1,:),   [x - 0.5_rk, x + 0.5_rk, -2*x], tol),      "LIN 2 first derivatives")
  call expect(close(c(1,1,:), [1.0_rk, 1.0_rk, -2.0_rk], tol),           "LIN 2 second derivatives")

  ! 8-node serendipity quadrilateral: corner 1 and midside 5
  xi(1:2) = [0.3_rk, -0.6_rk]
  s = shapefunc(cubature(CUB_QUA, 1), 2, family=SHP_SERENDIPITY)
  call evaluate(s, xi(1:2), f, d, c)
  associate (x => xi(1), y => xi(2))
    call expect(close(f([1, 5]), [(1 - x)*(1 - y)*(-x - y - 1)/4, (1 - x**2)*(1 - y)/2], tol), &
                "QUA 8 values")
    call expect(close(d(:,1), [(1 - y)*(2*x + y)/4, (1 - x)*(x + 2*y)/4], tol), &
                "QUA 8 first derivatives")
    call expect(close(reshape(c(:,:,5), [4]), [-(1 - y), x, x, 0.0_rk], tol), &
                "QUA 8 second derivatives")
  end associate

  ! Linear tetrahedron: barycentric coordinates
  xi = [0.1_rk, 0.2_rk, 0.3_rk]
  s = shapefunc(cubature(CUB_TET, 1), 1)
  call evaluate(s, xi, f, d, c)
  call expect(close(f, [xi, 1 - sum(xi)], tol), "TET 1 values")
  call expect(close(reshape(d, [12]), [1,0,0, 0,1,0, 0,0,1, -1,-1,-1]*1.0_rk, tol), &
              "TET 1 first derivatives")
  call expect(all(c == 0), "TET 1 second derivatives")

end subroutine test_closed_forms

!***********************************************************************

subroutine evaluate(s, x, f, d, c)
!! Allocate the outputs and evaluate at `x`
  type(shapefunc), intent(in) :: s
  real(rk), intent(in) :: x(:)
  real(rk), allocatable, intent(out) :: f(:), d(:,:), c(:,:,:)
  allocate(f(s%nnodes), d(s%dim,s%nnodes), c(s%dim,s%dim,s%nnodes))
  call s%eval(x, f, d, c)
end subroutine evaluate

!***********************************************************************

subroutine test_ordering()
!! A caller-defined numbering permutes every output consistently

  call check_ordering(CUB_QUA, 2, SHP_LAGRANGE, [SHP_FIN, SHP_INF])
  call check_ordering(CUB_HEX, 2, SHP_SERENDIPITY, [SHP_FIN])
  call check_ordering(CUB_TET, 2, SHP_LAGRANGE, [SHP_FIN])
  call check_ordering(CUB_WED, 3, SHP_LAGRANGE, [SHP_FIN, SHP_FIN, SHP_CHP])

end subroutine test_ordering

!***********************************************************************

subroutine check_ordering(elm, p, family, infin)
!! Compare a reversed and an identity numbering with the default one

  integer, intent(in) :: elm, p, family, infin(:)

  type(cubature) :: q
  type(shapefunc) :: s, r, id
  integer :: k
  integer, allocatable :: perm(:)
  real(rk), allocatable :: f(:), d(:,:), c(:,:,:), fr(:), dr(:,:), cr(:,:,:)
  character(20) :: label

  write(label,"(I0,' p=',I0,' fam=',I0)") elm, p, family
  q = cubature(elm, 2*p)
  s = shapefunc(q, p, infin, family)
  perm = [(s%nnodes + 1 - k, k = 1, s%nnodes)]
  r = shapefunc(q, p, infin, family, nodes=perm)
  id = shapefunc(q, p, infin, family, nodes=[(k, k = 1, s%nnodes)])

  call expect(r%is_valid() .and. r%nnodes == s%nnodes, trim(label)//" reordered is_valid")
  call expect(all(r%lattice == s%lattice(:,perm)) .and. all(r%coords == s%coords(:,perm)), &
              trim(label)//" reordered lattice and coordinates")
  call expect(all(abs(r%func - s%func(perm,:)) <= tol) &
              .and. all(abs(r%derv - s%derv(:,perm,:)) <= tol*max(1.0_rk, maxval(abs(s%derv)))) &
              .and. all(abs(r%curv - s%curv(:,:,perm,:)) <= tol*max(1.0_rk, maxval(abs(s%curv)))), &
              trim(label)//" reordered values and derivatives")
  call evaluate(s, q%abscissae(:,1), f, d, c)
  call evaluate(r, q%abscissae(:,1), fr, dr, cr)
  call expect(all(abs(fr - f(perm)) <= tol), trim(label)//" reordered eval")
  call expect(all(id%func == s%func) .and. all(id%lattice == s%lattice), trim(label)//" identity ordering")

end subroutine check_ordering

!***********************************************************************

subroutine test_integration()
!! Integrals of the shape functions with a cubature: the sum is the
!! element measure, and linear simplex functions integrate to the same value

  real(rk), parameter :: measure(6) = [2.0_rk, 0.5_rk, 4.0_rk, 1.0_rk/6, 8.0_rk, 1.0_rk]
  type(cubature) :: q
  type(shapefunc) :: s
  integer :: elm

  do elm = CUB_LIN, CUB_WED
    q = cubature(elm, 4)
    s = shapefunc(q, 2)
    call expect(close([sum(matmul(s%func, q%weights))], [measure(elm)], tol), &
                "measure of element "//achar(iachar('0') + elm))
  end do

  q = cubature(CUB_TRI, 1)
  s = shapefunc(q, 1)
  call expect(close(matmul(s%func, q%weights), [1, 1, 1]/6.0_rk, tol), "TRI 1 integrals")
  q = cubature(CUB_TET, 1)
  s = shapefunc(q, 1)
  call expect(close(matmul(s%func, q%weights), [1, 1, 1, 1]/24.0_rk, tol), "TET 1 integrals")

end subroutine test_integration

!***********************************************************************

subroutine test_interface()
!! Constructor against `set`, broadcast of `infin`, reuse, `eval`, `destroy`

  type(cubature) :: q
  type(shapefunc) :: a, b
  real(rk) :: f(9), d(2,9), c(2,2,9)

  q = cubature(CUB_QUA, 4)
  a = shapefunc(q, 2, [SHP_INF])
  call b%set(q, 2, [SHP_INF, SHP_INF])
  call expect(all(a%infin == b%infin), "infin broadcast")
  call expect(all(a%func == b%func) .and. all(a%derv == b%derv) .and. all(a%curv == b%curv), &
              "constructor equals set")

  call a%eval(q%abscissae(:,3), f, d, c)
  call expect(all(f == a%func(:,3)) .and. all(d == a%derv(:,:,3)) .and. all(c == a%curv(:,:,:,3)), &
              "eval equals stored values")

  call b%set(cubature(CUB_TET, 2), 3)
  call expect(b%is_valid() .and. b%elm == CUB_TET .and. b%nnodes == 20 .and. all(b%infin == SHP_FIN), &
              "reuse with another element")

  call b%destroy()
  call expect(.not. b%is_valid() .and. b%elm == 0 .and. b%nnodes == 0, "destroy")

  a = shapefunc(cubature(CUB_HEX, 1), 0)
  call expect(a%nnodes == 1 .and. all(a%coords == 0) .and. all(a%func == 1) .and. all(a%derv == 0), &
              "order 0")

end subroutine test_interface

!***********************************************************************

subroutine test_output()
!! `summary`, `show` and `numbering` write to a unit, also when unset

  type(shapefunc) :: s
  integer :: u, n, ios
  character(80) :: line

  open(newunit=u, status="scratch", action="readwrite")
  call s%summary(u)
  s = shapefunc(cubature(CUB_WED, 2), 1, [SHP_FIN, SHP_FIN, SHP_CHP])
  call s%summary(u)
  call s%show(u)
  call s%numbering(u)
  rewind(u)
  read(u,"(A)") line
  call expect(line == "shapefunc: not set", "summary of unset")
  n = 1
  do
    read(u,"(A)", iostat=ios) line
    if (ios /= 0) exit
    n = n + 1
  end do
  close(u)
  ! Unset summary, summary, show, numbering
  call expect(n == 1 + 7 + (7 + 1 + s%npoints*(s%nnodes + 1)) + (7 + 1 + s%nnodes), &
              "output line count")

end subroutine test_output

!***********************************************************************

end program shapefuncs_test
