module shapefuncs
!! Lagrange shape functions and their first and second derivatives on the
!! reference elements of the finite element method, including mapped
!! infinite elements.
!!
!! A [[shapefunc]] evaluates the \(n\) shape functions \(N_i(\xi)\) of an
!! element at the points of a `cubature`, together with
!! \(\partial N_i/\partial\xi_j\) and \(\partial^2 N_i/\partial\xi_j\partial\xi_k\).
!! The reference elements and coordinates are those of the `cubatures`
!! package.
!!
!! ## Elements and orders
!!
!! | Constant  | Element       | Reference domain            | Nodes, order \(p\)          | Orders |
!! |-----------|---------------|-----------------------------|-----------------------------|--------|
!! | `CUB_LIN` | line          | \([-1,1]\)                  | \(p+1\)                     | \(\ge 0\) |
!! | `CUB_QUA` | quadrilateral | \([-1,1]^2\)                | \((p+1)^2\)                 | \(\ge 0\) |
!! | `CUB_HEX` | hexahedron    | \([-1,1]^3\)                | \((p+1)^3\)                 | \(\ge 0\) |
!! | `CUB_TRI` | triangle      | \(x,y \ge 0,\ x+y \le 1\)     | \((p+1)(p+2)/2\)            | \(\ge 1\) |
!! | `CUB_TET` | tetrahedron   | \(x,y,z \ge 0,\ x+y+z \le 1\) | \((p+1)(p+2)(p+3)/6\)       | \(\ge 1\) |
!! | `CUB_WED` | wedge / prism | triangle \(\times\ [-1,1]\) | \((p+1)^2(p+2)/2\)          | \(\ge 1\) |
!!
!! Lines, quadrilaterals and hexahedra use tensor products of Lagrange
!! polynomials on equispaced nodes; order 0 is the constant function with
!! one node at the centre. Triangles and tetrahedra use the complete
!! polynomials of total degree \(p\) in the barycentric coordinates
!! (Silvester's polynomials [5]); wedges use the triangle functions times
!! Lagrange polynomials in the axial coordinate. There is no upper limit
!! on the order, but equispaced interpolation becomes ill-conditioned at
!! high orders.
!!
!! ## Factor coordinates
!!
!! Every shape function is a product of one-dimensional functions of
!! *factor coordinates* \(s_k\), which are affine in \(\xi\):
!!
!! - line, quadrilateral, hexahedron: \(s = \xi\);
!! - triangle: \(s = (\xi_1, \xi_2, 1-\xi_1-\xi_2)\), the barycentric coordinates;
!! - tetrahedron: \(s = (\xi_1, \xi_2, \xi_3, 1-\xi_1-\xi_2-\xi_3)\);
!! - wedge: \(s = (\xi_1, \xi_2, 1-\xi_1-\xi_2, \xi_3)\).
!!
!! Node \(i\) has one integer *lattice index* \(a_{ki}\) per factor
!! coordinate, and \(N_i(\xi) = \prod_k \phi_{a_{ki}}(s_k)\). Derivatives
!! follow from the product rule and the constant \(\partial s/\partial\xi\).
!!
!! ## Node numbering
!!
!! Nodes are numbered hierarchically: vertices, then edge nodes, then face
!! nodes, then interior nodes.
!!
!! - **Vertices** of a quadrilateral run counter-clockwise from
!!   \((-1,-1)\). A hexahedron has the bottom face \(\xi_3=-1\) counter-clockwise from
!!   \((-1,-1,-1)\), then the top face in the same order. Triangle and
!!   tetrahedron vertex \(k\) is where barycentric coordinate \(k\) equals
!!   1, so vertex 1 is \(\xi = (1,0,\dots)\) and the last vertex is the
!!   origin. Wedge vertices are the triangle vertices at \(\xi_3=-1\), then
!!   at \(\xi_3=1\).
!! - **Edges** are listed below. Edge nodes run from the first vertex of
!!   the edge to the second.
!!
!!   | Element       | Edges |
!!   |---------------|-------|
!!   | quadrilateral | 1-2, 2-3, 3-4, 4-1 |
!!   | hexahedron    | 1-2, 2-3, 3-4, 4-1, 5-6, 6-7, 7-8, 8-5, 1-5, 2-6, 3-7, 4-8 |
!!   | triangle      | 1-2, 2-3, 3-1 |
!!   | tetrahedron   | 1-2, 2-3, 3-1, 1-4, 2-4, 3-4 |
!!   | wedge         | 1-2, 2-3, 3-1, 4-5, 5-6, 6-4, 1-4, 2-5, 3-6 |
!!
!! - **Faces** of a hexahedron are \(\xi_3=-1\), \(\xi_3=1\), \(\xi_1=-1\),
!!   \(\xi_1=1\), \(\xi_2=-1\), \(\xi_2=1\). Tetrahedron faces are 1-2-3,
!!   1-2-4, 2-3-4, 1-3-4. Wedge faces are the bottom and top triangles,
!!   then the quadrilaterals on edges 1-2, 2-3, 3-1.
!! - **Interiors** of quadrilateral faces and elements, and of triangles
!!   and tetrahedra, are numbered recursively as an element of the same
!!   type of order \(p-2\) (quadrilateral, hexahedron), \(p-3\) (triangle)
!!   or \(p-4\) (tetrahedron). A quadrilateral face in \((\xi_j, \xi_k)\),
!!   \(j<k\), is numbered as a quadrilateral in those coordinates; a
!!   triangular face as a triangle with the vertices in the order listed.
!!   Wedge interior nodes are triangle interiors, layer by layer from
!!   \(\xi_3=-1\).
!!
!! These rules reproduce the classical numbering of the 4-, 9- and 16-node
!! quadrilateral, 8- and 27-node hexahedron, 3-, 6- and 10-node triangle,
!! 4-, 10- and 20-node tetrahedron and 6- and 18-node wedge.
!!
!! ## Infinite elements
!!
!! Along a direction flagged [[SHP_INF]] or [[SHP_CHP]], the element
!! extends to infinity at \(\xi=1\), and the \(p+1\) nodes sit at
!! \(t_a = -1 + 2a/(p+1)\), \(a = 0,\dots,p\).
!!
!! - [[SHP_INF]] gives the *mapping* functions of Zienkiewicz, Emson and
!!   Bettess [2], generalised to order \(p\) [3]:
!!   \[ M_a(\xi) = \ell_a(\xi)\,\frac{1-t_a}{1-\xi}, \]
!!   where \(\ell_a\) are the Lagrange polynomials on the nodes \(t_a\).
!!   They satisfy \(M_a(t_b)=\delta_{ab}\) and \(\sum_a M_a = 1\), and
!!   map \(\xi \to 1\) to infinity. Use them to interpolate coordinates.
!! - [[SHP_CHP]] gives the *chopped* functions: the Lagrange polynomials
!!   of order \(p+1\) on the nodes \(t_0,\dots,t_{p+1}=1\), without the
!!   function of the node at infinity. They interpolate a field that
!!   vanishes at infinity, and so do not sum to 1. Because the mapping
!!   makes \(1-\xi\) proportional to \(1/r\), they represent the decay
!!   \(\sum_{k=1}^{p+1} \alpha_k/r^k\) [2, 4]. Use them to interpolate
!!   the solution.
!!
!! ## References
!!
!! 1. Bettess, P. (1977). Infinite elements. *International Journal for
!!    Numerical Methods in Engineering*, 11(1), 53–64. doi:10.1002/nme.1620110107
!! 2. Zienkiewicz, O. C., Emson, C., Bettess, P. (1983). A novel boundary
!!    infinite element. *International Journal for Numerical Methods in
!!    Engineering*, 19(3), 393–404. doi:10.1002/nme.1620190307
!! 3. Marques, J. M. M. C., Owen, D. R. J. (1984). Infinite elements in
!!    quasi-static materially nonlinear problems. *Computers & Structures*,
!!    18(4), 739–751. doi:10.1016/0045-7949(84)90019-1
!! 4. Bettess, P. (1992). *Infinite Elements*. Penshaw Press.
!! 5. Silvester, P. (1969). High-order polynomial triangular finite
!!    elements for potential problems. *International Journal of
!!    Engineering Science*, 7(8), 849–861. doi:10.1016/0020-7225(69)90065-2
!! 6. Zienkiewicz, O. C., Taylor, R. L., Zhu, J. Z. (2005). *The Finite
!!    Element Method: Its Basis and Fundamentals*, 6th ed.
!!    Butterworth-Heinemann.

  use iso_fortran_env, only: output_unit
  use cubatures, only: cubature, rk, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED

  implicit none

  private

  public :: shapefunc
  public :: SHP_FIN, SHP_INF, SHP_CHP

  integer, parameter :: SHP_FIN = 1 !! Finite direction: Lagrange polynomials on \([-1,1]\)
  integer, parameter :: SHP_INF = 2 !! Infinite direction, mapping functions for coordinates
  integer, parameter :: SHP_CHP = 3 !! Infinite direction, chopped functions for the solution

  integer, parameter :: SIMPLEX = 4 !! Factor family: Silvester polynomials on \([0,1]\)

  character(3), parameter :: names(6) = ["LIN", "TRI", "QUA", "TET", "HEX", "WED"]
    !! Element names, indexed by element type
  character(3), parameter :: infnames(3) = ["FIN", "INF", "CHP"]
    !! Infinitude names, indexed by `SHP_*`

  type :: shapefunc
    !! Shape functions and derivatives of one element type and order at
    !! the points of a cubature.
    !!
    !! Create with the constructor, `s = shapefunc(q, 2)`, or in place with
    !! `call s%set(q, 2)`, where `q` is a `cubature`. Infinite directions
    !! are flagged with the optional `infin` argument.

    integer :: elm = 0                          !! Element type, `CUB_*`
    integer :: dim = 0                          !! Spatial dimension
    integer :: order = 0                        !! Polynomial order
    integer :: nnodes = 0                       !! Number of nodes
    integer :: npoints = 0                      !! Number of cubature points
    integer :: infin(3) = SHP_FIN               !! Infinitude per direction, `SHP_*`
    integer, allocatable :: lattice(:,:)        !! Lattice indices, shape `[nfactors, nnodes]`, see module documentation
    real(rk), allocatable :: coords(:,:)        !! Nodal coordinates, shape `[dim, nnodes]`
    real(rk), allocatable :: func(:,:)          !! Values \(N_i\), shape `[nnodes, npoints]`
    real(rk), allocatable :: derv(:,:,:)        !! First derivatives \(\partial N_i/\partial\xi_j\), shape `[dim, nnodes, npoints]`
    real(rk), allocatable :: curv(:,:,:,:)      !! Second derivatives \(\partial^2 N_i/\partial\xi_j\partial\xi_k\), shape `[dim, dim, nnodes, npoints]`

    integer, private :: nfactors = 0            ! Number of factor coordinates
    integer, private :: family(4) = 0           ! 1D function family per factor, `SHP_*` or `SIMPLEX`
    real(rk), private :: s0(4) = 0              ! Factor coordinates at \(\xi = 0\)
    real(rk), private :: dsdxi(4,3) = 0         ! Constant \(\partial s_k/\partial\xi_j\)

  contains

    procedure :: set
    procedure :: eval
    procedure :: is_valid
    procedure :: summary
    procedure :: show
    procedure :: numbering
    procedure :: destroy

  end type shapefunc

  interface shapefunc
    !! Construct shape functions, e.g. `shapefunc(q, 2)` or
    !! `shapefunc(q, 1, [SHP_FIN, SHP_INF])`
    module procedure new_shapefunc
  end interface shapefunc

contains

!***********************************************************************

pure function new_shapefunc(q, order, infin) result(self)
!! Construct shape functions of order `order` at the points of `q`

  type(cubature), intent(in) :: q            !! Cubature on the reference element
  integer, intent(in) :: order               !! Polynomial order
  integer, intent(in), optional :: infin(:)  !! Infinitude, `SHP_*`, size 1 or `q%dim`
  type(shapefunc) :: self

  call self%set(q, order, infin)

end function new_shapefunc

!***********************************************************************

pure subroutine set(self, q, order, infin)
!! Build shape functions of order `order` at the points of `q`.
!!
!! `infin` has size 1 (same in every direction) or `q%dim`, and defaults to
!! [[SHP_FIN]]. [[SHP_INF]] and [[SHP_CHP]] are allowed along the
!! directions of lines, quadrilaterals and hexahedra, and along the axial
!! direction 3 of wedges, with order \(\ge 1\). Stops on invalid input.

  class(shapefunc), intent(inout) :: self
  type(cubature), intent(in) :: q            !! Cubature on the reference element
  integer, intent(in) :: order               !! Polynomial order
  integer, intent(in), optional :: infin(:)  !! Infinitude, `SHP_*`, size 1 or `q%dim`

  integer :: g, k, p

  if (.not. q%is_valid()) error stop "shapefunc%set: cubature is not set"
  if (order < 0) error stop "shapefunc%set: order must be non-negative"

  call self%destroy()
  self%elm   = q%elm
  self%dim   = q%dim
  self%order = order
  p = order

  if (present(infin)) then
    if (size(infin) == 1) then
      self%infin(1:self%dim) = infin(1)
    else if (size(infin) == self%dim) then
      self%infin(1:self%dim) = infin
    else
      error stop "shapefunc%set: size(infin) must be 1 or the dimension"
    end if
    if (any(self%infin < SHP_FIN .or. self%infin > SHP_CHP)) &
      error stop "shapefunc%set: infin must be SHP_FIN, SHP_INF or SHP_CHP"
    if (any(self%infin /= SHP_FIN) .and. p < 1) &
      error stop "shapefunc%set: infinite directions need order >= 1"
  end if

  ! Factor coordinates s = s0 + dsdxi xi, and the 1D family of each
  select case (self%elm)
  case (CUB_LIN, CUB_QUA, CUB_HEX)
    self%nfactors = self%dim
    self%family(1:self%dim) = self%infin(1:self%dim)
    do k = 1, self%dim
      self%dsdxi(k,k) = 1
    end do
  case (CUB_TRI, CUB_TET)
    if (any(self%infin /= SHP_FIN)) error stop "shapefunc%set: simplices must be finite"
    self%nfactors = self%dim + 1
    self%family = SIMPLEX
    do k = 1, self%dim
      self%dsdxi(k,k) = 1
    end do
    self%dsdxi(self%dim+1,1:self%dim) = -1
    self%s0(self%dim+1) = 1
  case (CUB_WED)
    if (any(self%infin(1:2) /= SHP_FIN)) &
      error stop "shapefunc%set: wedges may be infinite only along direction 3"
    self%nfactors = 4
    self%family(1:3) = SIMPLEX
    self%family(4)   = self%infin(3)
    self%dsdxi(1,1) = 1
    self%dsdxi(2,2) = 1
    self%dsdxi(3,1:2) = -1
    self%dsdxi(4,3) = 1
    self%s0(3) = 1
  end select

  if (any(self%family == SIMPLEX) .and. p < 1) &
    error stop "shapefunc%set: simplex and wedge orders must be >= 1"

  ! Node numbering and coordinates
  select case (self%elm)
  case (CUB_LIN)
    self%lattice = line_lattice(p)
  case (CUB_QUA)
    self%lattice = quad_lattice(p)
  case (CUB_HEX)
    self%lattice = hex_lattice(p)
  case (CUB_TRI)
    self%lattice = tri_lattice(p)
  case (CUB_TET)
    self%lattice = tet_lattice(p)
  case (CUB_WED)
    self%lattice = wedge_lattice(p)
  end select
  self%nnodes = size(self%lattice, 2)

  allocate(self%coords(self%dim, self%nnodes))
  select case (self%elm)
  case (CUB_LIN, CUB_QUA, CUB_HEX)
    do k = 1, self%dim
      self%coords(k,:) = node_coordinate(self%family(k), p, self%lattice(k,:))
    end do
  case (CUB_TRI, CUB_TET)
    self%coords = real(self%lattice(1:self%dim,:), rk)/p
  case (CUB_WED)
    self%coords(1:2,:) = real(self%lattice(1:2,:), rk)/p
    self%coords(3,:)   = node_coordinate(self%family(4), p, self%lattice(4,:))
  end select

  ! Values and derivatives at the cubature points
  self%npoints = q%npoints
  allocate(self%func(self%nnodes, self%npoints))
  allocate(self%derv(self%dim, self%nnodes, self%npoints))
  allocate(self%curv(self%dim, self%dim, self%nnodes, self%npoints))
  do g = 1, self%npoints
    call self%eval(q%abscissae(:,g), self%func(:,g), self%derv(:,:,g), self%curv(:,:,:,g))
  end do

end subroutine set

!***********************************************************************

pure subroutine eval(self, xi, func, derv, curv)
!! Evaluate the shape functions and their derivatives at any point `xi`
!! of the reference element. The element must be set. Infinite
!! directions are singular at \(\xi = 1\).

  class(shapefunc), intent(in) :: self
  real(rk), intent(in) :: xi(:)                                       !! Point, size `dim`
  real(rk), intent(out) :: func(self%nnodes)                          !! Values \(N_i\)
  real(rk), intent(out) :: derv(self%dim, self%nnodes)                !! \(\partial N_i/\partial\xi_j\)
  real(rk), intent(out) :: curv(self%dim, self%dim, self%nnodes)      !! \(\partial^2 N_i/\partial\xi_j\partial\xi_k\)

  integer :: i, k, l, n
  real(rk) :: s(4)                                                    ! Factor coordinates
  real(rk), dimension(0:self%order, 4) :: f, d, c                     ! 1D values and derivatives per factor
  real(rk) :: fi(4), di(4), ci(4)                                     ! Factors of node i
  real(rk) :: grad(4), hess(4,4)                                      ! Derivatives with respect to s
  real(rk) :: jac(4, self%dim)                                        ! ds/dxi

  n = self%nfactors
  jac = self%dsdxi(1:4, 1:self%dim)
  s = self%s0 + matmul(jac, xi(1:self%dim))

  do k = 1, n
    call basis(self%family(k), self%order, s(k), f(:,k), d(:,k), c(:,k))
  end do

  do i = 1, self%nnodes
    do k = 1, n
      fi(k) = f(self%lattice(k,i), k)
      di(k) = d(self%lattice(k,i), k)
      ci(k) = c(self%lattice(k,i), k)
    end do

    ! N = prod f_k; product rule for the s-derivatives
    func(i) = product(fi(1:n))
    do k = 1, n
      grad(k) = di(k)*product(fi(1:n), mask=skip(n, k, k))
      hess(k,k) = ci(k)*product(fi(1:n), mask=skip(n, k, k))
      do l = k + 1, n
        hess(k,l) = di(k)*di(l)*product(fi(1:n), mask=skip(n, k, l))
        hess(l,k) = hess(k,l)
      end do
    end do

    ! Chain rule with constant ds/dxi, symmetrised exactly
    derv(:,i)   = matmul(grad(1:n), jac(1:n,:))
    curv(:,:,i) = matmul(transpose(jac(1:n,:)), matmul(hess(1:n,1:n), jac(1:n,:)))
    curv(:,:,i) = (curv(:,:,i) + transpose(curv(:,:,i)))/2
  end do

end subroutine eval

!***********************************************************************

pure function skip(n, k, l) result(mask)
!! Mask of `1:n` without `k` and `l`

  integer, intent(in) :: n, k, l
  logical :: mask(n)

  integer :: j

  mask = [(j /= k .and. j /= l, j = 1, n)]

end function skip

!***********************************************************************

pure logical function is_valid(self)
!! True if the shape functions are set and their arrays are consistent

  class(shapefunc), intent(in) :: self

  is_valid = allocated(self%lattice) .and. allocated(self%coords) .and. allocated(self%func) &
       .and. allocated(self%derv) .and. allocated(self%curv)
  if (.not. is_valid) return
  is_valid = self%nnodes > 0 .and. self%npoints > 0 &
       .and. all(shape(self%lattice) == [self%nfactors, self%nnodes]) &
       .and. all(shape(self%coords)  == [self%dim, self%nnodes]) &
       .and. all(shape(self%func)    == [self%nnodes, self%npoints]) &
       .and. all(shape(self%derv)    == [self%dim, self%nnodes, self%npoints]) &
       .and. all(shape(self%curv)    == [self%dim, self%dim, self%nnodes, self%npoints])

end function is_valid

!***********************************************************************

subroutine summary(self, unit)
!! Write element type, dimension, order, infinitude and numbers of nodes
!! and points

  class(shapefunc), intent(in) :: self
  integer, intent(in), optional :: unit   !! Output unit, default `output_unit`

  integer :: u

  u = output_unit
  if (present(unit)) u = unit

  if (self%elm == 0) then
    write(u,"(A)") "shapefunc: not set"
    return
  end if

  write(u,"(A,A)")              "Element:    ", names(self%elm)
  write(u,"(A,I0)")             "Dimension:  ", self%dim
  write(u,"(A,I0)")             "Order:      ", self%order
  write(u,"(A,*(A,:,', '))")    "Infinitude: ", infnames(self%infin(1:self%dim))
  write(u,"(A,I0)")             "Nodes:      ", self%nnodes
  write(u,"(A,I0)")             "Points:     ", self%npoints

end subroutine summary

!***********************************************************************

subroutine show(self, unit)
!! Write the summary, then for every point and node the value and first
!! derivatives, and per point the sums over the nodes

  class(shapefunc), intent(in) :: self
  integer, intent(in), optional :: unit   !! Output unit, default `output_unit`

  integer :: u, g, i

  u = output_unit
  if (present(unit)) u = unit

  call self%summary(u)
  if (.not. self%is_valid()) return

  write(u,"(A)") "Point, node, value, first derivatives"
  do g = 1, self%npoints
    do i = 1, self%nnodes
      write(u,"(2I5,*(1X,ES23.15E3))") g, i, self%func(i,g), self%derv(:,i,g)
    end do
    write(u,"(A10,*(1X,ES23.15E3))") "Sum", sum(self%func(:,g)), sum(self%derv(:,:,g), dim=2)
  end do

end subroutine show

!***********************************************************************

subroutine numbering(self, unit)
!! Write the node numbers and nodal coordinates

  class(shapefunc), intent(in) :: self
  integer, intent(in), optional :: unit   !! Output unit, default `output_unit`

  integer :: u, i

  u = output_unit
  if (present(unit)) u = unit

  call self%summary(u)
  if (.not. self%is_valid()) return

  write(u,"(A)") "Node, coordinates"
  do i = 1, self%nnodes
    write(u,"(I5,*(1X,F10.6))") i, self%coords(:,i)
  end do

end subroutine numbering

!***********************************************************************

pure subroutine destroy(self)
!! Deallocate and reset to the unset state

  class(shapefunc), intent(inout) :: self

  self%elm      = 0
  self%dim      = 0
  self%order    = 0
  self%nnodes   = 0
  self%npoints  = 0
  self%infin    = SHP_FIN
  self%nfactors = 0
  self%family   = 0
  self%s0       = 0
  self%dsdxi    = 0
  if (allocated(self%lattice)) deallocate(self%lattice)
  if (allocated(self%coords))  deallocate(self%coords)
  if (allocated(self%func))    deallocate(self%func)
  if (allocated(self%derv))    deallocate(self%derv)
  if (allocated(self%curv))    deallocate(self%curv)

end subroutine destroy

!***********************************************************************
! One-dimensional functions
!***********************************************************************

pure subroutine basis(family, p, s, f, d, c)
!! The \(p+1\) one-dimensional functions of a family, with first and
!! second derivatives, at factor coordinate `s`

  integer, intent(in) :: family      !! `SHP_FIN`, `SHP_INF`, `SHP_CHP` or `SIMPLEX`
  integer, intent(in) :: p           !! Order
  real(rk), intent(in) :: s          !! Factor coordinate
  real(rk), intent(out) :: f(0:p)    !! Values
  real(rk), intent(out) :: d(0:p)    !! First derivatives
  real(rk), intent(out) :: c(0:p)    !! Second derivatives

  integer :: a, j
  real(rk) :: t(0:p+1)               ! Nodes
  real(rk) :: h                      ! 1/(1-s)

  select case (family)

  ! Lagrange polynomials on equispaced nodes in [-1, 1]
  case (SHP_FIN)
    t(0:p) = node_coordinate(SHP_FIN, p, [(j, j = 0, p)])
    do a = 0, p
      call lagrange(t(0:p), a, s, f(a), d(a), c(a))
    end do

  ! Mapping functions: Lagrange polynomials times (1-t_a)/(1-s)
  case (SHP_INF)
    t(0:p) = node_coordinate(SHP_INF, p, [(j, j = 0, p)])
    h = 1/(1 - s)
    do a = 0, p
      call lagrange(t(0:p), a, s, f(a), d(a), c(a))
      call multiply(f(a), d(a), c(a), (1 - t(a))*h, (1 - t(a))*h**2, 2*(1 - t(a))*h**3)
    end do

  ! Chopped: order p+1 Lagrange polynomials without the node at s = 1
  case (SHP_CHP)
    t = node_coordinate(SHP_INF, p, [(j, j = 0, p + 1)])
    do a = 0, p
      call lagrange(t, a, s, f(a), d(a), c(a))
    end do

  ! Silvester polynomials on [0, 1]: prod_{k<a} (p s - k)/(k + 1)
  case (SIMPLEX)
    do a = 0, p
      f(a) = 1
      d(a) = 0
      c(a) = 0
      do j = 0, a - 1
        call multiply(f(a), d(a), c(a), (p*s - j)/(j + 1), real(p, rk)/(j + 1), 0.0_rk)
      end do
    end do

  end select

end subroutine basis

!***********************************************************************

pure subroutine lagrange(t, a, s, f, d, c)
!! Lagrange polynomial of node `a` on the nodes `t`, with first and second
!! derivatives, at `s`

  real(rk), intent(in) :: t(0:)      !! Nodes
  integer, intent(in) :: a           !! Node index
  real(rk), intent(in) :: s          !! Point
  real(rk), intent(out) :: f, d, c   !! Value, first and second derivative

  integer :: j

  f = 1
  d = 0
  c = 0
  do j = 0, ubound(t, 1)
    if (j == a) cycle
    call multiply(f, d, c, (s - t(j))/(t(a) - t(j)), 1/(t(a) - t(j)), 0.0_rk)
  end do

end subroutine lagrange

!***********************************************************************

pure subroutine multiply(f, d, c, g, dg, cg)
!! Multiply a function by a factor, updating the value and derivatives by
!! the product rule

  real(rk), intent(inout) :: f, d, c   !! Value, first and second derivative
  real(rk), intent(in) :: g, dg, cg    !! Factor, first and second derivative

  c = c*g + 2*d*dg + f*cg
  d = d*g + f*dg
  f = f*g

end subroutine multiply

!***********************************************************************

pure function node_coordinate(family, p, a) result(x)
!! Coordinates of lattice indices `a` along a tensor direction: equispaced
!! on \([-1,1]\) for finite directions (0 for order 0), and
!! \(-1 + 2a/(p+1)\) for infinite directions

  integer, intent(in) :: family      !! `SHP_*`
  integer, intent(in) :: p           !! Order
  integer, intent(in) :: a(:)        !! Lattice indices
  real(rk) :: x(size(a))

  if (family == SHP_FIN) then
    if (p == 0) then
      x = 0
    else
      x = -1 + 2*real(a, rk)/p
    end if
  else
    x = -1 + 2*real(a, rk)/(p + 1)
  end if

end function node_coordinate

!***********************************************************************
! Node numbering
!***********************************************************************

pure function line_lattice(p) result(a)
!! Lattice of a line of order `p`: vertices, then interior nodes

  integer, intent(in) :: p
  integer, allocatable :: a(:,:)

  integer :: j

  if (p == 0) then
    a = reshape([0], [1, 1])
  else
    a = reshape([0, p, (j, j = 1, p - 1)], [1, p + 1])
  end if

end function line_lattice

!***********************************************************************

pure recursive function quad_lattice(p) result(a)
!! Lattice of a quadrilateral of order `p`: vertices, edges, interior
!! (recursively of order `p-2`)

  integer, intent(in) :: p
  integer, allocatable :: a(:,:)

  integer :: v(2,4), e, n

  if (p == 0) then
    a = reshape([0, 0], [2, 1])
    return
  end if

  v = reshape([0,0, p,0, p,p, 0,p], [2, 4])
  allocate(a(2, (p + 1)**2))
  a(:,1:4) = v
  n = 4
  do e = 1, 4
    call add_edge(a, n, v(:,e), v(:,mod(e, 4) + 1), p)
  end do
  if (p >= 2) a(:,n+1:) = quad_lattice(p - 2) + 1

end function quad_lattice

!***********************************************************************

pure recursive function hex_lattice(p) result(a)
!! Lattice of a hexahedron of order `p`: vertices, edges, faces, interior
!! (recursively of order `p-2`)

  integer, intent(in) :: p
  integer, allocatable :: a(:,:)

  integer, parameter :: edges(2,12) = reshape([1,2, 2,3, 3,4, 4,1, 5,6, 6,7, 7,8, 8,5, &
                                               1,5, 2,6, 3,7, 4,8], [2, 12])
  integer, parameter :: fixed(6) = [3, 3, 1, 1, 2, 2]             ! Direction normal to face
  integer, parameter :: plane(2,6) = reshape([1,2, 1,2, 2,3, 2,3, 1,3, 1,3], [2, 6])
  integer :: v(3,8), e, n, m
  integer, allocatable :: b(:,:)

  if (p == 0) then
    a = reshape([0, 0, 0], [3, 1])
    return
  end if

  v = reshape([0,0,0, p,0,0, p,p,0, 0,p,0, 0,0,p, p,0,p, p,p,p, 0,p,p], [3, 8])
  allocate(a(3, (p + 1)**3))
  a(:,1:8) = v
  n = 8
  do e = 1, 12
    call add_edge(a, n, v(:,edges(1,e)), v(:,edges(2,e)), p)
  end do
  if (p >= 2) then
    b = quad_lattice(p - 2) + 1
    m = size(b, 2)
    do e = 1, 6
      a(fixed(e), n+1:n+m)    = merge(0, p, mod(e, 2) == 1)
      a(plane(1,e), n+1:n+m)  = b(1,:)
      a(plane(2,e), n+1:n+m)  = b(2,:)
      n = n + m
    end do
    a(:,n+1:) = hex_lattice(p - 2) + 1
  end if

end function hex_lattice

!***********************************************************************

pure recursive function tri_lattice(p) result(a)
!! Barycentric lattice of a triangle of order `p`: vertices, edges,
!! interior (recursively of order `p-3`)

  integer, intent(in) :: p
  integer, allocatable :: a(:,:)

  integer :: v(3,3), e, n

  if (p == 0) then
    a = reshape([0, 0, 0], [3, 1])
    return
  end if

  v = p*reshape([1,0,0, 0,1,0, 0,0,1], [3, 3])
  allocate(a(3, (p + 1)*(p + 2)/2))
  a(:,1:3) = v
  n = 3
  do e = 1, 3
    call add_edge(a, n, v(:,e), v(:,mod(e, 3) + 1), p)
  end do
  if (p >= 3) a(:,n+1:) = tri_lattice(p - 3) + 1

end function tri_lattice

!***********************************************************************

pure recursive function tet_lattice(p) result(a)
!! Barycentric lattice of a tetrahedron of order `p`: vertices, edges,
!! faces, interior (recursively of order `p-4`)

  integer, intent(in) :: p
  integer, allocatable :: a(:,:)

  integer, parameter :: edges(2,6) = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4], [2, 6])
  integer, parameter :: faces(3,4) = reshape([1,2,3, 1,2,4, 2,3,4, 1,3,4], [3, 4])
  integer :: v(4,4), e, n, m
  integer, allocatable :: b(:,:)

  if (p == 0) then
    a = reshape([0, 0, 0, 0], [4, 1])
    return
  end if

  v = p*reshape([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1], [4, 4])
  allocate(a(4, (p + 1)*(p + 2)*(p + 3)/6))
  a(:,1:4) = v
  n = 4
  do e = 1, 6
    call add_edge(a, n, v(:,edges(1,e)), v(:,edges(2,e)), p)
  end do
  if (p >= 3) then
    b = tri_lattice(p - 3) + 1
    m = size(b, 2)
    do e = 1, 4
      a(:,n+1:n+m) = 0
      a(faces(:,e), n+1:n+m) = b
      n = n + m
    end do
  end if
  if (p >= 4) a(:,n+1:) = tet_lattice(p - 4) + 1

end function tet_lattice

!***********************************************************************

pure function wedge_lattice(p) result(a)
!! Lattice of a wedge of order `p`: three barycentric indices and one
!! axial index. Vertices, edges, triangular faces, quadrilateral faces,
!! interior.

  integer, intent(in) :: p
  integer, allocatable :: a(:,:)

  integer :: v(4,6), e, n, m, k, z
  integer, allocatable :: b(:,:)

  v = 0
  do k = 1, 3
    v(k,k)   = p
    v(k,k+3) = p
  end do
  v(4,4:6) = p

  allocate(a(4, (p + 1)**2*(p + 2)/2))
  a(:,1:6) = v
  n = 6
  do e = 1, 3                                                      ! Bottom edges
    call add_edge(a, n, v(:,e), v(:,mod(e, 3) + 1), p)
  end do
  do e = 1, 3                                                      ! Top edges
    call add_edge(a, n, v(:,e+3), v(:,mod(e, 3) + 4), p)
  end do
  do e = 1, 3                                                      ! Vertical edges
    call add_edge(a, n, v(:,e), v(:,e+3), p)
  end do

  if (p >= 3) then                                                 ! Bottom and top faces
    b = tri_lattice(p - 3) + 1
    m = size(b, 2)
    do z = 0, p, p
      a(1:3,n+1:n+m) = b
      a(4,n+1:n+m) = z
      n = n + m
    end do
  end if

  if (p >= 2) then                                                 ! Quadrilateral faces
    b = quad_lattice(p - 2) + 1
    m = size(b, 2)
    do e = 1, 3
      k = mod(e, 3) + 1
      a(:,n+1:n+m) = 0
      a(e,n+1:n+m) = p - b(1,:)
      a(k,n+1:n+m) = b(1,:)
      a(4,n+1:n+m) = b(2,:)
      n = n + m
    end do
  end if

  if (p >= 3) then                                                 ! Interior, layer by layer
    b = tri_lattice(p - 3) + 1
    m = size(b, 2)
    do z = 1, p - 1
      a(1:3,n+1:n+m) = b
      a(4,n+1:n+m) = z
      n = n + m
    end do
  end if

end function wedge_lattice

!***********************************************************************

pure subroutine add_edge(a, n, va, vb, p)
!! Append the `p-1` interior lattice points of the edge from `va` to `vb`

  integer, intent(inout) :: a(:,:)   !! Lattice
  integer, intent(inout) :: n        !! Number of points so far
  integer, intent(in) :: va(:), vb(:)!! Edge vertices
  integer, intent(in) :: p           !! Order

  integer :: t

  do t = 1, p - 1
    n = n + 1
    a(:,n) = va + t*(vb - va)/p
  end do

end subroutine add_edge

!***********************************************************************

end module shapefuncs
