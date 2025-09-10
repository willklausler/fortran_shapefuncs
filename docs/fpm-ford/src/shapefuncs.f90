module shapefuncs
!! Procedures for evaluating shape and map functions and derivatives for
!! triangles, quadrilateral, tetrahedron, and hexahedron elements
!! Use cubature module with derived type

  use iso_fortran_env, only: rk => real64, stdout => output_unit
  use cubatures, only: cubature, counter

  implicit none

  private

  integer, parameter :: shaporders(*) = [4, 3, 4, 3, 4, 2]

  public :: shaporders

  abstract interface
    !! Abstract interface to allow procedure pointers
    pure module subroutine lagrpoly(order,x, func,derv,curv)
      integer, intent(in) :: order
      real(rk), intent(in) :: x
      real(rk), intent(out) :: func(order+1), derv(order+1), curv(order+1)
    end subroutine lagrpoly
  end interface

  type :: Shapefunc
    !! Includes Lagrange shape functions, linear - quartic, over the domains
    !!          [-1, 1] and [0, 1]
    !! First and second order mapping functions for [-1, 1] with infinity at 1

    !! Note: Cubatures for tetrahedra use body coordinates, and the condition
    !!       x1 + x2 + x3 + x4 = 1 allows the elimination of the final
    !!       coordinate (likewise with triangles)

    character(3) :: elmtype = ""    !! Element type
    integer :: dime = 0         !! Dimension
    integer :: order = 0        !! Order
    integer :: nodes = 0        !! Nodes
    integer :: points = 0       !! Cubature points
    character(3) :: inf(3) = ["FIN", "FIN", "FIN"]    !! Infinitude
    real(rk), allocatable :: func(:,:)      !! Shape function values (nodes,points)
    real(rk), allocatable :: derv(:,:,:)    !! Shape function first derivatives (dimension,nodes,points)
    real(rk), allocatable :: curv(:,:,:,:)  !! Shape function second derivatives (dimension,dimension,nodes,points)

    ! Internally useful

    ! Function pointers
    procedure(lagrpoly), pointer, nopass :: lagrange1, lagrange2, lagrange3

    integer, allocatable :: seq(:,:)    !! Nodal numbering sequence

    integer :: ounit = stdout
  contains

    ! Fulfill object requirements
    procedure :: check
    procedure :: summary
    procedure :: show
    procedure :: numbering
    procedure :: destroy

    ! Unique procedures
    procedure :: set
    procedure :: eval
    final :: destroy_final
  end type Shapefunc

  public :: Shapefunc

contains

!***********************************************************************

pure function check(self)
!! Check integrity of instance

  class(shapefunc), intent(in) :: self
  logical :: check

  check = (allocated(self%func) .neqv. allocated(self%derv))

end function check

!***********************************************************************

subroutine summary(self)
!! Summary contents to file

  class(shapefunc), intent(in) :: self

  associate (u => self%ounit)
    write(u,"(A20,A)") "Element type: ", self%elmtype
    write(u,"(A20,I0)") "Element dimension: ", self%dime
    write(u,"(A20,I0)") "Element orders: ", self%order
    write(u,"(A20,I0)") "Element nodes: ", self%nodes
    write(u,"(A20,I0)") "Element points: ", self%points
    write(u,"(A20,3(A3,:,', '))") "Element infinitude: ", self%inf
  end associate

end subroutine summary

!***********************************************************************

subroutine show(self)
!! Print shape function values and derivatives to terminal

  class(Shapefunc), intent(in) :: self

  integer :: g, n

  character(11) :: fmt1 = "*(x,F20.15)"

  ! Information
  write(self%ounit,"(A)") "SHAPE FUNCTION DATA"
  write(self%ounit,"(A20,A)")  "Element type: "   , self%elmtype
  write(self%ounit,"(A20,I0)") "Nodes: "          , self%nodes
  write(self%ounit,"(A20,I0)") "Cubature points: ", self%points

  if (any(self%inf /= "FIN")) then
    do n = 1,self%dime
      select case (self%inf(n))
      case ("FIN")
        write(self%ounit,"(A,I0,A)") "Direction ",n," is finite"
      case ("INF")
        write(self%ounit,"(A,I0,A)") "Direction ",n," is infinite"
      case ("CHP")
        write(self%ounit,"(A,I0,A)") "Direction ",n," is chopped"
      end select
    end do ! n
  else
    write(self%ounit,"(A)") "All directions have standard domain"
  end if

  select case (self%elmtype)

  ! Information on hexahedra and quadrilaterals
  case ("QUA","HEX")
    write(self%ounit,"(A,I1,A)") "Point, node, func, derv (",self%dime,")"
    do g = 1,self%points
      do n = 1,self%nodes
        write(self%ounit,"(2I4,"//fmt1//")") g, n, self%func(n,g), self%derv(:,n,g)
      end do ! n
      write(self%ounit,"(A8,"//fmt1//")") "Check: ", sum(self%func(  :,g)), &
                                            sum(self%derv(:,:,g),2)
    end do ! g

  ! Information on tetrahedra and triangles
  case ("TRI","TET")
    write(self%ounit,"(A,I1,A)") "Point, node, func, derv (",self%dime,")"
    do g = 1,self%points
      do n = 1,self%nodes
        write(self%ounit,"(2I4,"//fmt1//")") g, n, self%func(n,g), self%derv(:,n,g)
      end do ! n
      write(self%ounit,"(A8,"//fmt1//")") "Check: ", sum(self%func(:,g)), &
                                            sum(self%derv(:,:,g),2)
    end do ! g

  ! Information on wedges
  case ("WEJ")
    write(self%ounit,"(A,I1,A)") "Point, node, func, derv (",self%dime,")"
    do g = 1,self%points
      do n = 1,self%nodes
        write(self%ounit,"(2I4,"//fmt1//")") g, n, self%func(n,g), self%derv(:,n,g)
      end do ! n
      write(self%ounit,"(A8,"//fmt1//")") "Check: ", sum(self%func(:,g)), &
                                            sum(self%derv(:,:,g),2)
    end do ! g

  end select

  write(self%ounit,*)

end subroutine show

!***********************************************************************

subroutine numbering(self)
!! Print node numbers and coordinates

  class(shapefunc), intent(in) :: self

  integer :: i, j, k
  integer :: m, n

  real(rk) :: coords(self%dime,self%nodes)

  character(*), parameter :: fmt1 = "(I3,') ',*(F8.4,:,', '))"

  select case (self%elmtype)

  case ("LIN")

    ! Assemble coordinates
    do i = 1,self%order+1
      n = counter(i,1,1,[self%order+1,1,1])
      m = self%seq(1,n)
      coords(1,m) = -1.0_rk + 2.0_rk*(i-1)/self%order
    end do ! j

    ! Print coordinates
    do i = 1,self%nodes
      write(self%ounit,fmt1) i, coords(:,i)
    end do ! i

  case ("QUA")

    ! Assemble coordinates
    do i = 1,self%order+1
      do j = 1,self%order+1
        n = counter(i,j,1,[self%order+1,self%order+1,1])
        m = self%seq(1,n)
        coords(1,m) = -1.0_rk + 2.0_rk*(j-1)/self%order
        coords(2,m) = -1.0_rk + 2.0_rk*(i-1)/self%order
      end do ! j
    end do ! j

    ! Print coordinates
    do i = 1,self%nodes
      write(self%ounit,fmt1) i, coords(:,i)
    end do ! i

  case ("TRI")

    ! Assemble coordinates
    do i = 1,self%nodes
      coords(1,i) = 1.0_rk*(self%seq(1,i)-1)/self%order
      coords(2,i) = 1.0_rk*(self%seq(2,i)-1)/self%order
    end do ! j

    ! Print coordinates
    do i = 1,self%nodes
      write(self%ounit,fmt1) i, coords(:,i), 1-sum(coords(:,i))
    end do ! i

  case ("HEX")

    ! Assemble coordinates
    do i = 1,self%order+1
      do j = 1,self%order+1
        do k = 1,self%order+1
          n = counter(i,j,k,[self%order+1,self%order+1,self%order+1])
          m = self%seq(1,n)
          coords(1,m) = -1.0_rk + 2.0_rk*(k-1)/self%order
          coords(2,m) = -1.0_rk + 2.0_rk*(j-1)/self%order
          coords(3,m) = -1.0_rk + 2.0_rk*(i-1)/self%order
        end do ! k
      end do ! j
    end do ! j

    ! Print coordinates
    do i = 1,self%nodes
      write(self%ounit,fmt1) i, coords(:,i)
    end do ! i

  case ("TET")

    ! Assemble coordinates
    do i = 1,self%nodes
      coords(1,i) = 1.0_rk*(self%seq(1,i)-1)/self%order
      coords(2,i) = 1.0_rk*(self%seq(2,i)-1)/self%order
      coords(3,i) = 1.0_rk*(self%seq(3,i)-1)/self%order
    end do ! j

    ! Print coordinates
    do i = 1,self%nodes
      write(self%ounit,fmt1) i, coords(:,i), 1-sum(coords(:,i))
    end do ! i

  case ("WEJ")

    write(self%ounit,"(A)") "TODO"

  end select

end subroutine numbering

!***********************************************************************

pure subroutine destroy(self)
!! Reset type to initial state

  class(Shapefunc), intent(inout) :: self

  self%elmtype = ""
  self%order   = 0
  self%nodes   = 0
  self%points  = 0
  self%inf     = "FIN"
  if (allocated(self%func)) deallocate(self%func)
  if (allocated(self%derv)) deallocate(self%derv)
  if (allocated(self%curv)) deallocate(self%curv)
  if (allocated(self%seq )) deallocate(self%seq )

end subroutine destroy

!***********************************************************************

pure subroutine destroy_final(self)
!! Destroy derived type
  type(Shapefunc), intent(inout) :: self
  call self%destroy()
end subroutine destroy_final

!***********************************************************************

pure subroutine set(self,cube,order,infin)
!! Set shape functions given a cubature and a shape function code
!! - cubature and shape functions must be compatible

  class(Shapefunc), intent(inout) :: self
  class(Cubature), intent(in) :: cube             !! Cubature derived type
  integer, intent(in) :: order                !! Shape function order
  character(3), intent(in), optional :: infin(:)  !! Infinitude or chopitude along each direction

  integer :: i
  integer :: mem
  integer, parameter :: tri(4) = [3,  6, 10, 15]
  integer, parameter :: tet(4) = [4, 10, 20, 35]
  integer, parameter :: wej(4) = [6, 18, 40, 75]

  ! Wipe any existing data
  if (self%elmtype /= "") then
    call self%destroy()
  end if

  ! Read inputs and prepare memory

  ! Read for infinite domains
  if (present(infin)) then
    self%inf = "FIN"
    self%inf(1:size(infin)) = infin
  else
    self%inf = "FIN"
  end if

  ! Set function pointers to default
  self%lagrange1 => lagrfull
  self%lagrange2 => lagrfull
  self%lagrange3 => lagrfull

  ! Set parameters
  self%elmtype = cube%elmtype
  self%dime    = cube%dime
  self%order   = order
  mem = 0

  ! Set particular pointers
  select case (self%elmtype)
  case ("HEX","QUA","LIN")
    self%nodes = (self%order+1)**self%dime
    if ("INF" == self%inf(1)) self%lagrange1 => lagrinf
    if ("INF" == self%inf(2)) self%lagrange2 => lagrinf
    if ("INF" == self%inf(3)) self%lagrange3 => lagrinf
    mem = 1
  case ("TET")
    self%nodes = tet(self%order)
    self%lagrange1 => lagrhalf
    self%lagrange2 => lagrhalf
    self%lagrange3 => lagrhalf
    mem = 4
  case ("TRI")
    self%nodes = tri(self%order)
    self%lagrange1 => lagrhalf
    self%lagrange2 => lagrhalf
    self%lagrange3 => lagrhalf
    mem = 3
  case ("WEJ")
    self%nodes = wej(self%order)
    self%lagrange1 => lagrhalf
    self%lagrange2 => lagrhalf
    mem = 4
  end select

  ! Take cubature points
  self%points = cube%points

  ! Allocate function and derivative memory, and working sequence memory
  allocate(self%func(                        1:self%nodes,1:self%points))
  allocate(self%derv(            1:self%dime,1:self%nodes,1:self%points))
  allocate(self%curv(1:self%dime,1:self%dime,1:self%nodes,1:self%points))
  allocate(self%seq(mem,self%nodes))

  !!! Set order and sequences

  select case (self%elmtype)

  ! Line
  case ("LIN")

    select case (self%order)

    case (0)
      self%seq(1,:) = [1]

    case (1)
      self%seq(1,:) = [1, 2]

    case (2)
      self%seq(1,:) = [1, 3, 2]

    case (3)
      self%seq(1,:) = [1, 3, 4, 2]

    case (4)
      self%seq(1,:) = [1, 3, 4, 5, 2]

    end select

  ! Quadrilateral
  case ("QUA")

    select case (self%order)

    case (0)
      self%seq(1,:) = [1]

    case (1)
      self%seq(1,:) = [1, 2, 4, 3]

    case (2)
      self%seq(1,:) = [1, 5, 2, 8, 9, 6, 4, 7, 3]

    case (3)
      self%seq(1,:) = [1, 5, 6, 2, 12, 13, 14, 7, 11, 16, 15, 8, 4, 10, 9, 3]

    case (4)
      self%seq(1,:) = [(i, i=1,25)]

    case default
      error stop "Shapefunc%set: Invalid order for quadrilateral"
    end select

  ! Triangle
  case ("TRI")

    select case (self%order)

    case (1)
      self%seq(:,1) = [1, 0, 0]
      self%seq(:,2) = [0, 1 ,0]
      self%seq(:,3) = [0, 0, 1]
      self%seq = self%seq + 1

    case (2)

      ! Vertices
      self%seq(:, 1) = [2, 0, 0]
      self%seq(:, 2) = [0, 2, 0]
      self%seq(:, 3) = [0, 0, 2]

      ! Edges
      self%seq(:, 4) = [1, 1, 0]
      self%seq(:, 5) = [0, 1, 1]
      self%seq(:, 6) = [1, 0, 1]

      self%seq = self%seq + 1

    case (3)

      ! Vertices
      self%seq(:, 1) = [3, 0, 0]
      self%seq(:, 2) = [0, 3, 0]
      self%seq(:, 3) = [0, 0, 3]

      ! Edges
      self%seq(:, 4) = [2, 1, 0]
      self%seq(:, 5) = [1, 2, 0]
      self%seq(:, 6) = [0, 2, 1]
      self%seq(:, 7) = [0, 1, 2]
      self%seq(:, 8) = [1, 0, 2]
      self%seq(:, 9) = [2, 0, 1]

      ! Face
      self%seq(:,10) = [1, 1, 1]

      self%seq = self%seq + 1

    case default
      error stop "Shapefunc%set: Invalid order for triangle"
    end select

  ! Hexahedron
  case ("HEX")

    select case (self%order)

    case (0)
      self%seq(1,:) = [1]

    ! Linear full Lagrange polynomials - 2x2x2 = 8 nodes
    case (1)
      self%seq(1,:) = [1, 2, 4, 3, 5, 6, 8, 7]

    ! Quadratic full Lagrange polynomials - 3x3x3 = 27 nodes
    case (2)
      self%seq(1,:) = [ 1,  9,  2, 12, 21, 10,  4, 11,  3, &
                       17, 25, 18, 23, 27, 24, 20, 26, 19, &
                        5, 13,  6, 16, 22, 14,  8, 15,  7]

    ! Cubic full Lagrange polynomials - 4x4x4 = 64 nodes
    case (3)
      self%seq(1,:) = [(i, i=1,64)]

    ! Quartic full Lagrange polynomials - 5x5x5 = 125 nodes
    case (4)
      self%seq(1,:) = [(i, i=1,125)]

    case default
      error stop "Shapefunc%set: Invalid order for hexahedron"
    end select

  ! Tetrahedron
  case ("TET")

    select case (self%order)

    ! Linear half Lagrange polynomials - 4 nodes
    case (1)
      self%seq(:,1) = [1, 0, 0, 0]
      self%seq(:,2) = [0, 1, 0, 0]
      self%seq(:,3) = [0, 0, 1, 0]
      self%seq(:,4) = [0, 0, 0, 1]

      self%seq = self%seq + 1

    ! Quadratic half Lagrange polynomials - 10 nodes
    case (2)

      ! Vertices
      self%seq(:, 1) = [2, 0, 0, 0]
      self%seq(:, 2) = [0, 2, 0, 0]
      self%seq(:, 3) = [0, 0, 2, 0]
      self%seq(:, 4) = [0, 0, 0, 2]

      ! Edges
      self%seq(:, 5) = [1, 1, 0, 0]
      self%seq(:, 6) = [0, 1, 1, 0]
      self%seq(:, 7) = [1, 0, 1, 0]
      self%seq(:, 8) = [1, 0, 0, 1]
      self%seq(:, 9) = [0, 1, 0, 1]
      self%seq(:,10) = [0, 0, 1, 1]

      self%seq = self%seq + 1

    ! Cubic half Lagrange polynomials - 10 nodes
    case (3)

      ! Vertices
      self%seq(:, 1) = [3, 0, 0, 0]
      self%seq(:, 2) = [0, 3, 0, 0]
      self%seq(:, 3) = [0, 0, 3, 0]
      self%seq(:, 4) = [0, 0, 0, 3]

      ! Edges
      self%seq(:, 5) = [2, 1, 0, 0]
      self%seq(:, 6) = [1, 2, 0, 0]
      self%seq(:, 7) = [0, 2, 1, 0]
      self%seq(:, 8) = [0, 1, 2, 0]
      self%seq(:, 9) = [1, 0, 2, 0]
      self%seq(:,10) = [2, 0, 1, 0]
      self%seq(:,11) = [2, 0, 0, 1]
      self%seq(:,12) = [1, 0, 0, 2]
      self%seq(:,13) = [0, 2, 0, 1]
      self%seq(:,14) = [0, 1, 0, 2]
      self%seq(:,15) = [0, 0, 2, 1]
      self%seq(:,16) = [0, 0, 1, 2]

      ! Faces
      self%seq(:,17) = [1, 1, 1, 0]
      self%seq(:,18) = [1, 1, 0, 1]
      self%seq(:,19) = [0, 1, 1, 1]
      self%seq(:,20) = [1, 0, 1, 1]

      self%seq = self%seq + 1

    case default
      error stop "Shapefunc%set: Invalid order for tetrahedron"
    end select

  ! Wedge
  case ("WEJ")

    select case (self%order)

    ! Linear polynomials - 6 nodes
    case (1)

      ! Bottom
      self%seq(:,1) = [1, 0, 0, 0]
      self%seq(:,2) = [0, 1, 0, 0]
      self%seq(:,3) = [0, 0, 1, 0]

      ! Top
      self%seq(:,4) = [1, 0, 0, 1]
      self%seq(:,5) = [0, 1, 0, 1]
      self%seq(:,6) = [0, 0, 1, 1]

      self%seq = self%seq + 1

    ! Quadratic polynomials - 18 nodes
    case (2)

      ! Bottom
      self%seq(:,  1) = [2, 0, 0, 0]
      self%seq(:,  2) = [0, 2, 0, 0]
      self%seq(:,  3) = [0, 0, 2, 0]
      self%seq(:,  4) = [2, 0, 0, 2]
      self%seq(:,  5) = [0, 2, 0, 2]
      self%seq(:,  6) = [0, 0, 2, 2]

      ! Top
      self%seq(:,  7) = [1, 1, 0, 0]
      self%seq(:,  8) = [0, 1, 1, 0]
      self%seq(:,  9) = [1, 0, 1, 0]
      self%seq(:, 10) = [1, 1, 0, 2]
      self%seq(:, 11) = [0, 1, 1, 2]
      self%seq(:, 12) = [1, 0, 1, 2]

      ! Middle
      self%seq(:, 13) = [2, 0, 0, 1]
      self%seq(:, 14) = [0, 2, 0, 1]
      self%seq(:, 15) = [0, 0, 2, 1]
      self%seq(:, 16) = [1, 1, 0, 1]
      self%seq(:, 17) = [0, 1, 1, 1]
      self%seq(:, 18) = [1, 0, 1, 1]

      self%seq = self%seq + 1

    case default
      error stop "Shapefunc%set: Invalid order for wedge"
    end select
  end select

  !!! Calculate shape functions and derivatives

  ! Loop over cubature points
  do i = 1,self%points
    call self%eval(cube%abscissae(:,i), &
      self%func(:,i),self%derv(:,:,i),self%curv(:,:,:,i))
  end do ! i

end subroutine set

!***********************************************************************

pure subroutine eval(self,cord, func,derv,curv)
!! Evaluate shape functions and derivatives at a point

  class(Shapefunc), intent(in) :: self
  real(rk), intent(in) :: cord(self%dime)                       !! Point at which to evaluate
  real(rk), intent(out) :: func(self%nodes)                     !! Shape function values
  real(rk), intent(out) :: derv(self%dime,self%nodes)           !! Shape function derivatives
  real(rk), intent(out) :: curv(self%dime,self%dime,self%nodes) !! Shape function second derivatives

  integer :: i, j, k, m, n
  integer :: a(4)

  real(rk) :: r
  real(rk) :: f(5,4) ! Directional shape function values      (node,coordinate)
  real(rk) :: d(5,4) ! Directional shape function derivatives (node,coordinate)
  real(rk) :: c(5,4) ! Directional shape function second derivatives

  real(rk) :: coord(3)

  ! Copy abscissae
  coord(3) = 0
  coord(1:self%dime) = cord(1:self%dime)

  ! Calculate lagrange polynomials and derivatives at each abscissa and for each node
  f = 0
  d = 0
  c = 0
  if ("CHP" == self%inf(1)) then
    call self%lagrange1(self%order+1,coord(1), &
      f(1:self%order+2,1),d(1:self%order+2,1),c(1:self%order+2,1))
  else
    call self%lagrange1(self%order  ,coord(1), &
      f(1:self%order+1,1),d(1:self%order+1,1),c(1:self%order+1,1))
  end if
  if (self%dime >= 2) then
    if ("CHP" == self%inf(2)) then
      call self%lagrange2(self%order+1,coord(2), &
        f(1:self%order+2,2),d(1:self%order+2,2),c(1:self%order+2,2))
    else
      call self%lagrange2(self%order  ,coord(2), &
        f(1:self%order+1,2),d(1:self%order+1,2),c(1:self%order+1,2))
    end if
  end if
  if (self%dime == 3) then
    if ("CHP" == self%inf(3)) then
      call self%lagrange3(self%order+1,coord(3), &
        f(1:self%order+2,3),d(1:self%order+2,3),c(1:self%order+2,3))
    else
      call self%lagrange3(self%order  ,coord(3), &
        f(1:self%order+1,3),d(1:self%order+1,3),c(1:self%order+1,3))
    end if
  end if

  select case (self%elmtype)

  ! Calculate lineal shape functions and derivatives
  case ("LIN")

    do i = 1,self%order+1
      n = counter(i,1,1,[self%order+1,1,1])
      m = self%seq(1,n)
      func(m)     = f(i,1)
      derv(1,m)   = d(i,1)
      curv(1,1,m) = c(i,1)
    end do ! i, j

  ! Calculate quadrilateral shape functions and derivatives
  case ("QUA")

    do i = 1,self%order+1
      do j = 1,self%order+1
        n = counter(i,j,1,[self%order+1,self%order+1,1])
        m = self%seq(1,n)
        func(m)   = f(j,1)*f(i,2)
        derv(1,m) = d(j,1)*f(i,2)
        derv(2,m) = f(j,1)*d(i,2)
        curv(1,1,m) = c(j,1)*f(i,2)
        curv(2,1,m) = d(j,1)*d(i,2)
        curv(1,2,m) = curv(2,1,m)
        curv(2,2,m) = f(j,1)*c(i,2)
      end do ! j
    end do ! i

  ! Calculate triangle shape functions
  case ("TRI")

    r = 1 - coord(1) - coord(2)
    call self%lagrange1(self%order,r, f(:,3),d(:,3),c(:,3))

    do m = 1,self%nodes
      a(1:3) = self%seq(:,m)
      func(m)   = f(a(1),1)*f(a(2),2)*f(a(3),3)
      derv(1,m) = d(a(1),1)*f(a(2),2)*f(a(3),3) - &
                  f(a(1),1)*f(a(2),2)*d(a(3),3)
      derv(2,m) = f(a(1),1)*d(a(2),2)*f(a(3),3) - &
                  f(a(1),1)*f(a(2),2)*d(a(3),3)
      curv(1,1,m) =   c(a(1),1)*f(a(2),2)*f(a(3),3) - &
                    2*d(a(1),1)*f(a(2),2)*d(a(3),3) + &
                      f(a(1),1)*f(a(2),2)*c(a(3),3)
      curv(2,1,m) =   d(a(1),1)*d(a(2),2)*f(a(3),3) - &
                      d(a(1),1)*f(a(2),2)*d(a(3),3) - &
                      f(a(1),1)*d(a(2),2)*d(a(3),3) + &
                      f(a(1),1)*f(a(2),2)*c(a(3),3)
      curv(1,2,m) = curv(2,1,m)
      curv(2,2,m) =   f(a(1),1)*c(a(2),2)*f(a(3),3) - &
                    2*f(a(1),1)*d(a(2),2)*d(a(3),3) + &
                      f(a(1),1)*f(a(2),2)*c(a(3),3)
    end do ! i

  ! Calculate hexahedron shape functions
  case ("HEX")

    do i = 1,self%order+1
      do j = 1,self%order+1
        do k = 1,self%order+1
          n = counter(i,j,k,[self%order+1,self%order+1,self%order+1])
          m = self%seq(1,n)
          func(m)   = f(k,1)*f(j,2)*f(i,3)
          derv(1,m) = d(k,1)*f(j,2)*f(i,3)
          derv(2,m) = f(k,1)*d(j,2)*f(i,3)
          derv(3,m) = f(k,1)*f(j,2)*d(i,3)
          curv(1,1,m) = c(k,1)*f(j,2)*f(i,3)
          curv(2,1,m) = d(k,1)*d(j,2)*f(i,3)
          curv(3,1,m) = d(k,1)*f(j,2)*d(i,3)
          curv(1,2,m) = curv(2,1,m)
          curv(2,2,m) = f(k,1)*c(j,2)*f(i,3)
          curv(3,2,m) = f(k,1)*d(j,2)*d(i,3)
          curv(1,3,m) = curv(3,1,m)
          curv(2,3,m) = curv(3,2,m)
          curv(3,3,m) = f(k,1)*f(j,2)*c(i,3)
        end do ! k
      end do ! j
    end do ! i

  ! Calculate tetrahedron shape functions
  ! Note the chain rule for handling the fourth body coordinate
  case ("TET")

    r = 1 - coord(1) - coord(2) - coord(3)
    call self%lagrange1(self%order,r, f(:,4),d(:,4),c(:,4))

    do m = 1,self%nodes
      a = self%seq(:,m)
      func(m)   = f(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4)

      derv(1,m) = f(a(1),1)*d(a(2),2)*f(a(3),3)*f(a(4),4) - &
                  d(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4)
      derv(2,m) = f(a(1),1)*f(a(2),2)*d(a(3),3)*f(a(4),4) - &
                  d(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4)
      derv(3,m) = f(a(1),1)*f(a(2),2)*f(a(3),3)*d(a(4),4) - &
                  d(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4)

!      curv(1,1,m) =   c(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4) - &
!                    2*d(a(1),1)*f(a(2),2)*f(a(3),3)*d(a(4),4) + &
!                      f(a(1),1)*f(a(2),2)*f(a(3),3)*c(a(4),4)
!      curv(2,1,m) =   d(a(1),1)*d(a(2),2)*f(a(3),3)*f(a(4),4) - &
!                      d(a(1),1)*f(a(2),2)*f(a(3),3)*d(a(4),4) - &
!                      f(a(1),1)*d(a(2),2)*f(a(3),3)*d(a(4),4) + &
!                      f(a(1),1)*f(a(2),2)*f(a(3),3)*c(a(4),4)
!      curv(3,1,m) =   d(a(1),1)*f(a(2),2)*d(a(3),3)*f(a(4),4) - &
!                      d(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4) - &
!                      f(a(1),1)*f(a(2),2)*d(a(3),3)*d(a(4),4) + &
!                      f(a(1),1)*f(a(2),2)*f(a(3),3)*c(a(4),4)
!      curv(1,2,m) =   curv(2,1,m)
!      curv(2,2,m) =   f(a(1),1)*c(a(2),2)*f(a(3),3)*f(a(4),4) - &
!                    2*f(a(1),1)*d(a(2),2)*f(a(3),3)*f(a(4),4) + &
!                      f(a(1),1)*f(a(2),2)*f(a(3),3)*c(a(4),4)
!      curv(3,2,m) =   f(a(1),1)*d(a(2),2)*d(a(3),3)*f(a(4),4) - &
!                      f(a(1),1)*d(a(2),2)*f(a(3),3)*d(a(4),4) - &
!                      f(a(1),1)*f(a(2),2)*d(a(3),3)*d(a(4),4) + &
!                      f(a(1),1)*f(a(2),2)*f(a(3),3)*c(a(4),4)
!      curv(1,3,m) =   curv(3,1,m)
!      curv(2,3,m) =   curv(3,2,m)
!      curv(3,3,m) =   f(a(1),1)*f(a(2),2)*c(a(3),3)*f(a(4),4) - &
!                    2*f(a(1),1)*f(a(2),2)*d(a(3),3)*d(a(4),4) + &
!                      f(a(1),1)*f(a(2),2)*f(a(3),3)*c(a(4),4)
    end do !

    curv = 0

  ! Calculate wedge shape functions
  case ("WEJ")

    r = 1 - coord(1) - coord(2)
    f(:,4) = f(:,3)
    d(:,4) = d(:,3)
    c(:,4) = c(:,3)
    call self%lagrange1(self%order,r, f(:,3),d(:,3),c(:,3))

    do m = 1,self%nodes
      a(1:4) = self%seq(:,m)
      func(m)   = f(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4)
      derv(1,m) = d(a(1),1)*f(a(2),2)*f(a(3),3)*f(a(4),4) - &
                  f(a(1),1)*f(a(2),2)*d(a(3),3)*f(a(4),4)
      derv(2,m) = f(a(1),1)*d(a(2),2)*f(a(3),3)*f(a(4),4) - &
                  f(a(1),1)*f(a(2),2)*d(a(3),3)*f(a(4),4)
      derv(3,m) = f(a(1),1)*f(a(2),2)*f(a(3),3)*d(a(4),4)
    end do ! m

    curv = 0

  end select

end subroutine eval

!***********************************************************************

pure subroutine lagrfull(order,x, func,derv,curv)
!! Provide lagrangian polynomials and derivatives for the domain [-1, 1]

  integer, intent(in) :: order            !! Polynomial order
  real(rk), intent(in) :: x                   !! Point at which to evaluate
  real(rk), intent(out) :: func(abs(order)+1) !! Polynomial values
  real(rk), intent(out) :: derv(abs(order)+1) !! Polynomial derivatives
  real(rk), intent(out) :: curv(abs(order)+1) !! Polynomial second derivatives

  select case (order)

  case (0)
    func = [1.0_rk]
    derv = [0.0_rk]
    curv = [0.0_rk]

  case (1)
    func = [   -x+1,    x+1]/2
    derv = [-1.0_rk, 1.0_rk]/2
    curv = [ 0.0_rk, 0.0_rk]

  case (2)
    func = [x*(x-1)/2, 1.0_rk-x**2, x*(x+1)/2]
    derv = [    x-0.5,        -2*x,     x+0.5]
    curv = [   1.0_rk,     -2.0_rk,    1.0_rk]

  case (3)
    func = [        -(3*x+1)*(3*x-1)*(x-1),&
             9*(x+1)*        (3*x-1)*(x-1),&
            -9*(x+1)*(3*x+1)*        (x-1),&
               (x+1)*(3*x+1)*(3*x-1)      ]/16

    derv = [   -27*x**2+18*x+1 ,&
            9*(  9*x**2- 2*x-3),&
            9*( -9*x**2- 2*x+3),&
                27*x**2+18*x-1]/16

    curv = [   -54*x+18,&
               182*x-18,&
              -182*x-18,&
                54*x+18]/16

  case (4)
    func = [         (2*x+1)*x*(2*x-1)*(x-1)/6,&
            -4*(x+1)*        x*(2*x-1)*(x-1)/3,&
                 (x+1)*(2*x+1)*  (2*x-1)*(x-1)  ,&
            -4*(x+1)*(2*x+1)*x*        (x-1)/3,&
                 (x+1)*(2*x+1)*x*(2*x-1)      /6]

    derv = [   (16*x**3-12*x**2- 2*x+1)/6,&
            -4*( 8*x**3- 6*x**2- 4*x+1)/3,&
                16*x**3        -10*x     ,&
            -4*( 8*x**3+ 6*x**2- 4*x-1)/3,&
               (16*x**3+12*x**2- 2*x-1)/6]

    curv = [( 24*x**2-12*x- 1)/3,&
            (-96*x**2+24*x+16)/3,&
            ( 48*x**2     -10)  ,&
            (-96*x**2-24*x+16)/3,&
            ( 48*x**2+24*x- 2)/6]

  case default
    error stop "shapefuncs: lagrfull invalid order"
  end select

end subroutine lagrfull

!***********************************************************************

pure subroutine lagrhalf(order,x, func,derv,curv)
!! Provide lagrangian polynomials and derivatives for the domain [0, 1]

  integer, intent(in) :: order            !! Polynomial order
  real(rk), intent(in) :: x                   !! Point at which to evaluate
  real(rk), intent(out) :: func(abs(order)+1) !! Polynomial values
  real(rk), intent(out) :: derv(abs(order)+1) !! Polynomial derivatives
  real(rk), intent(out) :: curv(abs(order)+1) !! Polynomial second derivatives

  select case (order)

  case (1)
    func = [1.0_rk,      x]
    derv = [0.0_rk, 1.0_rk]
    curv = [0.0_rk, 0.0_rk]

  case (2)
    func = [1.0_rk,   2* x, x*(2*x-1)]
    derv = [0.0_rk, 2.0_rk,    4*x-1 ]
    curv = [0.0_rk, 0.0_rk,    4.0_rk]

  case (3)
    func = [1.0_rk,    3*x, 3*x*(3*x-1)/2, x*(3*x-1)*(3*x-2)/2]
    derv = [0.0_rk, 3.0_rk,       9*x-1.5,     27*x**2/2-9*x+1]
    curv = [0.0_rk, 0.0_rk,        9.0_rk,              27*x-9]

  case (4)
    func = [1.0_rk,   4* x, 2*x*(4*x-1), 4*x*(4*x-1)*(2*x-1)/3, x*(4*x-1)*(2*x-1)*(4*x-3)/3]
    derv = [0.0_rk, 4.0_rk,      16*x-2, 32*x**2-16*x+4.0_rk/3, 128*x**3/3+48*x**2+44*x/3-1]
    curv = [0.0_rk, 0.0_rk,     16.0_rk,               64*x-16,     128*x**2+48*x+44.0_rk/3]

  end select

end subroutine lagrhalf

!***********************************************************************

pure subroutine lagrinf(order,x, func,derv,curv)
!! Provide lagrangian polynomials and derivatives for the domain [-1, 1],
!! for infinity at 1

  integer, intent(in) :: order            !! Polynomial order
  real(rk), intent(in) :: x                   !! Point at which to evaluate
  real(rk), intent(out) :: func(abs(order)+1) !! Polynomial values
  real(rk), intent(out) :: derv(abs(order)+1) !! Polynomial derivatives
  real(rk), intent(out) :: curv(abs(order)+1) !! Polynomial second derivatives

  select case (order)

  ! To complement linear shape functions
  case (1)
    func = [   -2*x,    1+x]/(1-x)
    derv = [-2.0_rk, 2.0_rk]/(1-x)**2
    curv = [-2.0_rk, 2.0_rk]/(1-x)**3

  ! To complement quadratic shape functions
  case (2)
    func = [(3*x+1)*(3*x-1)/4, -(x+1)*(3*x-1),  (x+1)*(3*x+1)/4]/(1-x)
    derv = [(-9*x*x+18*x-1)/4,    3*x*x-6*x-1, (-3*x*x+6*x+5)/4]/(1-x)**2
    curv = [           4.0_rk,        -8.0_rk,           4.0_rk]/(1-x)**3

  ! To complement cubic shape functions
  case (3)
    func = [-2*      (2*x+1)*x*(2*x-1)/3, &
             3*(x+1)*        x*(2*x-1)  , &
            -  (x+1)*(2*x+1)*  (2*x-1)  , &
               (x+1)*(2*x+1)*x        /3]/(1-x)
    derv = [ 2*( 8*x**3-12*x*x    +1)/3, &
            -3*( 4*x**3- 5*x*x-2*x+1),   &
                 8*x**3- 8*x*x-8*x+2,    &
               (-4*x**3+ 3*x*x+6*x+1)/3]/(1-x)**2
    curv = [(-16*x**3+48*x**2-48*x+4)/3, &
            ( 12*x**3-36*x**2+36*x   ), &
            (- 8*x**3+24*x**2-24*x- 4), &
            (  4*x**3-12*x**2+12*x+ 8)/3]/(1-x)**3

  ! To complement quartic shape functions
  case (4)
    func = [        (5*x+3)*(5*x+1)*(5*x-1)*(5*x-3)/192, &
             -(x+1)*        (5*x+1)*(5*x-1)*(5*x-3)/12,  &
            3*(x+1)*(5*x+3)*        (5*x-1)*(5*x-3)/32,  &
             -(x+1)*(5*x+3)*(5*x+1)*        (5*x-3)/24,  &
              (x+1)*(5*x+3)*(5*x+1)*(5*x-1)        /192]/(1-x)
    derv = [  (-1875*x**4+2500*x**3+250*x*x-500*x+ 9)/192, &
              (  375*x**4- 400*x**3-230*x*x+160*x- 1)/12, &
            3*(- 375*x**4+ 300*x**3+370*x*x-140*x-27)/32, &
              (  375*x**4- 200*x**3-470*x*x+ 40*x+63)/24, &
              (- 375*x**4+ 100*x**3+530*x*x+140*x-11)/192]/(1-x)**2
    curv = [( 1875*x**4-5000*x**3+3750*x**2      -241)/96, &
            (- 375*x**4+ 950*x**3- 600*x**2-150*x+ 79)/6, &
            ( 1125*x**4-2700*x**3+1350*x**2+900*x-291)/16, &
            (- 375*x**4+ 850*x**3- 300*x**2-450*x+ 83)/12, &
            (  375*x**4- 800*x**3+ 150*x**2+600*x+ 59)/96]/(1-x)**3

  case default
    error stop "shapefuncs: lagrinf invalid order"
  end select

end subroutine lagrinf

!***********************************************************************

end module shapefuncs
