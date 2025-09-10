module cubatures
!! Define cubature derived type

  use iso_fortran_env, only: rk => real64, stdout => output_unit

  implicit none

  private

  character(3), parameter :: elmtypes(*) = &
    [ "LIN",   "TRI",   "QUA",    "TET",  "HEX",  "WEJ"]

  real(rk), parameter :: volumes(*) = &
    [2.0_rk, 1.0_rk/2, 4.0_rk, 1.0_rk/6, 8.0_rk, 1.0_rk]

  integer, parameter :: maxorders(*) = &
    [     5,        3,      5,        3,      5,      3]

  public :: elmtypes, volumes, maxorders

  type :: Cubature
    !! This derived type holds data, weights, and abscissae for numerical
    !! integration of 2D (triangle, quadrilaterals) and 3D bodies (tetrahedrons,
    !! hexahedrons, and wedges)
    !!
    !! Note: Cubatures for triangles and tetrahedrons use body coordinates,
    !! and the condition x1 + x2 + x3 (+ x4) = 1 permits the elimination of
    !! the final coordinate. For wedges, the first two coordinates are body,
    !! the last is ordinary.
    !!
    !! Domains:
    !! LINear - [-1, 1]
    !! QUAdrilateral - [-1, 1] x [-1, 1]
    !! TRIangle - [0, 1] x [0, 1-x]
    !! HEXahedron - [-1, 1] x [-1, 1] x [-1, 1]
    !! TETrahedron - [0, 1] x [0, 1-x] x [0, 1-x-y]
    !! WEJ/prism - [0, 1] x [0, 1-x] x [-1, 1]

    character(3) :: elmtype = ""              !! Element type: TRI, QUA, TET, HEX, WEJ
    integer :: dime = 0                       !! Element dimension
    integer :: orders(3) = 0                  !! Order of cubature
    integer :: points = 0                     !! Number of cubature points
    real(rk), allocatable :: abscissae(:,:)   !! Abscissae (coordinates)
    real(rk), allocatable :: weights(:)       !! Weights

    integer :: ounit = stdout

  contains

    ! Fulfill object requirements
    procedure :: check
    procedure :: summary
    procedure :: show
    procedure :: destroy

    ! Unique procedures
    procedure :: set
    final :: destroy_final
  end type Cubature

  public :: Cubature

  public :: counter

contains

!***********************************************************************

pure function check(self)
!! Check integrity of derived type

  class(cubature), intent(in) :: self
  logical :: check

  check = (allocated(self%abscissae) .neqv. allocated(self%weights))

end function check

!***********************************************************************

subroutine summary(self)
!! Summarize cubature scheme

  class(cubature), intent(in) :: self

  associate (u => self%ounit)
    write(u,"(A20,A)") "Element type: ", self%elmtype
    write(u,"(A20,I0)") "Element dimension: ", self%dime
    write(u,"(A20,3(I1,:,','))") "Element orders: ", self%orders
    write(u,"(A20,I0)") "Element points: ", self%points
  end associate

end subroutine summary

!***********************************************************************

subroutine show(self)
!! Print cubature scheme to output

  class(Cubature), intent(in) :: self

  integer :: i

  character(*), parameter :: fmt0 = "*(G0,:,', ')"
  character(*), parameter :: fmt1 = "(I3,')',1x,"//fmt0//")"
  character(*), parameter :: fmt2 = "('Chk:',1x,"//fmt0//")"

  write(self%ounit,"(A)") "CUBATURE DATA"

  select case (self%elmtype)

  ! Information on lines
  case ("LIN")
    write(self%ounit,"(A,1I2,1x,A)") "LIN", self%orders(1), "Abscissae (1), Weight"
    do i = 1,self%points
      write(self%ounit,fmt1) i, self%abscissae(:,i), self%weights(i)
    end do ! i
    write(self%ounit,fmt2) sum(self%abscissae(1,:)), sum(self%weights)

  ! Information on quadrilaterals
  case ("QUA")
    write(self%ounit,"(A,2I2,1x,A)") "QUA", self%orders(1:2), "Abscissae (2), Weight"
    do i = 1,self%points
      write(self%ounit,fmt1) i, self%abscissae(:,i), self%weights(i)
    end do ! i
    write(self%ounit,fmt2) (sum(self%abscissae(i,:)), i=1,2), sum(self%weights)

  ! Information on triangles
  case ("TRI")
    write(self%ounit,"(A,1I2,1x,A)") "TRI", self%orders(1), "Abscissae (2), Coord 3, Weight"
    do i = 1,self%points
      write(self%ounit,fmt1) i, &
                             self%abscissae(:,i), &
                             1-sum(self%abscissae(:,i)), self%weights(i)
    end do ! i
    write(self%ounit,fmt2) (sum(self%abscissae(i,:))/self%points, i=1,2), &
                           1-sum(self%abscissae(:,:))/self%points, &
                           sum(self%weights)

  ! Information on hexahedra
  case ("HEX")
    write(self%ounit,"(A,3I2,1x,A)") "HEX", self%orders, "Abscissae (3), Weight"
    do i = 1,self%points
      write(self%ounit,fmt1) i, self%abscissae(:,i), self%weights(i)
    end do ! i
    write(self%ounit,fmt2) (sum(self%abscissae(i,:)), i=1,3), sum(self%weights)

  ! Information on tetrahedra
  case ("TET")
    write(self%ounit,"(A,1I2,1x,A)") "TET", self%orders(1), "Abscissae (3), Coord 4, Weight"
    do i = 1,self%points
      write(self%ounit,fmt1) i, &
                             self%abscissae(:,i), &
                             1-sum(self%abscissae(:,i)), self%weights(i)
    end do ! i
    write(self%ounit,fmt2) (sum(self%abscissae(i,:))/self%points, i=1,3), &
                           1-sum(self%abscissae(:,:))/self%points, &
                           sum(self%weights)

  ! Information on wedges
  case ("WEJ")
    write(self%ounit,"(A,1I2,1x,A)") "WEJ", self%orders(1), "Abscissae (1-2), Coord 3, Abscissa 3, Weight"
    do i = 1,self%points
      write(self%ounit,fmt1) i, &
                             self%abscissae(1:2,i), &
                             1-sum(self%abscissae(1:2,i)), &
                             self%abscissae(3,i), &
                             self%weights(i)
    end do ! i
    write(self%ounit,fmt2)   sum(self%abscissae(1  ,:))/self%points, &
                             sum(self%abscissae(2  ,:))/self%points, &
                           1-sum(self%abscissae(1:2,:))/self%points, &
                             sum(self%abscissae(3  ,:))/self%points, &
                             sum(self%weights)

  end select

write(self%ounit,*)

end subroutine show

!***********************************************************************

pure subroutine destroy(self)
!! wipe data and memory

  class(Cubature), intent(inout) :: self

  self%elmtype = ""
  self%orders  = 0
  self%points  = 0
  if (allocated(self%abscissae)) deallocate(self%abscissae)
  if (allocated(self%weights  )) deallocate(self%weights)
  self%ounit = stdout

end subroutine destroy

!***********************************************************************

pure subroutine destroy_final(self)
!! Destroy derived type
  type(Cubature), intent(inout) :: self
  call self%destroy()
end subroutine destroy_final

!***********************************************************************

pure subroutine set(self, elmtype,order)
!! Set cubature scheme by inputting shape and order

  class(Cubature), intent(inout) :: self
  character(3), intent(in) :: elmtype       !! Element type: TRI, QUA, TET, HEX, WEJ
  integer, intent(in) :: order(:)    !! Cubature order: size 1 or dimension

  integer :: i, j, k     !! Iterators
  integer :: n           !! Index of abscissae

  ! Wipe any existing data
  call self%destroy()

  ! Set data
  self%elmtype = elmtype
  select case (self%elmtype)
  case ("LIN")
    self%dime = 1
  case ("TRI","QUA")
    self%dime = 2
  case ("HEX","TET","WEJ")
    self%dime = 3
  end select

  ! Read order
  select case (size(order,1))
  case (1)
    self%orders(1) = order(1)
    if (self%orders(1) > 0) then
      self%orders(1:self%dime) = order(1)
      self%orders(self%dime+1:3) = 1
    end if
  case (2)
    if ((self%elmtype == "HEX") .or. (self%elmtype == "TET") .or. (self%elmtype == "WEJ")) then
      error stop "cubatures%set: Invalid order array for 3D"
    end if
    if (self%elmtype == "LIN") then
      error stop "cubatures%set: Invalid order array for 1D"
    end if
    self%orders(1:2) = order
    self%orders(3)   = 1
  case (3)
    if ((self%elmtype == "QUA") .or. (self%elmtype == "TRI")) then
      error stop "cubatures%set: Invalid order array for 2D"
    end if
    if (self%elmtype == "LIN") then
      error stop "cubatures%set: Invalid order array for 1D"
    end if
    self%orders = order
  end select

  select case (self%elmtype)

  case ("LIN")

    self%points = self%orders(1)

    ! Allocate memory
    allocate(self%abscissae(1:1,1:self%points))
    allocate(self%weights(1:self%points))

    block
      real(rk) :: abscissae(maxval(self%orders))
      real(rk) :: weights(maxval(self%orders))

      ! Retrieve lineature
      call gauss(self%orders(1), abscissae,weights)

      ! Accumulate cubature
      do i = 1, self%orders(1)
        n = counter(i,1,1,self%orders)
        self%abscissae(1,n) = abscissae(i)
        self%weights(n)     = weights(i)
      end do ! i

    end block

  case ("QUA")

    self%points = self%orders(1)*self%orders(2)

    ! Allocate memory
    allocate(self%abscissae(1:2,1:self%points))
    allocate(self%weights(1:self%points))

    block
      real(rk) :: abscissae(maxval(self%orders),2)
      real(rk) :: weights(maxval(self%orders),2)

      ! Retrieve lineature
      call gauss(self%orders(1), abscissae(:,1),weights(:,1))
      call gauss(self%orders(2), abscissae(:,2),weights(:,2))

      ! Accumulate cubature
      do i = 1,self%orders(1)
        do j = 1,self%orders(2)
          n = counter(i,j,1,self%orders)
          self%abscissae(:,n) = [abscissae(i,1), abscissae(j,2)]
          self%weights(n)     =    weights(i,1)*   weights(j,2)
        end do ! j
      end do ! i

    end block

  case ("TRI")

    select case (self%orders(1))

    ! 1-Point triangle quadrature
    case (1)

      self%points = 1

      allocate(self%abscissae(1:2,1:self%points))
      allocate(self%weights(1:self%points))

      self%abscissae(:,1) = 1.0_rk/3
      self%weights(1) = 1.0_rk/2

    ! 3-Point triangle quadrature - interior (see Zienkiewicz & Taylor)
    case (2)

      self%points = 3

      allocate(self%abscissae(1:2,1:self%points))
      allocate(self%weights(1:self%points))

      self%abscissae      = 1.0_rk/6
      self%abscissae(1,1) = 2.0_rk/3
      self%abscissae(2,2) = 2.0_rk/3

      self%weights = 1.0_rk/6

    ! 4-Point triangle quadrature - interior
    case (3)

      self%points = 4

      allocate(self%abscissae(1:2,1:self%points))
      allocate(self%weights(1:self%points))

      self%abscissae(:,1)   = 1.0_rk/3
      self%abscissae(:,2:4) = 0.2_rk
      self%abscissae(1,2)   = 0.6_rk
      self%abscissae(2,3)   = 0.6_rk

      self%weights(1)   = -27.0_rk/96
      self%weights(2:4) =  25.0_rk/96

    case default
      error stop "Cubature%set: Invalid order for triangle"
    end select

  case ("HEX")

    if (self%orders(2) > 0) then

      self%points = product(self%orders)

      ! Allocate memory
      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      block
        real(rk) :: abscissae(maxval(self%orders),3)
        real(rk) :: weights(maxval(self%orders),3)

        ! Retrieve lineature
        do i = 1,3
          call gauss(self%orders(i), abscissae(:,i),weights(:,i))
        end do ! i

        ! Accumulate cubature
        do i = 1, self%orders(1)
          do j = 1, self%orders(2)
            do k = 1, self%orders(3)
              n = counter(i,j,k,self%orders)
              self%abscissae(:,n) = [abscissae(i,1), abscissae(j,2), abscissae(k,3)]
              self%weights(n)     =    weights(i,1)*   weights(j,2)*   weights(k,3)
            end do ! k
          end do ! j
        end do ! i

      end block

    ! 4-Point hexahedron cubature
    else if (self%orders(1) == -4) then

      self%points = 4

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      block
        real(rk), parameter :: sqrt13 = sqrt(1.0_rk/3)

        self%abscissae(:,1) = [-sqrt13, -sqrt13, -sqrt13]
        self%abscissae(:,2) = [ sqrt13,  sqrt13, -sqrt13]
        self%abscissae(:,3) = [ sqrt13, -sqrt13,  sqrt13]
        self%abscissae(:,4) = [-sqrt13,  sqrt13,  sqrt13]
        self%weights = 2
      end block

    ! 9-Point hexahedron cubature
    else if (self%orders(1) == -9) then

      self%points = 9

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      block
        real(rk), parameter :: sqrt35 = sqrt(3.0_rk/5)
        real(rk), parameter :: tick(1:2) = [-sqrt35, sqrt35]

        self%weights(1:8) = 5.0_rk/9
        self%weights(9)   = 32.0_rk/9

        do i = 1,2
          do j = 1,2
            do k = 1,2
              n = counter(i,j,k,[2, 2, 2])
              self%abscissae(:,n) = [tick(k), tick(j), tick(i)]
            end do ! k
          end do ! j
        end do ! i

        self%abscissae(:,9) = 0
      end block

    else
      error stop "Cubature%set: Invalid order for hexahedron"
    end if

  case ("TET")

    select case (self%orders(1))

    ! 1-Point tetrahedron cubature
    case (1)

      self%points = 1

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      self%abscissae(:,1) = 0.25_rk
      self%weights(1)     = 1.0_rk/6

    ! 4-Point tetrahedron cubature - interior (see Zienkiewicz & Taylor)
    case (2)

      self%points = 4

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      block
        real(rk), parameter :: alpha = 0.5854101966249685_rk
        real(rk), parameter ::  beta = 0.1381966011250105_rk

        self%abscissae = beta
        do i = 1,3
          self%abscissae(i,i) = alpha
        end do ! i
      end block

      self%weights = 1.0_rk/24

    ! 5-Point tetrahedron cubature - interior
    case (3)

      self%points = 5

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      self%abscissae(:,1:4) = 1.0_rk/6
      self%abscissae(1,1)   = 0.5_rk
      self%abscissae(2,2)   = 0.5_rk
      self%abscissae(3,3)   = 0.5_rk
      self%abscissae(:,5)   = 0.25_rk

      self%weights(1:4)     =  0.45_rk/6
      self%weights(5)       = -0.8_rk/6

    case default
      error stop "Cubature%set: Invalid order for tetrahedron"
    end select

  case ("WEJ")

    select case (self%orders(1))

    ! 2-Point wedge cubature
    case (1)

      self%points = 2

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      block
        real(rk) :: abscissae(2)
        real(rk) :: weights(2)

        call gauss(2, abscissae,weights)

        self%abscissae(1:2,:) = 1.0_rk/3
        self%abscissae(  3,:) = abscissae

        self%weights = 0.5_rk*weights

      end block

    ! 9-Point wedge cubature
    case (2)

      self%points = 9

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      block
        real(rk) :: abscissae(3)
        real(rk) :: weights(3)

        call gauss(3, abscissae,weights)

        self%abscissae(1:2,:) = 1.0_rk/6
        self%abscissae(1,1) = 2.0_rk/3
        self%abscissae(2,2) = 2.0_rk/3
        self%abscissae(1,4) = 2.0_rk/3
        self%abscissae(2,5) = 2.0_rk/3
        self%abscissae(1,7) = 2.0_rk/3
        self%abscissae(2,8) = 2.0_rk/3

        self%abscissae(3,1:3) = abscissae(1)
        self%abscissae(3,4:6) = abscissae(2)
        self%abscissae(3,7:9) = abscissae(3)

        self%weights(1:3) = 1.0_rk/6*weights(1)
        self%weights(4:6) = 1.0_rk/6*weights(2)
        self%weights(7:9) = 1.0_rk/6*weights(3)

      end block

    ! 16-Point wedge cubature
    case (3)

      self%points = 16

      allocate(self%abscissae(1:3,1:self%points))
      allocate(self%weights(1:self%points))

      block
        real(rk) :: abscissae(4)
        real(rk) :: weights(4)

        call gauss(4, abscissae,weights)

        do i = 1,4
          self%abscissae(1:2,4*(i-1)+1) = 1.0_rk/3
          self%abscissae(1:2,4*(i-1)+2) = [0.6_rk, 0.2_rk]
          self%abscissae(1:2,4*(i-1)+3) = [0.2_rk, 0.6_rk]
          self%abscissae(1:2,4*(i-1)+4) = 0.2_rk
        end do ! i

        self%abscissae(3, 1: 4) = abscissae(1)
        self%abscissae(3, 5: 8) = abscissae(2)
        self%abscissae(3, 9:12) = abscissae(3)
        self%abscissae(3,13:16) = abscissae(4)

        do i = 1,4
          self%weights(4*(i-1)+1) = -27.0_rk/96*weights(i)
          self%weights(4*(i-1)+2) =  25.0_rk/96*weights(i)
          self%weights(4*(i-1)+3) =  25.0_rk/96*weights(i)
          self%weights(4*(i-1)+4) =  25.0_rk/96*weights(i)
        end do ! i

      end block

    case default
      error stop "Cubature%set: Invalid order for wedge"
    end select

  case default
    error stop "Cubature%set: Invalid element type"
  end select

end subroutine set

!***********************************************************************

pure function counter(i,j,k,r) result(c)
!! Count in base r - indexing starts at 1

  integer, intent(in) :: i, j, k, r(3)
  integer :: c

  c = r(2)*r(3)*(i-1) + r(3)*(j-1) + (k-1) + 1

end function counter

!***********************************************************************

pure subroutine gauss(order, abscissae,weights)
!! Provide Gaussian lineature abscissae and weights given the order, for domain [-1, 1]

  integer, intent(in) :: order           !! Order of integration
  real(rk), intent(out) :: abscissae(order)  !! Abscissae (coordinates)
  real(rk), intent(out) :: weights(order)    !! Weights

  select case (order)

  ! 1 point
  case (1)
    abscissae(1) = 0
    weights(1)   = 2

  ! 2 point
  case (2)
    block
      real(rk), parameter :: sqrt13 = sqrt(1.0_rk/3)

      abscissae = [-sqrt13, sqrt13]
      weights   = [ 1.0_rk, 1.0_rk]
    end block

  ! 3 point
  case (3)
    block
      real(rk), parameter :: sqrt35 = sqrt(0.6_rk)
      real(rk), parameter :: frac59 = 5.0_rk/9
      real(rk), parameter :: frac89 = 8.0_rk/9

      abscissae = [-sqrt35, 0.0_rk, sqrt35]
      weights   = [ frac59, frac89, frac59]
    end block

  ! 4 point
  case (4)
    block
      real(rk), parameter :: sr4d8 = sqrt(4.8_rk)
      real(rk), parameter :: sr30  = sqrt(30.0_rk)
      real(rk), parameter :: root1 = sqrt((3.0_rk - sr4d8)/7)
      real(rk), parameter :: root2 = sqrt((3.0_rk + sr4d8)/7)
      real(rk), parameter :: frac1 = 0.5_rk + sr30/36
      real(rk), parameter :: frac2 = 0.5_rk - sr30/36

      abscissae = [-root2, -root1, root1, root2]
      weights   = [ frac2,  frac1, frac1, frac2]
    end block

  ! 5 point
  case (5)
    block
      real(rk), parameter :: sr107 = sqrt(10.0_rk/7)
      real(rk), parameter :: sr70  = sqrt(70.0_rk)
      real(rk), parameter :: root1 = sqrt(5.0_rk-2*sr107)/3
      real(rk), parameter :: root2 = sqrt(5.0_rk+2*sr107)/3
      real(rk), parameter :: frac1 = (322.0_rk+13*sr70)/900
      real(rk), parameter :: frac2 = (322.0_rk-13*sr70)/900

      abscissae = [-root2, -root1,       0.0_rk, root1, root2]
      weights   = [ frac2,  frac1, 128.0_rk/225, frac1, frac2]
    end block
  end select

end subroutine gauss

!***********************************************************************

end module cubatures
