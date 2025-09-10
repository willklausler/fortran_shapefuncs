program shapefuncs_example
!! Demonstrate shape functions

  use cubatures
  use shapefuncs

  implicit none

  type(Cubature)  :: cuba
  type(Shapefunc) :: shaper

  call cuba%set("HEX",[2])

  call shaper%set(cuba,2)

  call shaper%show()

  call shaper%numbering()

end program shapefuncs_example