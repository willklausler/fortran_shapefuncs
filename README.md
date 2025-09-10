# Fortran shape functions

## About

Interpolation function derived type for standard domains (LaGrange) and infinite domains for finite elements in

- 1D
  - line elements - standard and infinite+chopped domains,
- 2D
  - triangle elements - standard domain,
  - quadrilateral elements - standard and infinite+chopped domain,
- 3D
  - tetrahedron elements - standard domain,
  - hexahedron elements - standard and infinite+chopped domain,
  - wedge/prism elements - standard domain,

providing the functions as well as first and second derivatives with respect to the natural coordinates.

This derived type depends on the derived type 'cubatures' found in the repository 'fortran cubatures'. See the documentation for details on usage.

## Infinite domain finite elements

Usually, one interpolates the coordinate and solution fields of a finite element with identical shape functions. To interpolate over an infinite domain, one substitutes coordinate shape functions with so-called mapping functions so that the isoparametric domain [-1, 1] maps to [x(-1), infinity]. The provided mapping functions fulfill partition of unity and are suitable for isotropic nodal configurations. The corresponding shape functions for the solution field are modified to accommodate the assumption of a solution value of 0 at the infinite coordinate. These modified shape functions are the function of one greater order, and chopped to achieve nodal isotropy.

## Usage

```fortran
! Declare type
type(cubature) :: scheme
type(shapefunc) :: shapes
character(3) :: geo = "HEX"
integer :: order = 1

! Set cubature
call scheme%set(geo, order+1)

! Set shape functions
! Set domain with "FIN" = finite (default), "INF" = infinite, "CHP" = chopped
call shapes%set(scheme, order, infin=["FIN", "FIN", "INF"])

! Access number of integration points
print *, shapes%points

! Access number of nodes
print *, shapes%nodes

! Access dimension
print *, shapes%dime

! Access shape functions, shape = [nodes, points]
print *, shapes%func

! Access shape function derivatives, shape = [dime, nodes, points]
print *, shapes%derv

! Access shape function second derivatives, shape = [dime, dime, nodes, points]
print *, shapes%curv

! Display nodal coordinates in isoparametric configuration
call shapes%numbering()

! Display shape function values and derivatives at each node and integration point
call shapes%show()
```

## To do

- Infinite domain wedges

- Unit testing
