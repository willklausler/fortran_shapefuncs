---
title: Getting started
---

## Installation

### With fpm

Add the package to the `fpm.toml` of your project:

```toml
[dependencies]
fortran_shapefuncs = { git = "https://github.com/willklausler/fortran_shapefuncs" }
```

fpm fetches the `fortran_cubatures` dependency automatically. To install the
library and module files instead, run `fpm install --prefix <dir>`.

### With make

```sh
make                              # clones fortran_cubatures, builds the library
make test                         # runs the unit tests
make install PREFIX=/usr/local    # library, modules and pkg-config file
```

The library `libfortran_shapefuncs.a` bundles the `cubatures` module, so you
link only this one library:

```sh
gfortran main.f90 $(pkg-config --cflags --libs fortran_shapefuncs)
```

To build against a local copy of `fortran_cubatures`, pass its source file:
`make CUBATURES=../fortran_cubatures/src/cubatures.f90`. Other variables are
`FC`, `FFLAGS`, `PREFIX` and `DESTDIR` (for staged installs when packaging).
Module files are compiler-specific, so build and install with the compiler
that compiles your code.

### By hand

Compile `src/cubatures.f90` from `fortran_cubatures`, then `src/shapefuncs.f90`.

## A first program

```fortran
program first
  use cubatures, only: cubature, rk, CUB_QUA
  use shapefuncs, only: shapefunc

  implicit none

  type(cubature) :: q
  type(shapefunc) :: s

  q = cubature(CUB_QUA, 4)     ! Gauss rule exact for degree 4 per direction
  s = shapefunc(q, 2)          ! 9-node quadratic quadrilateral

  call s%summary()             ! element, order, nodes, points
  call s%numbering()           ! nodal coordinates

  print *, s%func(:,1)         ! N_i at the first point
  print *, s%derv(:,:,1)       ! dN_i/dxi_j at the first point, shape [2, 9]
end program first
```

The arrays are:

| Component | Shape | Contents |
| --- | --- | --- |
| `s%func` | `[nnodes, npoints]` | \(N_i(\xi_g)\) |
| `s%derv` | `[dim, nnodes, npoints]` | \(\partial N_i/\partial\xi_j\) at \(\xi_g\) |
| `s%curv` | `[dim, dim, nnodes, npoints]` | \(\partial^2 N_i/\partial\xi_j\partial\xi_k\) at \(\xi_g\) |
| `s%coords` | `[dim, nnodes]` | nodal coordinates in the reference element |

`s%eval(xi, func, derv, curv)` evaluates the same quantities at any other
point, for example for post-processing.

## The element loop

With nodal coordinates `xn(dim, nnodes)` of a physical element, the
isoparametric map gives the Jacobian \(J = \partial x/\partial\xi\) and the
physical gradients \(\partial N/\partial x = J^{-T}\,\partial N/\partial\xi\).
A stiffness matrix for the Laplacian of a quadrilateral element:

```fortran
q = cubature(CUB_QUA, 2*p)
s = shapefunc(q, p)
k = 0
do g = 1, q%npoints
  jac  = matmul(xn, transpose(s%derv(:,:,g)))         ! dx/dxi, [2, 2]
  detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
  jinv = reshape([jac(2,2), -jac(2,1), -jac(1,2), jac(1,1)], [2, 2])/detj
  dndx = matmul(transpose(jinv), s%derv(:,:,g))       ! dN/dx, [2, nnodes]
  k = k + matmul(transpose(dndx), dndx)*detj*q%weights(g)
end do
```

The shape functions depend only on the reference element, so build each
`shapefunc` once per element type and reuse it for every element of the mesh.
A code with mixed meshes can keep a table indexed by the element constants:

```fortran
type(shapefunc) :: table(6)
table(CUB_TET) = shapefunc(cubature(CUB_TET, 4), 2)
table(CUB_HEX) = shapefunc(cubature(CUB_HEX, 4), 2)
```

## Choosing the cubature

For an affine element of order \(p\), the mass matrix needs a rule of degree
\(2p\) and the stiffness matrix \(2p-2\). Curved and distorted elements need
more. Infinite elements are not polynomial: see
[Infinite elements](infinite-elements.html).

## Running the example

`example/shapefuncs_example.f90` prints a node numbering, computes the area of
a curved quadrilateral and integrates \(1/x^2\) over \([1,\infty)\) with one
infinite element:

```sh
fpm run --example     # or: make example
```
