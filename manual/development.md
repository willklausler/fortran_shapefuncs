---
title: Development
---

## Building and testing

```sh
fpm build
fpm test --profile debug     # bounds checks, sanitizers, FP traps
fpm test --profile release
fpm run --example
make test                    # the same tests without fpm
ford ford.md                 # API documentation and this manual in docs/
```

The `debug` and `release` profiles in `fpm.toml` set flags for gfortran, Intel
ifx and NVIDIA nvfortran.

## What the tests check

`test/shapefuncs_test.f90` sweeps every element type over a range of orders,
and every combination of finite, infinite and chopped directions. For each
case it checks that the shape functions

- are consistent (`is_valid`) and have the expected numbers of nodes and points;
- are 1 at their own node and 0 at every other node;
- reproduce every function of their space exactly at every cubature point,
  including its gradient and Hessian. For finite directions the space is the
  polynomials, for mapping functions \(q(\xi)/(1-\xi)\), and for chopped
  functions \((1-\xi)\,q(\xi)\). Together with the nodal check, this proves
  that every function and derivative is correct;
- have derivatives that agree with fourth-order finite differences, and
  symmetric second derivatives;
- sum to 1 with derivatives summing to 0, or vanish at infinity when chopped.

It also checks the classical numbering of orders 1 to 3 node by node, closed
forms from the literature, integrals over each element, the constructor, `set`,
reuse, `eval`, `destroy`, order 0 and the output routines.

## Continuous integration

GitHub Actions (`.github/workflows/ci.yml`) runs on every push and pull request:

- `fpm test` in both profiles and the example, with gfortran 13–15, Intel ifx
  and NVIDIA nvfortran on Linux, and gfortran on macOS and Windows;
- `make test`, `make example` and a staged `make install`, then a program
  compiled against the installed library;
- a FORD documentation build.

## Design

Every shape function is a product of one-dimensional functions of affine
*factor coordinates* (the tensor coordinates, or the barycentric
coordinates of simplices), so one routine evaluates all element types.
Derivatives follow from the product rule and a constant chain-rule matrix.
The one-dimensional functions are built as products of linear or rational
factors, which gives exact derivatives for any order without closed-form
polynomials. The node numbering is generated recursively (vertices, edges,
faces, interior), so any order is supported.

## Pedigree and support

The finite shape functions are the standard Lagrange elements of finite element
textbooks. The infinite elements follow Zienkiewicz, Emson and Bettess and
Marques and Owen (see [Infinite elements](infinite-elements.html)). The linear
mapping and chopped functions are checked against the published closed forms.
Every function is verified by the tests described above. The code has not been
independently reviewed. It is maintained by the author on a best-effort basis.

Please report bugs and request features through
[GitHub issues](https://github.com/willklausler/fortran_shapefuncs/issues).
Pull requests are welcome; please include a test.
