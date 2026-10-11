# fortran_shapefuncs

[![CI](https://github.com/willklausler/fortran_shapefuncs/actions/workflows/ci.yml/badge.svg)](https://github.com/willklausler/fortran_shapefuncs/actions/workflows/ci.yml)
[![Docs](https://github.com/willklausler/fortran_shapefuncs/actions/workflows/docs.yml/badge.svg)](https://willklausler.github.io/fortran_shapefuncs/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Shape functions for the finite element method: one derived type, `shapefunc`,
that holds the shape functions of an element and their first and second
derivatives at the points of a cubature from
[fortran_cubatures](https://github.com/willklausler/fortran_cubatures).

```fortran
use cubatures
use shapefuncs
type(cubature) :: q
type(shapefunc) :: s

q = cubature(CUB_HEX, 4)                     ! integration points
s = shapefunc(q, 2)                          ! 27-node hexahedron at those points
print *, s%nnodes, sum(s%func(:,1))          ! 27  1.0: a partition of unity

s = shapefunc(q, 2, family=SHP_SERENDIPITY)
print *, s%nnodes                            ! 20
s = shapefunc(cubature(CUB_LIN, 2), 1, [SHP_INF])
print *, s%coords                            ! -1.0  0.0, infinite at xi = 1
```

> Coding agents: the complete API reference is [llms.txt](llms.txt), and
> contributor instructions are in [AGENTS.md](AGENTS.md).

## Features

- Lines, triangles, quadrilaterals, tetrahedra, hexahedra and wedges (prisms).
- Lagrange elements of any order, with values, first and second derivatives.
- Serendipity quadrilaterals and hexahedra: 8-, 12-, 20- and 32-node elements.
- Mapped infinite elements for unbounded domains on lines, quadrilaterals,
  hexahedra, and wedges along their axis.
- Classical node numbering (vertices, edges, faces, interior), or any
  numbering given at run time.
- Pure procedures, no module state.

Conventions: the reference elements of fortran_cubatures, `s%func(node, point)`
with derivatives in the leading dimensions, and invalid input stops the program.

## Installation

With [fpm](https://fpm.fortran-lang.org) 0.13 or later, add to your `fpm.toml`:

```toml
[dependencies]
fortran_shapefuncs = { git = "https://github.com/willklausler/fortran_shapefuncs", tag = "v0.2.0" }
```

It pulls in fortran_cubatures v0.2.0.

With make, which clones fortran_cubatures v0.2.0 and builds
`build/make/libfortran_shapefuncs.a` (both modules) and the module files:

```sh
make FC=gfortran && make test
make install PREFIX=$HOME/.local    # library, .mod files and man page
```

## Documentation

- [User manual and API documentation](https://willklausler.github.io/fortran_shapefuncs/), generated with FORD (`make doc`)
- [llms.txt](llms.txt): the API reference as one plain-text file, also the manual's source
- `man 3 shapefuncs` after `make install` ([man/shapefuncs.3](man/shapefuncs.3))
- Example programs in [example/](example/) (`fpm run --example --all`) and an app in [app/](app/) (`fpm run`) that tabulates the interpolation error against order for every element type

## Compilers

CI tests every push on:

| Compiler | Linux | macOS | Windows |
| --- | --- | --- | --- |
| gfortran 14 | yes | yes | yes |
| gfortran 13 | yes | | |
| ifx 2026 | yes | | yes |
| LLVM flang 21 | yes | | |
| nvfortran 26 | yes | | |

There are no known compiler limitations; see [llms.txt](llms.txt#compiler-notes).

## Why another shape-function package?

No other Fortran package gives finite element shape functions as a library.
[FEM_BASIS](https://people.sc.fsu.edu/~jburkardt/f_src/fem_basis/fem_basis.html)
evaluates Lagrange basis functions of any degree, but only on simplices and
with no node numbering for meshes.
[FEM2D_PACK](https://people.sc.fsu.edu/~jburkardt/f_src/fem2d_pack/fem2d_pack.html)
has low-order 2D triangles and quadrilaterals only.
[FELIPE](https://people.bath.ac.uk/mbr20/felipe/) is a Fortran 77 teaching
package whose programs include mapped infinite elements in 1D and 2D, built
into the solvers. Solvers such as [Elmer](https://github.com/ElmerCSC/elmerfem)
keep their shape functions inside. This package covers all six FE elements of
any order, serendipity and mapped infinite elements, through one derived type
that works with fortran_cubatures. Every element is verified by
function-space reproduction tests, and it is packaged for fpm and tested with
gfortran, ifx, LLVM flang and nvfortran on Linux, macOS and Windows.

## References

1. P. Bettess, Infinite elements, *International Journal for Numerical Methods in Engineering* 11(1), 53–64, 1977. [doi:10.1002/nme.1620110107](https://doi.org/10.1002/nme.1620110107)
2. O. C. Zienkiewicz, C. Emson, P. Bettess, A novel boundary infinite element, *International Journal for Numerical Methods in Engineering* 19(3), 393–404, 1983. [doi:10.1002/nme.1620190307](https://doi.org/10.1002/nme.1620190307)
3. J. M. M. C. Marques, D. R. J. Owen, Infinite elements in quasi-static materially nonlinear problems, *Computers & Structures* 18(4), 739–751, 1984. [doi:10.1016/0045-7949(84)90019-1](https://doi.org/10.1016/0045-7949(84)90019-1)
4. P. Bettess, *Infinite Elements*, Penshaw Press, 1992.
5. P. Silvester, High-order polynomial triangular finite elements for potential problems, *International Journal of Engineering Science* 7(8), 849–861, 1969. [doi:10.1016/0020-7225(69)90065-2](https://doi.org/10.1016/0020-7225(69)90065-2)
6. O. C. Zienkiewicz, R. L. Taylor, J. Z. Zhu, *The Finite Element Method: Its Basis and Fundamentals*, 6th ed., Butterworth-Heinemann, 2005.
7. D. N. Arnold, G. Awanou, The serendipity family of finite elements, *Foundations of Computational Mathematics* 11(3), 337–344, 2011. [doi:10.1007/s10208-011-9087-3](https://doi.org/10.1007/s10208-011-9087-3)

The finite elements are the standard Lagrange and serendipity elements [6, 7],
with Silvester's polynomials for simplices [5]. The infinite elements are those
of Bettess [1] and Zienkiewicz, Emson and Bettess [2], generalised to higher
orders [3, 4]. The code has not been independently reviewed.

## Status and support

Version 0.2.0. The API may still change before 1.0; changes are recorded in
[CHANGELOG.md](CHANGELOG.md). Bug reports, questions and pull requests are
welcome on [GitHub](https://github.com/willklausler/fortran_shapefuncs/issues);
see [CONTRIBUTING.md](CONTRIBUTING.md). To cite this package, use
[CITATION.cff](CITATION.cff).

## License

[MIT](LICENSE)
