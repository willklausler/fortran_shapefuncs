# fortran_shapefuncs

[![CI](https://github.com/willklausler/fortran_shapefuncs/actions/workflows/ci.yml/badge.svg)](https://github.com/willklausler/fortran_shapefuncs/actions/workflows/ci.yml)

Shape functions for the finite element method, packaged as one derived type,
`shapefunc`, that works with the cubatures of
[fortran_cubatures](https://github.com/willklausler/fortran_cubatures).

- Lines, triangles, quadrilaterals, tetrahedra, hexahedra and wedges (prisms)
- Lagrange elements of any order, with values, first and second derivatives
- Serendipity quadrilaterals and hexahedra: 8-, 12-, 20- and 32-node elements
- Mapped infinite elements for unbounded domains: lines, quadrilaterals,
  hexahedra, and wedges along their axis
- Classical node numbering (vertices, edges, faces, interior), or any
  numbering given at run time
- Pure procedures, Fortran 2018, fpm and make builds

## Installation

Add the package to your `fpm.toml`:

```toml
[dependencies]
fortran_shapefuncs = { git = "https://github.com/willklausler/fortran_shapefuncs" }
```

Without fpm, `make && make install PREFIX=...` builds and installs the library,
the module files and a pkg-config file. See
[Getting started](manual/getting-started.md).

## Quick start

```fortran
use cubatures,  only: cubature, rk, CUB_HEX
use shapefuncs, only: shapefunc, SHP_FIN, SHP_INF, SHP_CHP

type(cubature) :: q
type(shapefunc) :: s

q = cubature(CUB_HEX, 4)     ! integration points
s = shapefunc(q, 2)          ! 27-node quadratic hexahedron at those points

! s%func(i,g)       N_i at point g
! s%derv(:,i,g)     dN_i/dxi
! s%curv(:,:,i,g)   d2N_i/dxi2
```

In an element loop, the Jacobian is `matmul(xnodes, transpose(s%derv(:,:,g)))`.
Infinite directions are flagged per direction. Use mapping functions for the
geometry and chopped functions for the solution, which vanishes at infinity:

```fortran
geom  = shapefunc(q, 1, [SHP_FIN, SHP_FIN, SHP_INF])   ! infinite at xi_3 = 1
field = shapefunc(q, 1, [SHP_FIN, SHP_FIN, SHP_CHP])
```

Serendipity elements and a caller-defined node ordering are optional arguments:

```fortran
s = shapefunc(q, 2, family=SHP_SERENDIPITY)              ! 20-node hexahedron
s = shapefunc(q, 2, family=SHP_SERENDIPITY, nodes=perm)  ! your node k is node perm(k)
```

The runnable example is [example/shapefuncs_example.f90](example/shapefuncs_example.f90):
`fpm run --example`.

## API

| Entity | Description |
| --- | --- |
| `shapefunc(q, order[, infin, family, nodes])` | Constructor. `q` is a `cubature`, `infin` has size 1 or `q%dim`, `nodes` is a permutation |
| `call s%set(q, order[, infin, family, nodes])` | Build in place, reusing the object |
| `call s%eval(xi, func, derv, curv)` | Evaluate at any point of the reference element |
| `s%elm`, `s%dim`, `s%order`, `s%family` | Element type (`CUB_*`), dimension, order, family |
| `s%nnodes`, `s%npoints` | Numbers of nodes and cubature points |
| `s%infin(3)` | Infinitude per direction, `SHP_*` |
| `s%func(nnodes, npoints)` | Values |
| `s%derv(dim, nnodes, npoints)` | First derivatives |
| `s%curv(dim, dim, nnodes, npoints)` | Second derivatives |
| `s%coords(dim, nnodes)` | Nodal coordinates |
| `s%lattice(:, nnodes)` | Integer lattice position of each node |
| `s%is_valid()` | `.true.` if set and consistent |
| `call s%summary([unit])`, `show`, `numbering` | Print a summary, values, or nodes |
| `call s%destroy()` | Deallocate and reset |
| `SHP_FIN`, `SHP_INF`, `SHP_CHP` | Finite, infinite (mapping), infinite (chopped) |
| `SHP_LAGRANGE`, `SHP_SERENDIPITY` | Element family |

Invalid input stops with `error stop` and a message. The
[API reference](manual/api.md) lists every condition.

## Elements

| Constant | Element | Nodes, order p | Orders | Infinite directions |
| --- | --- | --- | --- | --- |
| `CUB_LIN` | line | p+1 | ≥ 0 | 1 |
| `CUB_QUA` | quadrilateral | (p+1)² | ≥ 0 | any |
| `CUB_HEX` | hexahedron | (p+1)³ | ≥ 0 | any |
| `CUB_TRI` | triangle | (p+1)(p+2)/2 | ≥ 1 | none |
| `CUB_TET` | tetrahedron | (p+1)(p+2)(p+3)/6 | ≥ 1 | none |
| `CUB_WED` | wedge | (p+1)²(p+2)/2 | ≥ 1 | axial (3) |
| `CUB_QUA`, serendipity | quadrilateral | 4p | 1–3 | none |
| `CUB_HEX`, serendipity | hexahedron | 8 + 12(p−1) | 1–3 | none |

The reference elements are those of `fortran_cubatures`. The
[user manual](manual/index.md) ([online](https://willklausler.github.io/fortran_shapefuncs/page/index.html),
with the [API documentation](https://willklausler.github.io/fortran_shapefuncs/)) describes the shape functions, the
[node numbering](manual/elements.md) and the
[infinite elements](manual/infinite-elements.md).

## Building, testing and documentation

```sh
fpm test --profile debug     # bounds checks, sanitizers, FP traps
fpm test --profile release
fpm run --example
make test                    # without fpm
ford ford.md                 # API documentation and manual in docs/
```

For every element and family, a sweep of orders and every combination of
infinite directions, the tests check nodal interpolation, exact reproduction of the
element's function space with first and second derivatives, finite
differences, and partition of unity. They also check the classical numbering
node by node, closed forms from the literature and caller-defined orderings. CI runs gfortran 13–15,
Intel ifx and NVIDIA nvfortran on Linux, gfortran on macOS and Windows, the
make build and the documentation build, and publishes the documentation
to GitHub Pages from `main`.

## Pedigree and support

The finite shape functions are the standard Lagrange and serendipity
elements [6, 7], with Silvester's formula for simplices [5]. The infinite elements are the mapped
infinite elements of Bettess [1] and Zienkiewicz, Emson and Bettess [2],
generalised to higher orders [3, 4]. Every function is verified by the tests
described above. The code has not been independently reviewed. It is
maintained by the author on a best-effort basis.

Please report bugs and request features through
[GitHub issues](https://github.com/willklausler/fortran_shapefuncs/issues).
Pull requests are welcome; please include a test.

## References

1. Bettess, P. (1977). Infinite elements. *International Journal for Numerical Methods in Engineering*, 11(1), 53–64. [doi:10.1002/nme.1620110107](https://doi.org/10.1002/nme.1620110107)
2. Zienkiewicz, O. C., Emson, C., Bettess, P. (1983). A novel boundary infinite element. *International Journal for Numerical Methods in Engineering*, 19(3), 393–404. [doi:10.1002/nme.1620190307](https://doi.org/10.1002/nme.1620190307)
3. Marques, J. M. M. C., Owen, D. R. J. (1984). Infinite elements in quasi-static materially nonlinear problems. *Computers & Structures*, 18(4), 739–751. [doi:10.1016/0045-7949(84)90019-1](https://doi.org/10.1016/0045-7949(84)90019-1)
4. Bettess, P. (1992). *Infinite Elements*. Penshaw Press.
5. Silvester, P. (1969). High-order polynomial triangular finite elements for potential problems. *International Journal of Engineering Science*, 7(8), 849–861. [doi:10.1016/0020-7225(69)90065-2](https://doi.org/10.1016/0020-7225(69)90065-2)
6. Zienkiewicz, O. C., Taylor, R. L., Zhu, J. Z. (2005). *The Finite Element Method: Its Basis and Fundamentals*, 6th ed. Butterworth-Heinemann.
7. Arnold, D. N., Awanou, G. (2011). The serendipity family of finite elements. *Foundations of Computational Mathematics*, 11(3), 337–344. [doi:10.1007/s10208-011-9087-3](https://doi.org/10.1007/s10208-011-9087-3)

## License

[MIT](LICENSE) © 2025-2026 Will Klausler
