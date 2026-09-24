---
title: API reference
---

```fortran
use shapefuncs, only: shapefunc, SHP_FIN, SHP_INF, SHP_CHP
use cubatures,  only: cubature, rk, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
```

The element constants `CUB_*`, the real kind `rk` (`real64`) and the type
`cubature` come from `fortran_cubatures`.

## Constants

| Constant | Meaning |
| --- | --- |
| `SHP_FIN` | Finite direction: Lagrange polynomials on \([-1,1]\) |
| `SHP_INF` | Infinite direction: mapping functions, for coordinates |
| `SHP_CHP` | Infinite direction: chopped functions, for the solution |

## Type `shapefunc`

### Components

| Component | Type | Contents |
| --- | --- | --- |
| `elm` | `integer` | Element type, one of `CUB_*`; 0 when unset |
| `dim` | `integer` | Spatial dimension |
| `order` | `integer` | Polynomial order \(p\) |
| `nnodes` | `integer` | Number of nodes |
| `npoints` | `integer` | Number of cubature points |
| `infin(3)` | `integer` | Infinitude per direction, `SHP_*`; unused entries are `SHP_FIN` |
| `lattice(:, nnodes)` | `integer`, allocatable | Lattice indices of each node, see [Elements](elements.html) |
| `coords(dim, nnodes)` | `real(rk)`, allocatable | Nodal coordinates in the reference element |
| `func(nnodes, npoints)` | `real(rk)`, allocatable | \(N_i(\xi_g)\) |
| `derv(dim, nnodes, npoints)` | `real(rk)`, allocatable | \(\partial N_i/\partial\xi_j\) at \(\xi_g\) |
| `curv(dim, dim, nnodes, npoints)` | `real(rk)`, allocatable | \(\partial^2 N_i/\partial\xi_j\partial\xi_k\) at \(\xi_g\), symmetric in \(j, k\) |

Treat the components as read-only: `set` fills them consistently.

### Constructor

```fortran
s = shapefunc(q, order [, infin])
```

Same arguments as `set`.

### `call s%set(q, order [, infin])`

Build the shape functions in place, discarding any previous contents.

| Argument | Type | Meaning |
| --- | --- | --- |
| `q` | `type(cubature)`, in | Cubature whose points are used; its element type sets the element |
| `order` | `integer`, in | Polynomial order |
| `infin(:)` | `integer`, in, optional | `SHP_*` per direction, size 1 (all directions) or `q%dim`; default `SHP_FIN` |

`set` stops with `error stop` and a message when

- `q` is not set;
- `order` is negative, or below 1 for triangles, tetrahedra, wedges or any
  infinite direction;
- `infin` has the wrong size or a value other than `SHP_*`;
- a triangle or tetrahedron has an infinite direction, or a wedge has one
  other than direction 3.

### `call s%eval(xi, func, derv, curv)`

Evaluate at any point `xi(dim)` of the reference element. The outputs have
the shapes `func(nnodes)`, `derv(dim, nnodes)` and `curv(dim, dim, nnodes)`.
The object must be set. Infinite directions are singular at \(\xi = 1\).
`set` calls `eval` at every cubature point.

### `s%is_valid()`

`.true.` if the object is set and its arrays have consistent shapes.

### `call s%summary([unit])`

Write the element type, dimension, order, infinitude and the numbers of nodes
and points. `unit` defaults to `output_unit`.

### `call s%show([unit])`

Write the summary, then for every point and node the value and first
derivatives, and for every point the sums over the nodes (1 and 0 unless a
direction is chopped).

### `call s%numbering([unit])`

Write the summary, then the number and reference coordinates of every node.

### `call s%destroy()`

Deallocate and reset to the unset state.

## Purity and threads

The constructor, `set`, `eval`, `is_valid` and `destroy` are `pure`. Once
set, a `shapefunc` can be shared read-only between threads.
