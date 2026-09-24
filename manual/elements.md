---
title: Elements and node numbering
---

## Reference elements

The reference elements and coordinates \(\xi\) are those of `fortran_cubatures`.

| Constant | Element | Reference domain | Nodes, order \(p\) | Orders |
| --- | --- | --- | --- | --- |
| `CUB_LIN` | line | \([-1,1]\) | \(p+1\) | \(\ge 0\) |
| `CUB_QUA` | quadrilateral | \([-1,1]^2\) | \((p+1)^2\) | \(\ge 0\) |
| `CUB_HEX` | hexahedron | \([-1,1]^3\) | \((p+1)^3\) | \(\ge 0\) |
| `CUB_TRI` | triangle | \(\xi_1,\xi_2 \ge 0,\ \xi_1+\xi_2 \le 1\) | \((p+1)(p+2)/2\) | \(\ge 1\) |
| `CUB_TET` | tetrahedron | \(\xi_i \ge 0,\ \xi_1+\xi_2+\xi_3 \le 1\) | \((p+1)(p+2)(p+3)/6\) | \(\ge 1\) |
| `CUB_WED` | wedge / prism | triangle \(\times\ [-1,1]\) | \((p+1)^2(p+2)/2\) | \(\ge 1\) |
| `CUB_QUA`, serendipity | quadrilateral | \([-1,1]^2\) | \(4p\) | 1–3 |
| `CUB_HEX`, serendipity | hexahedron | \([-1,1]^3\) | \(8 + 12(p-1)\) | 1–3 |

Common cases:

| Order | Line | Quad | Hex | Triangle | Tet | Wedge | Serendipity quad | Serendipity hex |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2 | 4 | 8 | 3 | 4 | 6 | 4 | 8 |
| 2 | 3 | 9 | 27 | 6 | 10 | 18 | 8 | 20 |
| 3 | 4 | 16 | 64 | 10 | 20 | 40 | 12 | 32 |
| 4 | 5 | 25 | 125 | 15 | 35 | 75 | | |

There is no upper limit on the order. The nodes are equispaced, so
interpolation becomes ill-conditioned at high orders (Runge's phenomenon).
The unit tests cover orders up to 8 for lines and triangles, 6 for
quadrilaterals and tetrahedra, 5 for wedges and 4 for hexahedra.

## Shape functions

- **Lines, quadrilaterals, hexahedra**: tensor products of Lagrange polynomials
  of degree \(p\) on the equispaced nodes \(-1 + 2a/p\). These are the full
  Lagrange elements [2]. Order 0 is one constant function with
  its node at the centre.
- **Triangles, tetrahedra**: the complete polynomials of total degree \(p\).
  With barycentric coordinates \(\lambda_k\) (for a triangle
  \(\lambda = (\xi_1, \xi_2, 1-\xi_1-\xi_2)\)), node \(i\) at
  \(\lambda = a_i/p\) has
  \[ N_i = \prod_k P_{a_{ki}}(\lambda_k), \qquad
     P_a(\lambda) = \prod_{j=0}^{a-1} \frac{p\lambda - j}{j+1}, \]
  Silvester's formula [1].
- **Wedges**: the triangle functions in \((\xi_1, \xi_2)\) times the Lagrange
  polynomials in \(\xi_3\).
- **Serendipity quadrilaterals and hexahedra** (`family=SHP_SERENDIPITY`,
  orders 1 to 3): nodes on the vertices and edges only. They span the
  monomials of *superlinear degree* \(\le p\), the degree counting only the
  variables of exponent 2 or more [3]. For \(p=2\) this is the classical
  8-node quadrilateral, with corner and midside functions
  \[ N_i = \tfrac14(1+\xi\xi_i)(1+\eta\eta_i)(\xi\xi_i+\eta\eta_i-1), \qquad
     N_i = \tfrac12(1-\xi^2)(1+\eta\eta_i), \]
  and the 20-node hexahedron. `set` computes each serendipity function once
  as a combination of the tensor-product Lagrange functions of the same
  order, so the derivatives are exact. Serendipity elements must be finite.

The derivatives are exact. They follow from the product rule and the constant
derivatives of the barycentric coordinates.

## Node numbering

Nodes are numbered hierarchically: vertices, then edges, then faces, then the
interior. Orders 1 to 3 give the classical numbering used by most finite
element codes. Higher orders follow the same rules. Serendipity elements
keep the vertex and edge nodes of the Lagrange element in the same order:
the 8-node quadrilateral is the 9-node one without its centre, and the
20-node hexahedron is the 27-node one without its face and centre nodes.

### Vertices

- Quadrilateral: counter-clockwise from \((-1,-1)\).
- Hexahedron: the bottom face \(\xi_3=-1\) counter-clockwise from
  \((-1,-1,-1)\), then the top face \(\xi_3=1\) in the same order.
- Triangle and tetrahedron: vertex \(k\) is where barycentric coordinate \(k\)
  is 1, so vertex 1 is \((1,0)\) or \((1,0,0)\) and the last vertex is the
  origin.
- Wedge: the triangle vertices at \(\xi_3=-1\), then at \(\xi_3=1\).

### Edges

Edge nodes run from the first vertex of the edge to the second.

| Element | Edges |
| --- | --- |
| quadrilateral | 1-2, 2-3, 3-4, 4-1 |
| hexahedron | 1-2, 2-3, 3-4, 4-1, 5-6, 6-7, 7-8, 8-5, 1-5, 2-6, 3-7, 4-8 |
| triangle | 1-2, 2-3, 3-1 |
| tetrahedron | 1-2, 2-3, 3-1, 1-4, 2-4, 3-4 |
| wedge | 1-2, 2-3, 3-1, 4-5, 5-6, 6-4, 1-4, 2-5, 3-6 |

### Faces and interiors

- Hexahedron faces: \(\xi_3=-1\), \(\xi_3=1\), \(\xi_1=-1\), \(\xi_1=1\),
  \(\xi_2=-1\), \(\xi_2=1\). A face in \((\xi_j,\xi_k)\), \(j<k\), is numbered
  as a quadrilateral in those two coordinates.
- Tetrahedron faces: 1-2-3, 1-2-4, 2-3-4, 1-3-4, each numbered as a triangle
  with its vertices in that order.
- Wedge faces: the bottom triangle, the top triangle, then the quadrilaterals
  on edges 1-2, 2-3 and 3-1. A quadrilateral face is numbered in (position
  along the edge, \(\xi_3\)).
- The interior of a quadrilateral, hexahedron, triangle or tetrahedron (and of
  each face) is numbered recursively as the same element of order \(p-2\)
  (quadrilateral, hexahedron), \(p-3\) (triangle) or \(p-4\) (tetrahedron),
  shrunk onto the interior nodes.
- Wedge interior nodes are the triangle interiors, layer by layer from
  \(\xi_3=-1\).

### Examples

The 16-node quadrilateral (\(p=3\)), with \(\xi_1\) to the right and
\(\xi_2\) up:

```text
 4 --10 -- 9 -- 3
 |              |
11   16   15    8
 |              |
12   13   14    7
 |              |
 1 -- 5 -- 6 -- 2
```

The 15-node triangle (\(p=4\)), with vertex 3 at the origin:

```text
 2
 | \
 7   6
 |     \
 8  14   5
 |         \
 9  15  13   4
 |             \
 3 --10--11--12-- 1
```

In the triangle, edge 1-2 holds nodes 4, 5, 6 (from vertex 1), edge 2-3 holds
7, 8, 9, edge 3-1 holds 10, 11, 12, and the interior 13, 14, 15 is a linear
triangle numbered from the corner nearest vertex 1. `call s%numbering()`
prints the coordinates of every node of any element.

### Caller-defined node ordering

Mesh generators and file formats number the nodes of higher-order elements
differently. Pass the `nodes` argument to get every output in your
numbering:

```fortran
s = shapefunc(q, 2, family=SHP_SERENDIPITY, nodes=perm)
```

Your node `k` is node `perm(k)` of the default numbering. The easiest way to
build `perm` is to match reference coordinates. Write down the reference
coordinates `xmine(:,k)` of your node `k`, then find the node of a default
`shapefunc` at the same place:

```fortran
s = shapefunc(q, 2, family=SHP_SERENDIPITY)
do k = 1, s%nnodes
  perm(k) = minloc(sum(abs(s%coords - spread(xmine(:,k), 2, s%nnodes)), dim=1), 1)
end do
s = shapefunc(q, 2, family=SHP_SERENDIPITY, nodes=perm)
```

Section 5 of the example program does this for a 20-node hexahedron whose
edge nodes are listed as in Gmsh.

### Lattice indices

`s%lattice(:, i)` holds the integer position of node \(i\):

- lines, quadrilaterals, hexahedra: one index \(0 \dots p\) per direction;
- triangles, tetrahedra: the barycentric indices \(a_k\), which sum to \(p\);
- wedges: three barycentric indices, then the axial index \(0 \dots p\).

They locate each node exactly, without comparing floating-point
coordinates. For example, a node lies on the boundary of a triangle where one of
its barycentric indices is 0.

## References

1. Silvester, P. (1969). High-order polynomial triangular finite elements for
   potential problems. *International Journal of Engineering Science*, 7(8),
   849–861. [doi:10.1016/0020-7225(69)90065-2](https://doi.org/10.1016/0020-7225(69)90065-2)
2. Zienkiewicz, O. C., Taylor, R. L., Zhu, J. Z. (2005). *The Finite Element
   Method: Its Basis and Fundamentals*, 6th ed. Butterworth-Heinemann.
3. Arnold, D. N., Awanou, G. (2011). The serendipity family of finite
   elements. *Foundations of Computational Mathematics*, 11(3), 337–344.
   [doi:10.1007/s10208-011-9087-3](https://doi.org/10.1007/s10208-011-9087-3)
