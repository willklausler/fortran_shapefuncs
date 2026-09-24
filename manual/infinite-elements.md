---
title: Infinite elements
---

Infinite elements model unbounded domains, such as soil under a foundation,
the far field of a potential flow, or a body in free space. One layer of
infinite elements surrounds the finite mesh and extends it to infinity.
`fortran_shapefuncs` implements the *mapped* infinite elements of Bettess [1]
in the form of Zienkiewicz, Emson and Bettess [2], generalised to any order as
by Marques and Owen [3]. Bettess [4] gives a full treatment.

## Directions and nodes

Any direction of a line, quadrilateral or hexahedron, and the axial direction
\(\xi_3\) of a wedge, can be infinite. The element extends to infinity at
\(\xi = 1\) in that direction. An infinite direction of order \(p\) has
\(p + 1\) nodes at

\[ t_a = -1 + \frac{2a}{p+1}, \qquad a = 0, \dots, p, \]

which is the grid of order \(p + 1\) without the node \(t_{p+1} = 1\) at
infinity. Linear (\(p=1\)): \(t = -1, 0\). Quadratic: \(t = -1, -1/3, 1/3\).

Choose the infinitude of each direction with the `infin` argument:

| Constant | Functions | Use for |
| --- | --- | --- |
| `SHP_FIN` | Lagrange polynomials on \([-1,1]\) (default) | finite directions |
| `SHP_INF` | mapping functions \(M_a\) | the coordinates |
| `SHP_CHP` | chopped functions \(\ell^{+}_a\) | the solution |

An infinite element therefore uses **two** `shapefunc` objects at the same
cubature points: one with `SHP_INF` for the geometry and one with `SHP_CHP`
for the field. Both have the same nodes.

```fortran
q     = cubature(CUB_QUA, 6)
geom  = shapefunc(q, 2, [SHP_FIN, SHP_INF])   ! infinite in xi_2
field = shapefunc(q, 2, [SHP_FIN, SHP_CHP])
```

## Mapping functions

With \(\ell_a\) the Lagrange polynomials on the nodes \(t_0, \dots, t_p\),

\[ M_a(\xi) = \ell_a(\xi)\,\frac{1 - t_a}{1 - \xi}. \]

They are 1 at their own node and 0 at the others, they sum to 1 (because
\(\ell_a\) interpolate the linear function \(1-\xi\) exactly), and they are
singular at \(\xi = 1\). The linear functions are those of [2]:

\[ M_0 = \frac{-2\xi}{1-\xi}, \qquad M_1 = \frac{1+\xi}{1-\xi}. \]

The map \(x(\xi) = \sum_a M_a(\xi)\,x_a\) sends \(\xi \to 1\) to infinity. In
one dimension with nodes \(x_0\) and \(x_1\), it is

\[ x(\xi) = x_P + \frac{2\,(x_1 - x_0)}{1 - \xi}, \qquad x_P = 2x_0 - x_1, \]

so the element behaves like a ray from the *pole* \(x_P\): place the nodes
along rays from a point inside the finite domain, with \(x_1 - x_0 = x_0 - x_P\).
Then \(1 - \xi = 2(x_1 - x_0)/(x - x_P)\), which is proportional to \(1/r\)
with \(r\) the distance from the pole.

The mapping functions of order \(p\) reproduce every function
\(q(\xi)/(1-\xi)\) with \(q\) of degree \(\le p\), including
\(x_P + c/(1-\xi)\). So for any order the map is the same ray from the pole
when node \(a\) sits at distance \(c/(1-t_a)\) from the pole. For
\(p = 2\), the three nodes are at distances in the ratio \(2 : 3 : 6\).

## Chopped functions

The chopped functions are the Lagrange polynomials of order \(p+1\) on the
nodes \(t_0, \dots, t_{p+1} = 1\), without the function of the node at
infinity:

\[ \ell^{+}_a(\xi) = \prod_{j \ne a,\ j \le p+1} \frac{\xi - t_j}{t_a - t_j},
   \qquad a = 0, \dots, p. \]

Interpolating with them sets the field to zero at infinity. Because
\(1 - \xi \propto 1/r\), a polynomial in \(\xi\) that vanishes at \(\xi = 1\)
is a polynomial in \(1/r\) without a constant term, so the field decays as

\[ u(r) = \frac{\alpha_1}{r} + \frac{\alpha_2}{r^2} + \dots + \frac{\alpha_{p+1}}{r^{p+1}}, \]

the expansion of many far fields [2, 4]. The chopped functions do **not** sum
to 1; the tests check that they vanish at \(\xi = 1\) instead. For a field that
tends to a nonzero constant, subtract that constant first.

## Integration

The integrands of infinite elements are rational, not polynomial, so no rule
integrates them exactly in general. Gauss–Legendre rules work well because
their points stay away from \(\xi = 1\). Raise the degree until the result
converges. Integrals such as \(\int u^2\,dx\) with \(u \sim 1/r\) are finite,
and the combination of the \(1/(1-\xi)^2\) Jacobian with \((1-\xi)^2\) from
\(u^2\) often gives a polynomial integrand.

The example program integrates \(u = 1/x\) over \([1, \infty)\) with one
linear infinite element. The nodes are \(x = 1, 2\), so the pole is at 0 and
\(x = 2/(1-\xi)\). The chopped functions represent \(1/x\) exactly, the
integrand \(u^2\,dx/d\xi = 1/2\) is constant, and the result 1 is exact.

## Wedges

Wedges can be infinite only along \(\xi_3\), for example to extend a mesh of
triangular prisms downwards to infinity:

```fortran
q = cubature(CUB_WED, [4, 6])
s = shapefunc(q, 2, [SHP_FIN, SHP_FIN, SHP_INF])
```

Triangles and tetrahedra cannot be infinite.

## References

1. Bettess, P. (1977). Infinite elements. *International Journal for Numerical
   Methods in Engineering*, 11(1), 53–64. [doi:10.1002/nme.1620110107](https://doi.org/10.1002/nme.1620110107)
2. Zienkiewicz, O. C., Emson, C., Bettess, P. (1983). A novel boundary
   infinite element. *International Journal for Numerical Methods in
   Engineering*, 19(3), 393–404. [doi:10.1002/nme.1620190307](https://doi.org/10.1002/nme.1620190307)
3. Marques, J. M. M. C., Owen, D. R. J. (1984). Infinite elements in
   quasi-static materially nonlinear problems. *Computers & Structures*,
   18(4), 739–751. [doi:10.1016/0045-7949(84)90019-1](https://doi.org/10.1016/0045-7949(84)90019-1)
4. Bettess, P. (1992). *Infinite Elements*. Penshaw Press.
