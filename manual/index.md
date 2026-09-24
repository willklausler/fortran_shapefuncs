---
title: User manual
ordered_subpage: getting-started.md
                 elements.md
                 infinite-elements.md
                 api.md
                 development.md
---

`fortran_shapefuncs` provides one derived type, `shapefunc`, which holds the
shape functions of a finite element and their first and second derivatives at
the points of a cubature (numerical integration rule) from
[fortran_cubatures](https://github.com/willklausler/fortran_cubatures).

It supports lines, triangles, quadrilaterals, tetrahedra, hexahedra and wedges
of any polynomial order. Lines, quadrilaterals, hexahedra and the axial
direction of wedges can also be *infinite*: they extend to infinity and use the
mapped infinite elements of Zienkiewicz, Emson and Bettess.

## Contents

1. [Getting started](getting-started.html): installation, a first program, and
   the element loop of a finite element code
2. [Elements and node numbering](elements.html): reference elements, orders,
   the shape functions used, and node numbering
3. [Infinite elements](infinite-elements.html): mapping and chopped functions,
   and how to build an infinite element
4. [API reference](api.html): constants, type components and procedures
5. [Development](development.html): building, testing, documentation, support
   and pedigree

The pages are Markdown files in the `manual/` directory of the repository.
`ford ford.md` renders them, together with the API documentation, into `docs/`.
