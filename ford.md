---
project: fortran_shapefuncs
summary: Finite element shape functions and their derivatives at the points of a cubature, including serendipity and mapped infinite elements: the derived type `shapefunc`.
author: Will Klausler
github: https://github.com/willklausler
project_github: https://github.com/willklausler/fortran_shapefuncs
license: mit
src_dir: ./src
output_dir: ./build/doc
page_dir: ./build/pages
preprocess: false
display: public
source: true
search: true
graph: false
---

`use shapefuncs` gives the type `shapefunc`, the infinitude constants
`SHP_FIN`, `SHP_INF`, `SHP_CHP` and the families `SHP_LAGRANGE`,
`SHP_SERENDIPITY`. The cubature and the element constants come from
`use cubatures` (fortran_cubatures). Start with the
[user manual](page/index.html), which is also available as plain text in
[llms.txt](llms.txt).

```fortran
use cubatures
use shapefuncs
type(shapefunc) :: s
s = shapefunc(cubature(CUB_HEX, 4), 2)
print *, s%nnodes, sum(s%func(:,1))   ! 27  1.0000000000000000
```
