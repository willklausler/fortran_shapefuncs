# Changelog

All notable changes are listed here. The project follows [Semantic Versioning](https://semver.org).

## [0.2.0] - 2026-10-10

- Fixed: `call s%set(q, s%order)` built order 0, and `infin=s%infin` or `family=s%family` were likewise reset, because `set` cleared the object before it read its arguments.
- Works with fortran_cubatures 0.2.0, which no longer exports the kind `rk`. Take `real64` from `iso_fortran_env` instead.
- New layout following the maintainer's package standard: `llms.txt` API reference (replacing `manual/`), man page, `AGENTS.md`, `CITATION.cff`, `CONTRIBUTING.md`, six example programs and an app.
- fpm profiles with per-compiler flags; new make build that installs the library, module files and man page, without the pkg-config file.
- CI adds ifx on Windows and LLVM flang; documentation is deployed by a separate workflow.

## [0.1.0] - 2025-09-10

First public release.

- Type `shapefunc` with constructor and `set`: Lagrange shape functions with first and second derivatives at the points of a cubature, on lines, triangles, quadrilaterals, tetrahedra, hexahedra and wedges of any order.
- Mapped infinite elements: mapping (`SHP_INF`) and chopped (`SHP_CHP`) functions.
- Serendipity quadrilaterals and hexahedra of orders 1 to 3 (`SHP_SERENDIPITY`), and caller-defined node ordering (`nodes`).
- `eval`, `is_valid`, `summary`, `show`, `numbering` and `destroy`.
