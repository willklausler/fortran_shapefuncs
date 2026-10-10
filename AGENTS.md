# AGENTS.md

Instructions for coding agents working on this repository. The public API and
its conventions are in `llms.txt`; read it before changing behaviour.

## Layout

| Path | Role |
|---|---|
| `src/shapefuncs.f90` | **Source of truth**: the whole library, module `shapefuncs` |
| `test/test_shapefuncs.f90` | The tests: one program, `check(cond, name)` per assertion |
| `example/*.f90` | Demo programs; each is copied verbatim into `llms.txt` |
| `app/main.f90` | Interpolation error against order for every element type, the fpm app |
| `llms.txt` | API reference for agents, and the online user manual |
| `man/shapefuncs.3`, `README.md`, `ford.md` | Human documentation |
| `tools/check_docs.py` | CI check: every public name and example is in `llms.txt` |

The dependency `fortran_cubatures` (module `cubatures`) provides `cubature` and
the element constants. fpm fetches it; the Makefile clones it into
`build/make` unless `CUBATURES` points to its `src/cubatures.f90`.

## Commands

```sh
pip install ford                     # once
fpm update                           # fetch the latest fortran_cubatures
fpm test --profile debug             # fast check, gfortran
fpm test --profile release
fpm run --example --all
fpm run                              # the app
python3 tools/check_docs.py
./ci.sh                              # gfortran, ifx, nvfortran, flang if installed
make test FC=ifx                     # make build
make doc                             # FORD HTML in build/doc
```

Done means: the tests pass for every compiler you have in both profiles, the
examples and the app run, and `tools/check_docs.py` passes.

## Decided conventions (do not change without the maintainer)

- One type, `shapefunc`, with public components `elm`, `dim`, `order`, `family`, `nnodes`, `npoints`, `infin(3)`, `lattice`, `coords`, `func(nnodes, npoints)`, `derv(dim, nnodes, npoints)`, `curv(dim, dim, nnodes, npoints)`. Users read the components directly.
- Reference elements and the element constants `CUB_*` are those of fortran_cubatures; this package defines no element constants of its own.
- Nodes are equispaced and numbered hierarchically (vertices, edges, faces, interior), which reproduces the classical numbering for orders 1 to 3. The `nodes` argument permutes it once in `set`.
- Infinitude constants `SHP_FIN`, `SHP_INF`, `SHP_CHP` are 1 to 3; families `SHP_LAGRANGE`, `SHP_SERENDIPITY` are 1 and 2.
- Every shape function is a product of 1D functions of affine factor coordinates; serendipity functions are fixed combinations of the tensor-product functions. New element types should fit this design.
- Error policy: invalid arguments to the constructor or `set` stop with `error stop` and a message that starts with `shapefunc%set:`. The policy is also stated in `llms.txt` and the man page NOTES.
- Kinds are not exported. Internally `rk` aliases `real64`.
- No module state. Everything is pure except the output procedures `summary`, `show` and `numbering`.

## Adding or changing a procedure

1. Implement it in `src/shapefuncs.f90` with a `!!` line; helpers are private and `pure` (or `elemental`).
2. Bind it to the type (with a `!!` line) or add it to the `public` statement of its group.
3. Add checks to `test/test_shapefuncs.f90`; a new element or family goes through `check_element` over a sweep of orders.
4. Document it in `llms.txt` (the CI check requires the name), the man page and `CHANGELOG.md`.
5. Run the tests in both profiles, the examples and `tools/check_docs.py`.

## Test rules

- `check_element` proves each element correct: nodal interpolation (Kronecker delta), exact reproduction of the element's function space with first and second derivatives, fourth-order finite differences, symmetric second derivatives, and partition of unity (or vanishing at infinity for chopped directions). Classical numbering, closed forms, orderings, integration and the interface have their own checks.
- Tolerances are multiples of epsilon: `tol = 2**19*eps` for exact identities and `fdtol = 64*sqrt(eps)` for finite differences. Equispaced interpolation up to order 8 needs about `2**15*eps`, and the finite differences about `2**27*eps`; keep the margin when you add higher orders.
- Error stops cannot be tested.
- The tests must run under floating-point trapping and the sanitizers of the debug profile.

## Compiler quirks that shaped the code

None so far: the code needs no compiler-specific workaround.

## Commits

Commit as the maintainer, Will Klausler, with the GitHub no-reply address
`226649647+willklausler@users.noreply.github.com` (set in the repository's git config).
Do not add `Co-authored-by` trailers: the maintainer is the only contributor.
