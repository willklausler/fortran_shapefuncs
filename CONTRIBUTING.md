# Contributing

Bug reports, questions and pull requests are welcome.

- **Bugs and questions:** open an [issue](https://github.com/willklausler/fortran_shapefuncs/issues) with the compiler, its version, the operating system, and a minimal program.
- **Pull requests:** fork, branch, and open a PR against `main`. CI must pass.

Before opening a PR:

1. Add tests to `test/test_shapefuncs.f90` for every change.
2. Run `fpm test --profile debug` and `fpm test --profile release`.
3. Document new or changed public procedures in `llms.txt`, `man/shapefuncs.3` and `CHANGELOG.md`; `python3 tools/check_docs.py` must pass.

[AGENTS.md](AGENTS.md) describes the layout, conventions and compiler quirks in detail; it applies to human contributors too.
