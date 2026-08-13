# CLAUDE.md

Cross-repo rules live in `GitHub\CLAUDE.md` (4-core cap, C: drive policy,
git workflow, no Co-Authored-By, README/ALGORITHM templates, research-
library-first rule, etc.) — not repeated here, only what's specific to
KSROP. Fortran-specific gotchas (implicit-typing trap, no-adaptive-step,
gfortran loop-control bug) are in the shared, tree-walked
`.claude/rules/fortran-ks-gotchas.md` — also not repeated here.

## Project

Fortran KS-regularized orbit propagator (fpm package), the base library
that OREM, KSROP-Lunar, and KSROP-Mars all depend on via git+tag. The
general (n,m) tesseral geopotential, drag, SRP, and third-body force
models live here; the Moon/Mars-centered variants are separate repos that
consume this one, not forks of it.

## Build / Test

```bash
fpm build --compiler ifx        # or --compiler gfortran
fpm test  --compiler ifx
fpm run driver_KS --compiler ifx
```

Windows/Intel oneAPI needs both environments initialized in the same
shell session before `fpm`/`ifx`:
```bat
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"
```
CI (`.github/workflows/ci.yml`) matrix-tests both `ifx` and `gfortran`;
`release.yml` handles tagged releases.

## Key code

- `src/` — library (force models, KS transform, integrator, geopotential)
- `app/driver_KS.F` — main propagator entry point
- `test/` — fpm test suite
- `fpm.toml` — sets `[fortran] source-form = "fixed"` +
  `implicit-typing`/`implicit-external = true` (F77-style code; fpm's
  modern defaults otherwise reject it)

## Always / never

- Consumers pin KSROP by git tag (`{ git = "...", tag = "vX.Y.Z" }`) —
  never hand-copy KSROP source into a consumer repo. (The "fix
  KSROP-lineage bugs here first" rule is stated once, in
  `GitHub\CLAUDE.md` §5 — not repeated per-repo.)
