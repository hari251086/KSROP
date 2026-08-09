# KSROP — KS Regular Orbit Propagator

Orbit propagation using **Kustaanheimo–Stiefel (KS) regular elements** with a Runge–Kutta–Gill 4th-order numerical integrator for Earth-orbiting satellites.

[![CI](https://github.com/hari251086/KSROP/actions/workflows/ci.yml/badge.svg)](https://github.com/hari251086/KSROP/actions/workflows/ci.yml)

**Author:** Harishkumar Sellamuthu · hari251086@gmail.com  
**Copyright:** 2018, Harishkumar Sellamuthu, All Rights Reserved

---

## 1. Features

| Perturbation | Model | Notes |
|---|---|---|
| Earth oblateness (Jn) | EGM2008 | Configurable up to degree 2190 |
| Luni-solar gravity | Analytical ephemeris | Configurable Legendre degree |
| Atmospheric drag | Oblate co-rotating exponential | Tabulated density (ATM.DAT), perigee-referenced |
| Solar radiation pressure | Cannon-ball | Cylindrical or conical shadow models |

---

## 2. Project Structure

```
KSROP/
├── fpm.toml                              fpm package manifest (src/app/test layout below)
├── src/
│   ├── Subrouts.F                        Shared subroutines (transforms, I/O, force models)
│   ├── Legendre.F                        Zonal Legendre polynomial evaluation
│   ├── TLEread.F                         TLE reader + SGP4 + TEME→J2000 conversion
│   ├── swx.F                             Epoch-resolved space weather (F10.7/Kp, issue #26)
│   ├── Cunningham.F                      General (n,m) tesseral harmonics, Cartesian solid-harmonic
│   │                                     recursion (issue #30) -- standing cross-check, not called
│   │                                     from driver_KS.F (see LegendreTess.F)
│   ├── LegendreTess.F                    General (n,m) tesseral harmonics, classical Legendre
│   │                                     expansion (issue #30) -- THE active geopotential
│   └── jr71_profile.F                    Jacchia-71 atmosphere profile (shared by gen_atm* tools)
├── app/
│   ├── driver_KS.F                       Main propagator program
│   ├── tle2opm.F                         TLE-to-OPM converter tool
│   ├── gen_atm_jr71.F                    1-D Jacchia-71 table generator
│   └── gen_atm2d_jr71.F                  2-D T-inf-gridded generator (issue #26)
│
├── input/
│   ├── const_new.dat                    Physical constants (single source of truth)
│   ├── input.opm                        Initial state for driver_KS (CCSDS OPM v2.0)
│   ├── input.dat                        Simulation parameters
│   ├── tle2opm.cfg                      TLE-to-OPM configuration
│   ├── KSROP_20250501T200802_47944.opm  Generated OPM (NORAD 47944, 2025-05-01)
│   ├── ATM.DAT                          Atmosphere density table (60–630 km)
│   ├── EGM2008_to2190_TideFree          Geopotential coefficients (~231 MB)
│   ├── example_47944.tle.txt            TLE data: 3,977 entries, NORAD 47944
│   ├── example_multi.tle.txt            TLE data: 94,597 entries, multi-satellite
│   └── cdm_sample.kvn                   CCSDS CDM sample (508.0-B-1 Annex)
│
├── output/                              Runtime-generated files
│   ├── KSROP_YYYYMMDDTHHMMSS.oem       Trajectory (CCSDS OEM v2.0)
│   ├── KSROP_YYYYMMDDTHHMMSS_Regular.out  KS elements debug dump
│   └── ksrop.opm                        Initial Keplerian elements (OPM)
│
├── test/
│   ├── test_subrouts.F                   Unit tests (82 checks)
│   ├── test_tle.F                        TLE reader tests (147 checks)
│   ├── test_tle2sv.F                     SGP4/SDP4/frame tests (156 checks)
│   ├── test_tle2opm.F                    TLE-to-OPM pipeline tests (21 checks)
│   ├── test_bugs.F                       Bug regression tests (17 checks)
│   ├── test_sw.F                         Space weather tests (17 checks)
│   └── test_cunningham.F                 General tesseral harmonics tests (25 checks, issue #30)
├── test_driver.py                       Integration test (10 checks)
├── test_initial_conditions.py           Multi-orbit test (110 checks)
├── benchmark.py                         Performance profiler
├── lint_check.sh                        Automated F77 lint checks (line length, precision, common blocks)
├── test_all.sh                          Unified test runner (lint + all test suites)
├── Makefile                             Unix build
└── build.bat                            Windows build (Intel Fortran)
```

---

## 3. Quick Start (User Manual)

This section is a guided walkthrough for a first-time user. For full detail on any
step, follow the link to its reference section below.

### Prerequisites

- A Fortran compiler: **GNU Fortran** (`gfortran`, free, any OS) or **Intel oneAPI
  Fortran** (`ifx`, Windows/Linux). See [Building](#4-building).
- **Python 3** (standard library only, no packages to install) -- only needed to run
  the integration tests, not to build or run the propagator itself.
- Optional: the **EGM2008** geopotential coefficient file if you want oblateness
  degrees above point-mass -- see the note under [Step 2](#step-2-run-a-propagation)
  below. It is not bundled with this repository (~231 MB).

### Step 1: Build

```bash
make            # Unix/Linux/macOS, gfortran -- see section 4 for Windows/ifx
```

This produces `driver_KS` (the propagator) in the repository root. `make tools`
additionally builds `tle2opm` (the TLE-to-OPM converter).

### Step 2: Run a propagation

```bash
./driver_KS
```

The tracked `input/input.opm` (initial state), `input/input.dat` (simulation
parameters), and `input/const_new.dat` (physical constants) are read as-is, so this
runs immediately with no editing required -- **except** that the tracked
`input/const_new.dat` requests `ngeo_deg = 50` (degree-50 EGM2008 gravity), and the
EGM2008 coefficient file itself is not part of this repository. Before your first
run, either:

- **Download EGM2008** and place it at `input/EGM2008_to2190_TideFree` (see
  [Known Issues](#12-known-issues) for the expected format), or
- **Edit `input/const_new.dat`** and set the second line's first value (`ngeo_deg`)
  to `0` for point-mass gravity (fastest way to confirm the build works).

A successful run prints progress to the console and writes two files to `output/`:
a CCSDS OEM trajectory (`KSROP_<timestamp>.oem`) and a KS-elements debug dump
(`KSROP_<timestamp>_Regular.out`) -- see [Output Files](#7-output-files).

### Step 3: Verify with the test suite

```bash
make test       # build everything + run test_all.sh (585 checks + 5 lint rules)
```

If this passes, your build is behaving identically to the one verified in CI (see
the badge at the top of this file).

### Step 4: Common customizations

All simulation parameters live in `input/input.dat` (4 lines: revolutions/step count,
force-model on/off flags, drag parameters, SRP parameters) -- see
[Input Files](#6-input-files) for the full field-by-field reference. A few common
edits:

| I want to... | Edit |
|---|---|
| Change how many orbits are propagated | `input.dat` line 1, first value (`nrev`) |
| Turn atmospheric drag off | `input.dat` line 3, second value (`IDRAG`) to `0` |
| Turn solar radiation pressure off | `input.dat` line 4, third value (`IPSR`) to `0` |
| Change the ballistic coefficient | `input.dat` line 3, first value (`BN`, kg/m²) |
| Run a pure two-body case (no perturbations) | `const_new.dat` line 2 all zeros, `input.dat` line 3/4 `IDRAG=0 IPSR=0` |
| Start from a different orbit | Edit the `X/Y/Z/X_DOT/Y_DOT/Z_DOT` state vector in `input/input.opm` |
| Start from a real satellite's TLE | See below |

### Step 5: Propagate from a real TLE

```bash
./tle2opm     # reads input/tle2opm.cfg (TLE file, NORAD ID, target epoch)
./driver_KS   # tle2opm already wrote a fresh input/input.opm for you
```

See [Running](#5-running) for the full TLE-to-OPM pipeline description.

### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `forrtl`/runtime error opening `EGM2008_to2190_TideFree` | `ngeo_deg > 1` but the file isn't present | Set `ngeo_deg = 0` in `const_new.dat`, or supply the file (see Step 2) |
| `list-directed I/O syntax error` on startup | A stray BOM or non-ASCII byte at the start of an input file | Re-save the file as plain ASCII/UTF-8 without a byte-order mark |
| Build errors on Linux mentioning `Rank mismatch` / `Expecting a scalar INTEGER` | You're building an older checkout predating the gfortran portability fixes (v2.0.0+) | Pull the latest `master`/tag |
| Nothing happens / no `output/` files | You're running from the wrong directory | All paths are relative to the repository root -- run `driver_KS`/`tle2opm` from there |

---

## 4. Building

Requires **Intel oneAPI Fortran** (`ifx`) or **GNU Fortran** (`gfortran`), plus a C/C++ linker (MSVC on Windows).

### fpm (Fortran Package Manager) — recommended, and the only way to consume KSROP as a dependency

```bash
fpm build --compiler ifx        # or --compiler gfortran
fpm test  --compiler ifx
fpm run driver_KS --compiler ifx
```

To depend on KSROP from another fpm project (e.g. OREM), add to that project's `fpm.toml`:

```toml
[dependencies]
ksrop = { git = "https://github.com/hari251086/KSROP", tag = "v2.1.0" }
```

`fpm.toml` sets `[fortran] source-form = "fixed"` (fpm defaults `.F` to free-form, which breaks
this codebase's column-1 comments) and `implicit-typing = true` / `implicit-external = true`
(this is 50-year-old-style F77 that relies on implicit typing and external-function declarations
in several places — fpm's stricter modern defaults otherwise reject valid code).

### Windows — Intel oneAPI ifx 2025.0 (manual, without fpm)

The Intel compiler requires both the Visual Studio and Intel environments to be initialised in the same shell session:

```bat
:: Single-command build (from cmd.exe or PowerShell via cmd /c)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"
cd /d "C:\Users\hari2\OneDrive\Documents\GitHub\KSROP"

:: Propagator
ifx app\driver_KS.F src\Subrouts.F src\Legendre.F src\TLEread.F /exe:driver_KS.exe

:: TLE-to-OPM converter
ifx app\tle2opm.F src\Subrouts.F src\TLEread.F src\Legendre.F /exe:tle2opm.exe

:: Tests
ifx test\test_subrouts.F src\Subrouts.F src\Legendre.F /exe:test_subrouts.exe
ifx test\test_bugs.F src\Subrouts.F src\Legendre.F /exe:test_bugs.exe
ifx test\test_tle.F src\Subrouts.F src\TLEread.F src\Legendre.F /exe:test_tle.exe
ifx test\test_tle2sv.F src\Subrouts.F src\TLEread.F src\Legendre.F /exe:test_tle2sv.exe
ifx test\test_tle2opm.F src\Subrouts.F src\TLEread.F src\Legendre.F /exe:test_tle2opm.exe
```

> **Note:** Use `/exe:name.exe` (not `-o`) for the output flag with `ifx` on Windows.

### Unix / Linux / macOS

```bash
make            # build driver_KS (FC=gfortran by default)
make tools      # build tle2opm
make tests      # build all test executables
make test       # build everything + run test_all.sh (lint + all tests)
make clean      # remove build artefacts
```

CI builds and runs the full suite on every push using this exact path (`.github/workflows/ci.yml`) -- see the badge at the top of this file.

### Manual (gfortran)

```bash
gfortran app/driver_KS.F src/Subrouts.F src/Legendre.F src/TLEread.F -o driver_KS.exe
gfortran app/tle2opm.F src/Subrouts.F src/TLEread.F src/Legendre.F -o tle2opm.exe
gfortran test/test_subrouts.F src/Subrouts.F src/Legendre.F -o test_subrouts.exe
gfortran test/test_bugs.F src/Subrouts.F src/Legendre.F -o test_bugs.exe
gfortran test/test_tle.F src/Subrouts.F src/TLEread.F src/Legendre.F -o test_tle.exe
gfortran test/test_tle2sv.F src/Subrouts.F src/TLEread.F src/Legendre.F -o test_tle2sv.exe
gfortran test/test_tle2opm.F src/Subrouts.F src/TLEread.F src/Legendre.F -o test_tle2opm.exe
```

---

## 5. Running

### Propagation

```bash
./driver_KS.exe
```

Reads `input/input.opm` (initial state), `input/input.dat` (parameters), `input/const_new.dat` (constants). Writes trajectory to `output/KSROP_*.oem`.

### TLE-to-OPM conversion

```bash
./tle2opm.exe
```

Reads `input/tle2opm.cfg`, selects closest TLE entry by NORAD/epoch, converts via SGP4 to J2000, writes `input/KSROP_<EPOCH>_<NORAD>.opm` and `input/input.opm`.

### Tests

```bash
./test_subrouts.exe                          # 67 unit tests
./test_bugs.exe                              # 17 regression tests
./test_tle.exe                               # 147 TLE reader tests
./test_tle2sv.exe                            # 156 SGP4/SDP4/frame tests
./test_tle2opm.exe                           # 21 pipeline tests
python test_driver.py driver_KS.exe          # 10 integration tests
python test_initial_conditions.py driver_KS.exe  # 110 multi-orbit tests

# Lint + all tests in one command
bash test_all.sh                             # lint + 419 Fortran + 120 integration
bash lint_check.sh                           # lint only (line length, precision, common blocks)
```

**Total: 585 automated checks + 5 lint rules.**

---

## 6. Input Files

### `const_new.dat` — Physical constants

Single source of truth for all programs. Read by `driver_KS` directly and by all other programs via `init_constants()`.

```
mu  R_Earth  AU  mu_Sun  mu_Moon       ! Line 1: gravitational/geometric constants
ngeo_deg  nsun_deg  nmoon_deg          ! Line 2: force model degrees
PSR_srp                                ! Line 3: solar radiation pressure (N/m²)
```

| Parameter | Value | Units |
|---|---|---|
| `mu` | 3.986004415×10⁵ | km³/s² |
| `R_Earth` | 6378.1363 | km |
| `AU` | 1.495978707×10⁸ | km |
| `mu_Sun` | 1.32712440018×10¹¹ | km³/s² |
| `mu_Moon` | 4.902801076×10³ | km³/s² |
| `ngeo_deg` | 0–2190 | Geopotential degree (0 or 1 = point mass) |
| `nsun_deg` | 0–2190 | Solar gravity Legendre degree |
| `nmoon_deg` | 0–2190 | Lunar gravity Legendre degree |
| `PSR_srp` | 4.56×10⁻⁶ | Solar radiation pressure at 1 AU (N/m²) |

### `input.opm` — Initial state (CCSDS OPM v2.0)

```
CCSDS_OPM_VERS = 2.0
CREATION_DATE  = 2025-05-01T20:08:02.901
ORIGINATOR     = KSROP

META_START
OBJECT_NAME    = SATELLITE
CENTER_NAME    = EARTH
REF_FRAME      = EME2000
TIME_SYSTEM    = UTC
META_STOP

STATE_VECTOR
EPOCH          = 2025-05-01T20:08:02.901
X              =      5652.381088 [km]
Y              =      3686.474468 [km]
Z              =       -14.037824 [km]
X_DOT          =      0.548211943 [km/s]
Y_DOT          =     -0.840302841 [km/s]
Z_DOT          =      7.623865946 [km/s]
```

> Output `ksrop.opm` uses the same format — can be fed back as `input.opm` for chained runs.

### `input.dat` — Simulation parameters

```
nrev  istep  tole               ! Revolutions, steps/revolution, integrator tolerance
n_geo  n_sun  n_moon            ! Force model flags (non-zero = on, 0 = off)
BN  IDRAG  WE_rot  EPS_f  FR_rot   ! Drag parameters
CR  AM  IPSR  ISHAD             ! SRP parameters
```

| Parameter | Example | Description |
|---|---|---|
| `BN` | `50.0` | Ballistic coefficient CdA/m (kg/m²) |
| `IDRAG` | `1` | Drag on/off |
| `WE_rot` | `7.2921150d-5` | Earth rotation rate (rad/s) |
| `EPS_f` | `3.35281066d-3` | Earth flattening |
| `CR` | `1.2` | Radiation pressure coefficient |
| `AM` | `0.01` | Area-to-mass ratio (m²/kg) |
| `IPSR` | `1` | SRP on/off |
| `ISHAD` | `1` | Shadow: 0=none, 1=cylindrical, 2=conical |

### `tle2opm.cfg` — TLE converter configuration

```
input/example_47944.tle.txt       ! TLE filename
47944                             ! Target NORAD ID
2025 05 01 12 00 0.000            ! Target epoch (YYYY MM DD HH MM SS.sss)
```

### `ATM.DAT` — Atmosphere density table

Two-block table: scale heights (km) and densities (kg/m³) for 291 altitude levels (60–500 km, 1 km steps to 200 km then 2 km steps to 500 km). Drag is suppressed above 500 km.

Generated by `gen_atm_jr71.F` using the real **Jacchia-71** profile with the Roberts-1971 polynomial anchors: quartic temperature 90–125 km (T(90)=183 K with zero gradient, inflection at 125 km), exponential-asymptotic above 125 km, δᵢⱼ species anchors at 125 km, ζ total-density anchor at 100 km, Aa mean-molecular-mass barometric segment 90–100 km. Reference conditions: F10.7 = F10.7B = 72, Kp = 1.0 (nighttime-minimum T∞ = 626.3 K; a static table cannot carry J71's diurnal factor). The SCH column is the local density scale height −dz/d ln ρ evaluated from the generated profile — the quantity the perigee-referenced drag model actually needs. Validated against GMAT R2026a JacchiaRoberts at the same conditions: 0.80–0.95 density ratio across 102–300 km (residual = diurnal geometry).

Replaces `gen_atm_j70.F`, whose hand-rolled single-exponential temperature profile (forced 12 K/km gradient at 90 km) ran ~127 K too warm through 90–125 km and produced a table **3.3–3.5× denser than JacchiaRoberts across the 140–200 km perigee band** (OREM issues #12/#14, 2026-07-14). Note: the earlier "drag ratio 95–97% vs NPOE" and the Phase-3 "~20% vs GMAT at Hp≈243 km" validations were both consistent with the old table because at 243 km the old table's excess was only ~5% — the defect was concentrated in the lower perigee band those comparisons didn't probe. (An initially-suspected companion "2× drag-model deficit" — OREM #25 — turned out to be an incommensurate-duration test comparison, 35 revolutions vs GMAT's 7 days; the revolution-level drag model is validated to ~1% against exact integration. A real *arc-level* drag-phase defect was found and fixed instead — see issue #24 for the `driver_KS.F` port.)

### `SW-All.csv` + `ATM2D.DAT` — Epoch-resolved space weather (optional)

Both **optional**, auto-detected at startup (no config flag). When both are present, `driver_KS` uses the real historical F10.7/Kp for each revolution's own epoch instead of the static `ATM.DAT` table; when either is absent, the legacy static-table path runs unchanged. Prints `[SW] epoch-resolved density: ENABLED` or `DISABLED` at startup so a run's actual mode is never silently ambiguous.

- `SW-All.csv`: CelesTrak's full daily solar/geomagnetic history (`curl -sL https://celestrak.org/SpaceData/SW-All.csv -o input/SW-All.csv` to refresh). Observed rows are daily; the predicted tail (PRM/PRD types) is monthly — looked up by JD binary search, never index arithmetic.
- `ATM2D.DAT`: a 2-D (altitude × exospheric-temperature) density/scale-height table, generated by `gen_atm2d_jr71.F` (same Jacchia-71/Roberts-1971 physics as `ATM.DAT`'s generator, just gridded over T∞ instead of one fixed value).
- Ported from OREM (issue #26 there): the loaders (`sw_load`/`atm2d_load`) and runtime consumers (`sw_tinf`/`atm2d_interp`) live in `src/swx.F`. Covered by `test_sw.F` (11 checks, including a hand-verified exospheric-temperature value for the 2024-05-11 G5 geomagnetic storm).

---

## 7. Output Files

### `KSROP_YYYYMMDDTHHMMSS.oem` — Trajectory (CCSDS OEM v2.0)

Uniquely named per run (UTC wall-clock timestamp). Contains `nrev × istep + 1` state vectors (position km, velocity km/s) at every integration step.

**Early termination:** Propagation stops on re-entry (alt < 80 km) or divergence (NaN). The OEM is truncated at the last valid state.

### `ksrop.opm` — Initial Keplerian elements

Written once at run start with the epoch state vector and osculating elements (a, e, i, Ω, ω, M).

---

## 8. Force Model Configuration

| Flag | Effect |
|---|---|
| `n_force(1) = 0` | Geopotential off (point mass) |
| `n_force(1) ≠ 0` | Geopotential on — degree from `ngeo_deg` |
| `n_force(2) = 0` | Solar gravity off |
| `n_force(3) = 0` | Lunar gravity off |
| `IDRAG = 1` | Atmospheric drag on |
| `IPSR = 1` | Solar radiation pressure on |
| `ISHAD = 0/1/2` | No shadow / cylindrical / conical |

---

## 9. Method

The propagator uses **KS regularisation** to remove the 1/r singularity. The 3D equations of motion are lifted to a 4D harmonic oscillator in KS space, integrated with **Runge–Kutta–Gill 4th order** using the generalised eccentric anomaly E as the independent variable.

| KS Element | Rate equation |
|---|---|
| `z(j+1)`, `z(j+5)` | State — conservative (geo + third-body) + drag + SRP |
| `z(1)` | Time — geopotential + third-body + drag + SRP |
| `z(10)` | Energy — non-conservative: dω/ds = −½(u̇·(q_drag+q_srp))/ω |

Step size is scaled by Γ = ω/ω_Kep to maintain accuracy across eccentricities.

---

## 10. Subroutines Reference

### Constants and initialisation

| Subroutine | Description |
|---|---|
| `init_constants()` | Read `input/const_new.dat` → populate `/xy/` common block (pi, d2r, r2d, amue, AU, R_Earth) |

### Coordinate transforms

| Subroutine | Description |
|---|---|
| `car2ks(x,xd,u,us,w)` | Cartesian → KS |
| `ks2car(u,us,x,xd,w)` | KS → Cartesian |
| `ks2ksr(y,u,us,E,cse,sie)` | KS regular elements → u, us |
| `car2oe(x,xd,pek)` | Cartesian → Keplerian elements (degrees) |
| `oe2car(pek,x,xd,tol)` | Keplerian → Cartesian |
| `u2uu(u,uu)` | KS index swap |
| `u2qu(u,qu,j)` | KS rearrangement for luni-solar terms |

### Force models

| Subroutine | Description |
|---|---|
| `geo_coeff(n,c_j)` | Stream EGM2008 zonal harmonics up to degree n |
| `geo_coeff_body(n,c_j,fname)` | Same, from an arbitrary gravity-coefficient file (non-Earth central bodies) |
| `geo_coeff_tess22(fname,c22,s22)` | (2,2) tesseral coefficients, normalized→unnormalized (issue #29) |
| `tess22_force(u,R,A,B,q22,tau22)` | KS-element force + time-element contribution from the (2,2) tesseral term (issue #29) |
| `cunningham_Vnm(nwork,x,y,z,r2,Vr,Vi)` | General complex solid spherical harmonics (Cunningham 1970 recursion), any degree/order (issue #30) |
| `cunningham_dVnm(n,m,nwork,Vr,Vi,dxr,dxi,dyr,dyi,dzr,dzi)` | Cartesian derivatives of `V_{n,m}` (issue #30) |
| `geo_coeff_tess_general(fname,nmax,Cnm,Snm)` | General (n,m) tesseral coefficient table loader, any degree/order (issue #30) |
| `tess_general_force(u,R,nmax,Cnm,Snm,q_gen,tau_gen)` | General (n,m) KS-element force + time-element contribution, replacing `tess22_force`'s hardcoded (2,2)-only formula (issue #30) |
| `rotate_tess_coeffs(nmax,Cnm,Snm,theta_deg,Cnm_rot,Snm_rot)` | Rotates every loaded (n,m) coefficient by its own `m·θ` GMST angle, generalizing the (2,2)-specific `A_tess`/`B_tess` trick (issue #30) |
| `force_models(n_for,ngeo,s,amoon)` | Apply force model on/off flags |
| `shadfncyl(x1,x2,x3,xs,ys,zs)` | Cylindrical shadow factor ν ∈ [0,1] |
| `shadfncone(x1,x2,x3,xs,ys,zs)` | Conical shadow factor with penumbra |
| `INTPOL(XT,YT,M1,X1,Y1)` | Linear interpolation in atmosphere table |

### Ephemeris

| Subroutine | Description |
|---|---|
| `solarnpv(dj,s)` | Sun position (geocentric, km) |
| `sun_azimuth(ai,omega,raan,alpha_s,delta_s,lambda_s)` | Sun azimuth angle w.r.t. orbital plane (Cook 1962) |
| `lunarpv(dj,tm)` | Moon position (geocentric, km) |
| `third_body_aux(deg,x,tb,R1,amu_tb,...)` | Third-body distance ratios and coefficients |
| `aLegP(n,x,P)` | Zonal Legendre polynomials |
| `aLegendreP(l,m,x)` | Associated Legendre polynomial P_lm(x) |

### CCSDS file I/O

| Subroutine | Description |
|---|---|
| `read_opm(iunit,x,xd,cal)` | Parse OPM v2.0 → state + epoch |
| `write_opm(iunit,epochstr,x,xd,pek)` | Write OPM v2.0 (state + elements) |
| `read_oem(iunit,maxpts,traj_jd,traj_x,traj_xd,npts)` | Parse OEM v2.0 ephemeris block |
| `write_oem(iunit,creation,start,stop,traj_jd,traj_x,traj_xd,npts)` | Write OEM v2.0 |
| `read_cdm(iunit,tca_cal,miss_dist,pc,x1,xd1,cov1,x2,xd2,cov2)` | Parse CDM v1.0 |
| `write_cdm(iunit,creation,tca,miss_dist,rel_speed,pc,...)` | Write CDM v1.0 |

### Time utilities

| Subroutine | Description |
|---|---|
| `cal2jd(cal,djd)` | Calendar → Julian date |
| `gmst_deg(djd,theta)` | Julian date → Greenwich Mean Sidereal Time (deg), for tesseral-harmonic longitude (issue #29) |
| `jd2epoch(djd,epochstr)` | Julian date → CCSDS epoch string |
| `parse_epoch(estr,cal)` | CCSDS epoch string → calendar |
| `utc_now_epoch(epochstr,compact)` | Current UTC → epoch string |

### Vector utilities

| Function | Description |
|---|---|
| `dotp3(x,y)` | 3-vector dot product |
| `dotp4(x,y)` | 4-vector dot product |
| `vmn(x)` | 3-vector magnitude |
| `cross(x,y,z)` | 3-vector cross product |
| `atan3(a,b)` | atan2 mapped to [0, 2π) |

### TLE reader and SGP4 (`TLEread.F`)

| Subroutine | Description |
|---|---|
| `read_tle(iunit,maxtle,...,ntle)` | Read all TLE entries (2/3-line, BOM, error handling) |
| `tle_parse1(line,...,ierr)` | Parse TLE line 1 (catalog ID, epoch, ndot, B*) |
| `tle_parse2(line,...,ierr)` | Parse TLE line 2 (orbital elements) |
| `tle_expval(str,val)` | Decode implied-decimal exponent format |
| `tle_chksum(line,ick,iok)` | Verify modulo-10 checksum |
| `tle_isline(raw,llen,lnum,istle)` | Validate TLE line format |
| `tle2sv(iyr4,eday,...,r_j2k,v_j2k,ierr)` | Complete TLE → J2000 state vector (dispatches SGP4/SDP4) |
| `tle_epoch2jd(iyr4,eday,djd)` | TLE epoch → Julian date |
| `tle_sgp4_init(...)` | SGP4 initialisation — near-Earth (period < 225 min) |
| `tle_sgp4_prop(tsince,...,r_teme,v_teme,ierr)` | SGP4 near-Earth propagation |
| `tle_sdp4_init(...,dsstate,iresfl,ierr)` | SDP4 initialisation — deep-space (period >= 225 min) |
| `tle_sdp4_prop(tsince,...,dsstate,iresfl,r_teme,v_teme,ierr)` | SDP4 deep-space propagation |
| `tle_teme2j2k(r_teme,v_teme,djd,r_j2k,v_j2k)` | TEME → J2000 frame rotation |

### Integrator

| Function | Description |
|---|---|
| `rkgil(n,y,f,x,h,nt)` | Runge–Kutta–Gill 4th-order (4-stage, fixed step) |

---

## 11. Performance

| Configuration | Throughput | Wall time (10 rev × 360 steps) |
|---|---|---|
| Two-body | ~140,000 steps/s | ~31 ms |
| J2 geopotential | ~62,000 steps/s | ~57 ms |
| Degree-50 geopotential | ~55,000 steps/s | ~65 ms |
| With drag | ~108,000 steps/s | ~34 ms |

EGM2008 streaming read: O(n²) lines for degree n (J2 = 3 lines, not 2.4M).

---

## 12. Known Issues

- EGM2008 file (~231 MB) not included; set `ngeo_deg = 0` for point-mass gravity.
- SRP is cannonball only; no tesseral harmonics or geometry-dependent variations.
- `solarnpv`/`lunarpv` are analytic ephemeris series (Montenbruck & Gill
  low-precision, since 2026-07-12), not JPL DE405: ~0.10% (Sun) / ~0.11%
  (Moon) position error vs DE405 over a 60-day sweep
  (`scratch_gmat/ephemeris_check_gmat.script`). This bounds achievable
  third-body dynamical agreement at roughly the 0.5 km/rev level (GTO,
  Moon-dominated); tighter would require DE-series ephemerides.
- Third-body Legendre truncation: Sun `nsun_deg=2`, Moon `nmoon_deg=3`
  by convention (thesis Ch. 4); the n=4 lunar term is ~1% of the lunar
  effect at GTO apogee.

---

## 13. Revision History

| Date | Change |
|---|---|
| 2018-06-15 | Initial program, J2 only |
| 2018-06-16 | Revolution-by-revolution output; nth-degree Legendre polynomial |
| 2018-09-13 | nth-degree geopotential up to 2190 (EGM2008 Jn) |
| 2021-07-21 | Legendre polynomial subroutine (`aLegP`) |
| 2026-06-06 | Luni-solar third-body perturbations |
| 2026-06-06 | Atmospheric drag (oblate exponential model, ATM.DAT) |
| 2026-06-06 | CCSDS OEM/OPM v2.0 I/O; `input.opm` for initial conditions |
| 2026-06-06 | Performance optimisation: streaming geo_coeff, removed 115 MB static array |
| 2026-06-07 | Drag model ported from KSJLSDNP2.F (oblate, co-rotating, perigee-referenced) |
| 2026-06-07 | Early-termination on re-entry (< 80 km) or divergence (NaN) |
| 2026-06-07 | Code cleanup: removed dead subroutines, consolidated `third_body_aux` |
| 2026-06-07 | Fixed `cal2jd` ~84-day epoch offset bug; fixed input.DAT UTF-8 BOM crash |
| 2026-06-07 | CCSDS CDM v1.0 I/O (`read_cdm`/`write_cdm`) |
| 2026-06-08 | Solar radiation pressure: cannonball + cylindrical/conical shadow models |
| 2026-06-20 | Fixed double-Gam step-size bug in KS integrator |
| 2026-06-20 | Added `TLEread.F`: TLE reader + SGP4 + TEME→J2000 (147+78 tests) |
| 2026-06-20 | Added `test_initial_conditions.py` (110 checks, 10 orbital regimes) |
| 2026-06-21 | Fixed 11 bugs in Subrouts.F (type mismatch, dead code, precision, uninitialized vars) |
| 2026-06-21 | Added `tle2opm.F`: TLE-to-OPM converter; generated NORAD 47944 OPM (2025-05-01) |
| 2026-06-21 | Refactored constants: `init_constants()` reads `const_new.dat`; unified common block |
| 2026-06-22 | Added SDP4 deep-space theory to `TLEread.F`: `tle_sdp4_init`/`tle_sdp4_prop` with lunar-solar secular gravity, synchronous and half-day resonance handling, long-period periodics; `tle2sv` now dispatches SGP4 (period < 225 min) or SDP4 (period >= 225 min) transparently; deep-space state carried in `dsstate(50)` array with `iresfl` resonance flag; 78 new tests (156 total in `test_tle2sv.F`) covering GEO, Molniya, GPS orbit types, resonance detection, propagation, edge cases (backward prop, 30-day duration, e=0.85, retrograde, polar, critical inc, angle sweeps, energy/momentum conservation, SGP4/SDP4 boundary continuity), and full pipeline validation |
| 2026-06-22 | Added `sun_azimuth` subroutine to `Subrouts.F`: computes Sun azimuth angle (Λ_S) w.r.t. spacecraft orbital plane per Cook (1962); 10 tests in `test_subrouts.F` (apsidal alignment, quadrant sweep, inclination/RAAN/ω rotation, polar orbit, range validation) |
| 2026-06-22 | Added automated lint checks: `lint_check.sh` (line length, single-precision intrinsics, common block consistency, tabs, trailing whitespace) and `test_all.sh` (unified runner for lint + all test suites); fixed 12 source lines exceeding 72-column F77 limit |
| 2026-06-24 | Fixed atmospheric drag crash for HEO orbits (Issue #16): replaced hardcoded 500 km perigee guard with `ALT_atm` table-bounds check; added `H_dg≤0` safety after INTPOL; added exponential overflow clamp (`|arg|>500 → 0`). Prevents NaN cascade from INTPOL returning zero scale height for out-of-range altitudes |
| 2026-07-04 | Fixed `car2oe` NaN at perigee/apogee: all `dacos()` calls now clamp their argument to [-1, 1] via `dmax1(-1.d0,dmin1(1.d0,...))`. Floating-point dot-product can exceed ±1 by ε at apsides, causing `dacos(>1) = NaN` that propagated through eccentric anomaly → drag density → full state divergence. |
| 2026-07-07 | Replaced `ATM.DAT` with Jacchia-70 multi-species diffusive equilibrium table (F10.7=72, Kp=1.0, T∞=640 K). Added `gen_atm_j70.F90` generator. Drag ratio vs NPOE: 48% → 95–97%. |
| 2026-07-11 | Fixed `aLegP` (`Legendre.F`) buffer overflow: ignored its own degree argument, always computed a full 50×50 Legendre grid into a 50-element output array and a 36×36 scratch array — a ~50x out-of-bounds write on every call (including calls with `n` as small as 2), silently corrupting adjacent stack variables. Found via independent GMAT cross-validation (see `scratch_gmat/`): a J2-only propagation diverged from an independent RK89 Cartesian propagator by ~150 km over one orbit despite the pure two-body case agreeing to sub-mm/s. Rewrote to honor `n` and use correctly-sized buffers. |
| 2026-07-11 | Fixed `aleg`/`sleg`/`oleg` off-by-one degree (`driver_KS.F`): `aLegP(ngeo_deg,...)` fills `aleg(1..ngeo_deg+1)`, but the oblateness force and time-element formulas need `aleg(i+2)` up to `i=ngeo_deg`, i.e. `aleg(ngeo_deg+2)` — one degree beyond what was requested. That slot was uninitialized. The pre-fix `aLegP` buffer overflow masked this (it always overfilled the array regardless of `n`), so fixing the overflow first exposed this separate, pre-existing bug. Fixed by requesting `ngeo_deg+1`/`nsun_deg+1`/`nmoon_deg+1` at both call sites. Validated against GMAT (independent RK89 propagator, EGM96 zonal-only J2, `scratch_gmat/`): a near-polar test orbit's classical nodal-regression rate went from 10x too large (0.692 vs GMAT's 0.0693 deg/orbit) to matching within 0.15% (0.0692 deg/orbit), averaged over 20 orbits. An earlier attempted fix to the `qj` force-formula's algebra itself was investigated and found unnecessary — reverted (see git history) after this off-by-one turned out to be the actual root cause; `qj` matches its long-standing legacy-heritage form and independently validates correctly against GMAT once fed complete Legendre data. |
| 2026-07-11 | Fixed `Tau_geo` (time-element numerator, `driver_KS.F`): was missing the `amue` factor present on every other term in the same `z(1)` sum, and had the wrong sign. Derived the correct closed form from the KS regular-elements time-element rate equation (Sellamuthu 2018 PhD thesis, eq. 2.56: `τ* = [μ-2rU]/(8w³) - r/(16w³)[u·∂U/∂u] - ...`), which reduces per zonal term to `+amue*(n-1)*c_j(n)*Re(n)*aleg(n+1)*ObyR(n)` — confirmed against the pre-fix formula numerically to 9 significant digits (exactly `-1/amue` times the old value). Since `Tau_geo` only affects OEM epoch timestamps, not propagated states (verified by isolation test, see 2026-07-11 entries above), the fix shows up as corrected time-labeling rather than a position change: re-running the GMAT J2–J20 cross-check with the corrected epochs closed the remaining ~18 km gap down to **~80 m** (a ~225x improvement), consistent with the ~2 s of previously mislabeled time this fix removes. |
| 2026-07-11 | Fixed `qsun`/`qmoon` (third-body force, `driver_KS.F`): the `slambda`/`olambda`-based closed form gave results 8 orders of magnitude off the correct value, with a wildly inconsistent per-component ratio (not a clean sign/scale error) — traced to the third-body potential's radial power law being `+n` (grows with `r`, a multipole expansion about the third body) rather than zonal's `-(n+1)` (falls off from a point source), which flips which Legendre recursion identity the chain-rule derivation needs (`P_(n-1)`-based, not `P_(n+1)`-based). Re-derived from first principles for a general (non-axis-aligned) third-body direction using the KS bilinear identity `L(u)ᵀx = r·u`, and verified the new closed form against `L(u)ᵀ` of a finite-difference gradient of the thesis's own third-body potential (eq. 4.1, confirmed matching sign convention) to ~1e-9 relative precision at multiple points, multiple degrees, and both `car2ks` branches. New formula: `qsun(j) = RbyRs(i)/(R(1)·den_s)·[i·(sleg(i+1)-cphi_s·sleg(i))·u(j) + (i+1)·(cphi_s·sleg(i+1)-sleg(i+2))·sq/Rs(1)]` (same structure for `qmoon`/moon). Separately, `Tau_3body1`/`Tau_3body2` were checked against the same eq. 2.56 framework (via Euler's homogeneous-function theorem, `x·∇V=m·V`, avoiding a full re-derivation) and found to **already be correct** — no change needed there. Dynamical GMAT sanity check (Sun-only, GTO-like orbit, one ~10.5h orbit): KSROP-vs-GMAT gap grows smoothly to ~6.5 km, consistent in scale with the ~0.6% difference between KSROP's compact analytic solar ephemeris (`solarnpv`) and GMAT's DE405 (quantified separately, see Known Issues) — not attributable to the force formula itself, which is independently validated to 1e-9 precision decoupled from ephemeris source. |
| 2026-07-12 | GMAT validation campaign completed (issues #19-#23): **Phase 3 drag** — unit audit of the drag acceleration exact (the `1e-7` factor combines with `DEN_atm`'s `1e10` read-scaling to the textbook `ρv²/(2·BN)`); secular SMA decay at Hp≈243 km, BN=50: KSROP -2.01 km/day vs GMAT JacchiaRoberts (matched F10.7=72/Kp=1.0) -2.4 to -2.6 km/day — ~20% apart with both implied densities physical, within normal inter-atmosphere-model spread. Sampling must be phase-locked (periapsis) to avoid J2 osculating-SMA aliasing. **Phase 5 full force** (J2-J20 + Sun + Moon + drag + SRP conical, TC4 GTO, 2 revs ≈ 21 h): 0.3 km @ 1.4 h growing to 76 km @ 20.9 h, error concentrated along-track at perigee (≈7 s accumulated timing), consistent with the documented lunar-ephemeris (~3.6%) and drag-model (~20%) gaps — no evidence of further force-model defects. Scripts: `scratch_gmat/KSROP_drag_crosscheck.script`, `KSROP_fullforce_crosscheck.script`. |
| 2026-07-12 | **Found & fixed the two real third-body bugs** (supersedes the 2026-07-11 `qsun`/`qmoon` entry and the Phase-4/5 "ephemeris gap" attributions): (1) `third_body_aux`'s `deg` dummy argument fell under `implicit double precision (a-h,o-z)` while every caller passes an integer — the bit pattern read as ~1e-323, its power-series do-loop ran zero times, `RbyRtb(2..)` stayed 0.0, and **the entire third-body force was exactly zero from the 2026-06-07 refactor until now** (found via a GMAT ON/OFF sensitivity test: sun+moon on vs off was bit-identical). Fixed with `integer deg`. (2) The 2026-07-11 `qsun`/`qmoon` rewrite used the wrong EOM convention — `L(u)ᵀ(−∇shape)` alone, when the KS-elements equations (thesis eq. 2.55/2.59) consume `2[U/2·u + (r/4)∂U/∂u]`, i.e. `qsun = shape·u + r·L(u)ᵀ(∇shape)` — missing a factor ~r (~1e4 at GTO). Corrected form verified to machine precision (2e-16) against the independently validated KSJLSDNP.F n=2 sun term (CNTS1/CNTS2 structure) and to ~1e-8 vs finite-difference at n=2,3. `Tau_3body1/2` confirmed correct under the same convention (algebra + dynamics). |
| 2026-07-12 | Upgraded `solarnpv`/`lunarpv` to Montenbruck & Gill low-precision analytic series (J2000 frame; the M&G solar series is already J2000-referenced — do NOT add the −1.3972°·T precession term there, it belongs only in the lunar series). Sun: fixed-1-AU distance replaced with the elliptic form — 0.6% → **0.097%** vs DE405; Moon: 7-term Cartesian series replaced with truncated Brown theory — 3.6% → **0.109%** vs DE405 (60-day sweep, 13 epochs). |
| 2026-07-12 | Re-validated third-body + full force vs GMAT after the fixes: **Sun-only GTO 1 rev: 1.2 m** (was 6.5 km); **Moon-only 1 rev: 0.46 km** (bounded by lunar ephemeris ~0.11% + n=3 truncation); **full conservative (J2-J20 + Sun + Moon + SRP), 2 revs: 1.9 km** (was effectively ~76 km); full force with drag: 55.9 km at the Hp≈178 km perigee — ~96% of which is the Jacchia-70-vs-JacchiaRoberts drag-model difference (documented model choice, see issue #21), not a code defect. Regression: 528/528. Scripts: `scratch_gmat/KSROP_sunonly_crosscheck.script`, `KSROP_moononly_crosscheck.script`, `KSROP_nodrag_crosscheck.script`. |
| 2026-07-14 | **Issue #24 fixed**: ported OREM v1.18's drag-phase fix to `driver_KS.F` — the drag density's eccentric anomaly is now read from the true state (`pek(7)`, per-stage `car2oe`) instead of the analytic sweep `DE_dg=(VIPP·π−EA₀)/istep`, whose `VIPP=4` branch half-rated the density phase whenever a revolution started past EA=π, dephasing the density peak from the true perigee along long decay arcs. Bit-identical on perigee-anchored windows (why the Phase-3 GMAT drag validation never saw it); on OREM's re-entry arcs the same fix moved RPE from −70..−97% to bracketing zero. Full suite passes (67+17+147+156+21). |
| 2026-07-14 | Replaced `gen_atm_j70.F` with `gen_atm_jr71.F` and regenerated `input/ATM.DAT`: the old generator's single-exponential temperature profile ran ~127 K too warm through 90–125 km, making the table **3.3–3.5× denser than GMAT's JacchiaRoberts across the 140–200 km perigee band** (quantified by OREM's GMAT density probe; drove OREM's re-entry predictions ~4–5× early — OREM issues #12/#14). New generator implements the real Jacchia-71 profile with the Roberts-1971 polynomial anchors; generated table tracks GMAT JR at 0.80–0.95 over 102–300 km at matched static conditions. SCH column now = local −dz/d ln ρ. No propagator source changed; the KSROP test suites do not read ATM.DAT (only `driver_KS` does at runtime), so test counts are unaffected. Follow-up (updated same day): the suspected "~2× drag-model deficit" (OREM #25) was an incommensurate-duration test comparison — the revolution-level drag model is validated to ~1% vs exact integration. However a real **arc-level drag-phase defect** was found: the analytic EA sweep's `VIPP=4` branch half-rates the density phase whenever a revolution starts past EA=π, dephasing the density peak from the true perigee along decay arcs. Fixed in OREM's `propagate_ks` (RPE −70..−97% → ensemble +11%); `driver_KS.F` carries the identical code — port tracked as **issue #24**. |
| 2026-07-18 | **Packaged for production (v2.0.0)**: added GitHub Actions CI (`ci.yml`, gfortran build + full 528-check suite on every push) and a release workflow (`release.yml`, builds and publishes Linux binaries on `v*.*.*` tags). Bringing the build up on gfortran/Linux (previously only ever compiled with Intel `ifx` on Windows) surfaced and fixed 4 real portability/correctness bugs, all the same "implicit-typing trap" class as earlier fixes: `cn0` (array-size parameter silently DOUBLE PRECISION, `driver_KS.F`), `aLegP`'s scalar/rank-1 argument mismatch (`Legendre.F`), `tle_ds_lpper`'s `mp` (mean anomaly) silently INTEGER across a call boundary where the caller declares it `double precision` (`TLEread.F`), and `iE` silently INTEGER against a REAL-format debug-dump write that crashed `driver_KS.exe` outright on gfortran (`driver_KS.F`). Also fixed a real Linux-only bug: `input/const_new.DAT`/`input.DAT` were tracked with uppercase extensions while every `open()` call (and both Python integration tests) has always addressed them lowercase -- silently masked by Windows/NTFS case-insensitivity; fresh Linux checkouts renamed to match. Repo hygiene: fixed `.gitignore` (the old `output/regular.out` pattern never matched the actual `output/KSROP_*_Regular.out` naming), removed ~80 accumulated untracked run-output files and orphaned scratch work, fixed the `Makefile` (default `FC=gfortran`, added missing `tle2opm`/`test_bugs`/`test_tle2opm` targets), and fixed `test_all.sh`'s hardcoded `.exe` executable detection and Python-integration-test summary parsing, both of which silently reported 0 tests run as "ALL TESTS PASSED" prior to this fix. |
| 2026-07-24 | **Migrated to fpm (Fortran Package Manager) packaging, v2.1.0**: restructured into `src/`(library)/`app/`(executables)/`test/`(test programs) with a new `fpm.toml` (`[fortran] source-form = "fixed"` — fpm defaults `.F` to free-form, which breaks column-1 comments immediately). Old Makefile/`test_all.sh`/`lint_check.sh`/`ci.yml` paths updated to match, not removed — both build paths work. Added a second CI job building+testing via `fpm`+gfortran on Linux, since that's the actual path any consumer resolving KSROP as a dependency (e.g. OREM) exercises. Tagged `v2.1.0`; OREM now depends on this tag via `fpm.toml`'s git dependency instead of hand-copied files. |
| 2026-08-07 | Added `geo_coeff_body(n,c_j,fname)`: `geo_coeff` refactored into a thin Earth-default wrapper over a new file-parameterized entry point, so a consumer depending on KSROP as a library (e.g. a future lunar- or Mars-centered driver, see this repo's issues #26 and #27) can load its own gravity-coefficient file in the same EGM2008-row format without duplicating the parsing logic. No behavior change for existing callers. 3 new tests (`test/fixture_geo_coeff.dat`), 542/542 total. |
| 2026-08-08 | **Tesseral harmonics, Phase 1 of 6 (issue #29, v2.4.0)**: added `gmst_deg(djd,theta)` (Greenwich Mean Sidereal Time, IAU 1982/Meeus) — the body-fixed rotation angle needed for tesseral-harmonic longitude (`λ = atan2(x(2),x(1)) − θ`), evaluated directly per-step from the already-tracked Julian date rather than propagated from a stored reference epoch. Grounded in M. Xavier James Raj's 2007 PhD thesis (Mahatma Gandhi University, supervisor Dr. R.K. Sharma) §3.2-3.3 for the geopotential/coordinate-chain physics only — not its §3.4 force-law derivation, which is canonical/Hamiltonian and does not match KSROP's own variation-of-parameters KS regular elements (Sellamuthu thesis eq. 2.55/2.56/2.59). 3 new tests, 545/545 total. Remaining phases (associated-Legendre grid, tesseral coefficient loader, force/time-element derivation, validation, KSROP-Lunar/Mars rollout) tracked in issue #29. |
| 2026-08-08 | **Tesseral harmonics, Phases 2-4 of 6 (issue #29, v2.5.0)**: `aLegendreP` verified directly against closed forms for the (2,0)/(2,1)/(2,2) set (Phase 2, no code change — it was already correct); `geo_coeff_tess22(fname,c22,s22)` added, reading the (2,2) sectorial row and converting normalized to unnormalized coefficients via N₂,₂=√(5/12), cross-checked against known Earth literature values using the existing test fixture's real EGM2008 (2,2) row (Phase 3); `tess22_force(u,R,A,B,q22,tau22)` added — the KS-element force and time-element contributions from the (2,2) tesseral potential (Phase 4). The force-law derivation was genuinely hard: an initial recipe (from a paraphrased citation) was off by a large factor; reverse-engineering numerically from the existing zonal formula gave an inconsistent ratio; the published source (Sellamuthu & Sharma, *J. Spacecraft and Rockets* 55(5):1282-1288, 2018 — the actual primary source behind this repo's own KS regular-element equations) was located and read directly, but even its explicit J2 formula didn't reproduce cleanly by hand. Resolution: an exact symbolic (sympy) derivation, verified against KSROP's own existing zonal `q(j)`/`Tau_geo` at independent exact rational points for two different degrees, which also caught two real errors before they reached this formula — a wrong prefactor, and a missing `1/r` in the (2,2) potential itself (the thesis's own eq. 3.2.1 extraction was short one power of `r` relative to the standard geopotential form, caught by cross-checking against Alvarellos (2009), *J. Astronautical Sciences* 57(4):701-715, whose eq. 1 has an explicit `μ/r` leading factor). The resulting Fortran was further verified against a Fortran-native finite-difference gradient (~1e-9 relative agreement) before being trusted. Not yet wired into `driver_KS.F`'s propagation loop or GMAT-validated — that's Phases 5-6. 9 new tests, 554/554 total. |
| 2026-08-08 | **Tesseral harmonics, Phase 5 start of 6 (issue #29, v2.6.0)**: `tess22_force`/`geo_coeff_tess22`/`gmst_deg` wired into `driver_KS.F`'s live propagation loop — loaded once (auto-detect, defaults to inert if no (2,2) row is present), the GMST-derived `A_tess`/`B_tess` refreshed at the same per-stage cadence as the existing third-body ephemerides, and the (2,2) contribution added to both `q(j)` and `Tau_geo` alongside the existing zonal terms. Verified genuinely active (not a silent no-op) via a before/after comparison run — a 3-revolution LEO propagation with `ngeo_deg=2`, real Earth EGM2008 (2,2) coefficients vs. zeroed ones, all other force terms off — showing a ~27 m position divergence by the end of the run, a physically reasonable magnitude for a term two orders of magnitude weaker than J2. Full regression suite unaffected (156/156 fpm test programs, 554/554 checks, lint clean). Remaining: external validation checkpoints (thesis Table 3.5, Alvarellos Mars libration, GMAT cross-check) and KSROP-Lunar/Mars rollout — Phases 5 (rest)-6. |
| 2026-08-08 | **Tesseral harmonics, Phase 5 external validation, rest of 6 (issue #29)**: ran the three remaining checkpoints. XJR thesis Table 3.5 turned out not to be a valid direct comparison target — it reports the combined effect of ~20+ tesseral coefficient pairs up to degree 6/10/19, not (2,2) alone, so a smaller-magnitude, mixed-sign result on an isolated-(2,2) test was expected once this was recognized (the scope mismatch should have been caught from Table 3.1B before running the test, not after). Alvarellos (2009) Mars was a strong pass: verified symbolically (sympy) that `Vhat22` is algebraically identical to Alvarellos eq. (4)'s (2,2) term at the equator (up to an overall sign, consistent with KSROP's own already-validated force-function convention), and reproduced the paper's own derived λ₂₂=74.7° and J₂₂=6.31e-5 from its Table 1 Mars coefficients, exact to stated precision. The GMAT cross-check (two new scripts, `scratch_gmat/KSROP_tess22_{off,on}_crosscheck.script`) gave an order-of-magnitude-consistent but phase-differing result (~52m GMAT vs. ~35m KSROP over 3 LEO revs, opposite sign on the dominant component) — diagnosed as `gmst_deg` being a mean (IAU 1982) sidereal-time formula versus GMAT's precise Earth-orientation kernels, which matters for short/non-resonant windows but far less for the ~130-day Mars libration timescale this capability was actually built for; documented as a known precision limitation rather than fixed now. Net: the (2,2) force law itself is verified exactly twice (Phase 4 symbolic match to KSROP's own validated zonal formula; Alvarellos match here) and at the right order of magnitude against GMAT. Phase 6 (KSROP-Lunar/Mars rollout) remains. |
| 2026-08-08 | **General (n,m) tesseral/mascon harmonic support, issue #30 (v2.7.0)**: scoped separately from #29 (which shipped only the (2,2) term) per explicit instruction. `cunningham_Vnm`/`cunningham_dVnm` (`src/Cunningham.F`) implement Cunningham (1970)'s recursive complex solid spherical harmonics and their Cartesian derivatives — arbitrary degree/order, no polar singularity, the same method GMAT itself uses internally. Chain-ruled through KSROP's own closed-form KS position Jacobian (`∂x_k/∂u_j`, from the existing `car2ks`/`ks2car` bilinear map) into the already-general `q(j)=V·u(j)+(r/2)∂V/∂u_j` recipe — no new KS-space derivation needed per term, unlike #29's one-off symbolic approach. `geo_coeff_tess_general` generalizes the coefficient loader to a full (n,m) table; resolved a real, non-obvious sign-convention split in the process (`geo_coeff`'s zonal loader applies an extra minus — the traditional "J_n≡-C_n0" relabeling — that `geo_coeff_tess22`'s tesseral loader does not; the general loader replicates both, branching on m=0 vs m>0). `tess_general_force` combines all of the above; verified to reproduce `tess22_force` exactly (diff ~1e-14, floating-point noise) at n=m=2, the acceptance gate for trusting every other (n,m). 19 new tests (`test_cunningham.F`): 9 against Cunningham's own Table I closed forms, 3 finite-difference derivative checks, 2 coefficient-loader cross-checks, 5 end-to-end force-law matches against `tess22_force`. Not yet wired into `driver_KS.F` or rolled out to KSROP-Lunar/KSROP-Mars — tracked in issue #30. |
| 2026-08-08 | **General (n,m) tesseral support wired into driver_KS.F, issue #30 (v2.8.0)**: `rotate_tess_coeffs` added, generalizing the (2,2)-specific GMST-rotation trick (`A_tess`/`B_tess`) to every loaded term — `tess_general_force` operates on inertial coordinates, so each order-`m` coefficient pair needs its own `m·θ` rotation, not the fixed `2θ` the (2,2)-only case used; verified by direct trigonometric expansion and 4 new tests confirming it reduces exactly to the old formula at m=2 and correctly leaves m=0 unrotated. `driver_KS.F`'s #29 (2,2)-only wiring replaced outright with the general path, reusing `ngeo_deg` as the shared max degree (capped at a new `ntess_cap=10` constant, independent of `ngeo_deg`'s own much larger practical zonal range, since the Cunningham recursion is O(n²) per step vs. the zonal path's O(n) — revisit that constant if a future consumer needs higher tesseral degree). **Caught and fixed a real bug during the integration test**: the general coefficient loader was loading `m=0` (zonal) rows too, double-counting J2 against KSROP's existing separate zonal pathway — produced a perturbation orders of magnitude too large before being traced and fixed (loader now always skips `m=0`, which also simplified away the sign-convention split documented in the previous entry — that reasoning was correct as a `(2,2)`-isolated verification exercise, but the m=0 branch was never actually needed once zonal rows are properly out of scope for this loader). After the fix, the general wiring reproduces the old (2,2)-only wiring's trajectory almost exactly (cm-level differences, consistent with floating-point operation-order noise between the two independent code paths) on the same before/after LEO comparison from #29 Phase 5, both with real (2,2) data and with it zeroed. Full regression green (585/585 total, +6 tests since the previous entry), lint clean, both Python integration suites (`test_driver.py` 10/10, `test_initial_conditions.py` 110/110) pass unchanged. Remaining: rollout to KSROP-Lunar/KSROP-Mars — Phase 6. |
| 2026-07-25 | **Ported epoch-resolved space weather from OREM (issue #26 there)**: KSROP already had the `ATM2D.DAT` generator (`gen_atm2d_jr71.F`) but nothing in `driver_KS.F` ever consumed it at runtime — only OREM's downstream fork did. Added `src/swx.F` (`sw_load`/`atm2d_load`/`sw_tinf`/`atm2d_interp`, ported unchanged from OREM's already-validated code) and wired an auto-detect hook into `driver_KS.F`'s per-revolution drag setup: when `input/SW-All.csv` and `input/ATM2D.DAT` are both present, the reference density is evaluated at the real historical F10.7/Kp for that revolution's epoch instead of the static `ATM.DAT` table; absent either file, the legacy path runs unchanged. Kept the loader+consumer subroutines together in one `src/` file (unlike OREM's own split, which puts the consumers in `propagate_ks.F`) since `driver_KS.F` is a standalone `program` in `app/`, not a shared library file — fpm only auto-links `src/` into every target, so consumers placed in `app/driver_KS.F` would have been unreachable by `test_sw`'s separate executable. New `test/test_sw.F` (11 checks, ported from OREM's test suite, including a hand-verified exospheric-temperature value for the 2024-05-11 G5 storm), added to `fpm.toml`, `test_all.sh`, the `Makefile`, and `lint_check.sh`. 539/539 tests pass (528 existing + 11 new), confirmed on both `ifx` (full end-to-end run, `[SW] epoch-resolved density: ENABLED` fires correctly) and via lint. |
| 2026-08-09 | **Independent second derivation of the general (n,m) tesseral force (`src/LegendreTess.F`), cross-validated against Cunningham (issue #30)**: evaluated an external from-scratch KS-regularized geopotential derivation (`EarthGravityPotential_KS (4).wl`, `E:\Research\1. R&D\Mathematica`) as a candidate alternative/cross-check to `Cunningham.F`'s Cartesian solid-harmonic recursion — a classical spherical-coordinate (latitude/longitude) associated-Legendre expansion substituted directly into KS `u`-variables instead. Verification found and fixed **two real bugs in the source derivation** before it could be trusted: (1) the general `(n,m)` associated Legendre recursion's `0<m<n` branch was structurally wrong (confirmed against `scipy.special.lpmv`, correct only by coincidence at the sectorial/sub-sectorial edges, error growing from O(1) at (3,1) to O(1e5) by (8,1)); (2) longitude was computed from the wrong coordinate pair (Y,Z instead of X,Y), inconsistent with the same derivation's own latitude convention (Z as the polar axis, matching this repo's `aLegP(ZbyR)` usage) — confirmed via a pure (2,2) test case disagreeing with the already-validated Cunningham potential by 50-100+ units (even wrong sign) using the source formula, matching to 1e-14 once corrected. With both fixed, `tess_legendre_force` (same `q(j)`/`tau` interface as `tess_general_force`, drop-in comparable) reproduces it to machine precision (worst case 5.6e-15) across `nmax`=2..6, random full-triangle `Cnm`/`Snm`, random `u`. New `test/test_legendre_tess.F` (26 checks: 6 closed-form associated-Legendre checks including the exact interior case that exposed bug #1, 20 end-to-end force cross-checks), added to `fpm.toml` and the `Makefile`; 611/611 total, confirmed on both `ifx` and `gfortran`. Also corrected stale `ALGORITHM.md` §5/§8/§9 tesseral documentation along the way (it had never been updated past issue #29's (2,2)-only state despite #30's general-(n,m) rollout already being complete and shipped to KSROP-Lunar/KSROP-Mars). |
| 2026-08-09 | **Wired `tess_legendre_force` into `driver_KS.F` as a selectable alternative to Cunningham**: new optional 6th `input.dat` line, `ITESS_METHOD` (0=Cunningham, default — every pre-existing `input.dat` without this line keeps running exactly as before; 1=LegendreTess), guarded by its own `read(...,end=98)` so backward compatibility is unconditional. Verified genuinely active and correct via a real before/after propagation (`ngeo_deg=50`→`ntess_use=10`, real EGM2008 data, 3601 steps): OEM output identical to the file's own printed precision, and the higher-precision `10f20.9` KS-element debug dump byte-for-byte identical between the two independently-derived force paths. 611/611 regression tests unaffected (the library itself didn't change, only `app/driver_KS.F`'s own call site). |
| 2026-08-09 | **Simplified to a single active geopotential, per explicit user direction ("keep only one")**: removed the `ITESS_METHOD` switch added the same day — `driver_KS.F`'s general (n,m) tesseral call site now unconditionally calls `tess_legendre_force`, the classical-Legendre derivation. `Cunningham.F`/`tess_general_force` remain in the library, unchanged, as a standing independent cross-check (`test_cunningham.F` still runs, `test_legendre_tess.F` still cross-validates the two to machine precision) — not deleted, since KSROP-Lunar and KSROP-Mars both depend on `tess_general_force` directly and Cunningham's Cartesian solid-harmonic recursion has no polar-latitude singularity, a real robustness property LegendreTess's classical latitude/longitude formulation lacks (see `ALGORITHM.md` §9). Re-verified via a fresh before/after propagation that the unconditional call produces the same trajectory as before. 611/611 regression unaffected. |
| 2026-08-09 | **Raised `ntess_cap` 10→72, per explicit user direction**: the standard truncated-gravity-model resolution (e.g. EGM96 72x72) used throughout operational precision orbit propagation. Verified with a real propagation at `ngeo_deg=72` (real EGM2008 data, 3601 steps): completes in 5.4s (this recursion is $O(n^2)$ per force evaluation, so real cost, but far from a practical concern), output sane (no NaN, reasonable trajectory magnitudes). 611/611 regression unaffected (`ntess_cap` is a `driver_KS.F`-local constant, not touched by any library test). |
| 2026-08-09 | **Energy-element `dU/dt` term added for the general tesseral potential, `feature/tesseral-energy-time-dependence` branch (isolated, not merged)**: `w*` (`z(10)`) previously contained only non-conservative (drag/SRP) terms — correct while every conservative force was static in the inertial frame, but the general tesseral potential is not (its coefficients get re-rotated by Earth's rotation angle every step). Sellamuthu 2018 PhD thesis eq. 2.54 confirms the gap directly: *"The first term...vanishes when U does not explicitly depend on t"* — no longer true here. Derived `dU/dt` two independent ways before writing any Fortran (sympy, zero residual; finite difference holding `u` exactly fixed and perturbing only `theta`, matched to `1.26e-9` relative) and implemented it reusing the potential's already-computed `dV/dlambda` block — no new recursion. **Independent literature corroboration found during review**: a second reference (same §3.2-3.3 already cited in `driver_KS.F`'s header) confirmed `V`/`dH/dphi`/`dH/dlambda` match exactly and stated the *same* incorrect Legendre recursion already found and fixed in the `.wl` source, reverified with exact rational arithmetic — strong evidence the `.wl` bug is inherited from this literature, not a one-off. **Primary-source-verified that `Tau_geo` needs no direct edit**: thesis eq. 2.56's own final term matches `driver_KS.F`'s pre-existing `z(1)` structure exactly, so once `z(10)` includes `dU/dt`, `z(1)` picks it up automatically — locked in with a permanent regression test (`test_dvdt_tess.F`, 7 new checks, 618/618 total). **Empirical validation, with a corrected hypothesis**: a 100-revolution conservative-only propagation at `ngeo_deg=72`, compared against the pre-fix commit via a temporary git worktree, showed `w` frozen to machine precision before the fix (confirming `z(10)` genuinely had zero conservative contribution) and varying with a bounded, oscillatory character after it (25 sign changes across 100 revs, no secular drift) — the *correct* physical signature (periodic energy exchange with a rotating asymmetric field, the same mechanism behind GEO/areostationary resonance), not the "more constant" outcome originally hypothesized in the plan. Third-body (`qsun`/`qmoon`) has the identical gap by the same mechanism, explicitly scoped out by user decision — filed as a known follow-on, not started. |
| 2026-08-09 | **Third-body (Sun/Moon) `dU/dt` extended onto the same branch**: the tesseral-only scoping above was reversed per user direction — since KSROP does full numerical integration rather than an approximate/analytical theory, `qsun`/`qmoon`'s identical explicit-time-dependence gap (`ts`/`tm` move via `solarnpv`/`lunarpv`) was brought in rather than left incomplete. Unlike the tesseral case (a pure rotation), the third-body disturbing function depends on time through the third body's real orbital motion — both `Rtb(t)` and `cos(phi(x,t))` vary. Derived `third_body_dvdt` (`src/Subrouts.F`) symbolically via the product rule through `rho_dot`/`c_dot`, fed by a central-difference third-body velocity (`ts_dot`/`tm_dot`, computed once per step). Verified against an independent finite-difference ground truth using the actual production ephemeris (`test_dvdt_3body.F`, 9 new checks: 4 position/epoch cases × Sun+Moon + a degree-guard check), summed into `z(10)` alongside the tesseral term. **Found and fixed a real implicit-typing bug along the way** (repeat of the `deg`-typing class, this time in the new test file itself): `leg`/`leg_p`/`leg_m` all start with `l`, silently INTEGER under `implicit double precision (a-h,o-z)` (which never touches the `i`-`n` range), overflowing the stack when `aLegP` wrote real DOUBLE PRECISION values into them — reproduced only on `gfortran`, isolated via a minimal repro that initially (wrongly) implicated `aLegP` itself before the actual caller-side type defect was found. `aLegP` and its existing call sites were unaffected. 627/627 total, confirmed on both `ifx` and `gfortran`. |
| 2026-08-09 | **Third-body `dU/dt` reverted; geopotential unified into a single general (n,m) model, same branch, after reviewing a compiled equations-of-motion report**: two corrections. (1) Third-body `dU/dt` (previous entry) is removed — `third_body_dvdt`, its `z(10)` wiring, and `test_dvdt_3body.F` are all gone; third-body is now treated as static for the energy-element purpose. (2) `driver_KS.F`'s separate zonal code path (`aLegP`/`c_j`/`Re`/`ObyR`/`ZbyR`, computing `V_pot`/`q(j)`/`Tau_geo`) is removed entirely — `m=0` was never a distinct model, it's the general `(n,m)` formula's own `m=0` row (the identical three-term Legendre recursion). `geo_coeff_tess_general` now loads `m=0` rows too (previously skipped deliberately to avoid double-counting against the path that no longer exists), with `(2-delta_0m)`-aware normalization and the traditional `J_n=-C_n0` sign, verified to reproduce `geo_coeff_body`'s independent zonal loader exactly (`test_cunningham.F`, new CG3-CG5). **Decisive lossless-unification check**: a standalone program reimplementing the exact removed zonal formulas matched `tess_legendre_force`'s `m=0` output to relative `1e-16` on `V_pot`/`Tau_geo`/all four `q(j)`, with `dvdt_gen` confirmed exactly zero at `m=0` as expected. Real trade: the geopotential is now capped at `ntess_cap=72` uniformly, including zonal-only runs (was degree 2190 at O(n) cost for zonal-only before). 619/619 regression (627 − 9 removed + 1 new CG5 check), both `ifx` and `gfortran`, plus a clean end-to-end `driver_KS.exe` smoke run (`ngeo_deg=50`, 10 revs, no NaN). |
