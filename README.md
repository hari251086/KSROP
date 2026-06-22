# KSROP — KS Regular Orbit Propagator

Orbit propagation using **Kustaanheimo–Stiefel (KS) regular elements** with a Runge–Kutta–Gill 4th-order numerical integrator for Earth-orbiting satellites.

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
├── driver_KS.F                          Main propagator program
├── Subrouts.F                           Shared subroutines (transforms, I/O, force models)
├── Legendre.F                           Zonal Legendre polynomial evaluation
├── TLEread.F                            TLE reader + SGP4 + TEME→J2000 conversion
├── tle2opm.F                            TLE-to-OPM converter tool
│
├── input/
│   ├── const_new.dat                    Physical constants (single source of truth)
│   ├── input.opm                        Initial state for driver_KS (CCSDS OPM v2.0)
│   ├── input.DAT                        Simulation parameters
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
├── test_subrouts.F                      Unit tests (57 checks)
├── test_tle.F                           TLE reader tests (147 checks)
├── test_tle2sv.F                        SGP4/SDP4/frame tests (156 checks)
├── test_tle2opm.F                       TLE-to-OPM pipeline tests (21 checks)
├── test_bugs.F                          Bug regression tests (17 checks)
├── test_driver.py                       Integration test (10 checks)
├── test_initial_conditions.py           Multi-orbit test (110 checks)
├── benchmark.py                         Performance profiler
├── lint_check.sh                        Automated F77 lint checks (line length, precision, common blocks)
├── test_all.sh                          Unified test runner (lint + all test suites)
├── Makefile                             Unix build
└── build.bat                            Windows build (Intel Fortran)
```

---

## 3. Building

Requires **Intel oneAPI Fortran** (`ifx`) or **GNU Fortran** (`gfortran`), plus a C/C++ linker (MSVC on Windows).

### Windows — Intel oneAPI ifx 2025.0

The Intel compiler requires both the Visual Studio and Intel environments to be initialised in the same shell session:

```bat
:: Single-command build (from cmd.exe or PowerShell via cmd /c)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"
cd /d "C:\Users\hari2\OneDrive\Documents\GitHub\KSROP"

:: Propagator
ifx driver_KS.F Subrouts.F Legendre.F TLEread.F /exe:driver_KS.exe

:: TLE-to-OPM converter
ifx tle2opm.F Subrouts.F TLEread.F Legendre.F /exe:tle2opm.exe

:: Tests
ifx test_subrouts.F Subrouts.F Legendre.F /exe:test_subrouts.exe
ifx test_bugs.F Subrouts.F Legendre.F /exe:test_bugs.exe
ifx test_tle.F Subrouts.F TLEread.F Legendre.F /exe:test_tle.exe
ifx test_tle2sv.F Subrouts.F TLEread.F Legendre.F /exe:test_tle2sv.exe
ifx test_tle2opm.F Subrouts.F TLEread.F Legendre.F /exe:test_tle2opm.exe
```

> **Note:** Use `/exe:name.exe` (not `-o`) for the output flag with `ifx` on Windows.

### Unix / Linux / macOS

```bash
make            # build driver_KS
make tests      # build all test executables
make test       # build + run all tests
make clean      # remove build artefacts
```

### Manual (gfortran)

```bash
gfortran driver_KS.F Subrouts.F Legendre.F TLEread.F -o driver_KS.exe
gfortran tle2opm.F Subrouts.F TLEread.F Legendre.F -o tle2opm.exe
gfortran test_subrouts.F Subrouts.F Legendre.F -o test_subrouts.exe
gfortran test_bugs.F Subrouts.F Legendre.F -o test_bugs.exe
gfortran test_tle.F Subrouts.F TLEread.F Legendre.F -o test_tle.exe
gfortran test_tle2sv.F Subrouts.F TLEread.F Legendre.F -o test_tle2sv.exe
gfortran test_tle2opm.F Subrouts.F TLEread.F Legendre.F -o test_tle2opm.exe
```

---

## 4. Running

### Propagation

```bash
./driver_KS.exe
```

Reads `input/input.opm` (initial state), `input/input.DAT` (parameters), `input/const_new.dat` (constants). Writes trajectory to `output/KSROP_*.oem`.

### TLE-to-OPM conversion

```bash
./tle2opm.exe
```

Reads `input/tle2opm.cfg`, selects closest TLE entry by NORAD/epoch, converts via SGP4 to J2000, writes `input/KSROP_<EPOCH>_<NORAD>.opm` and `input/input.opm`.

### Tests

```bash
./test_subrouts.exe                          # 47 unit tests
./test_bugs.exe                              # 17 regression tests
./test_tle.exe                               # 147 TLE reader tests
./test_tle2sv.exe                            # 156 SGP4/SDP4/frame tests
./test_tle2opm.exe                           # 21 pipeline tests
python test_driver.py driver_KS.exe          # 10 integration tests
python test_initial_conditions.py driver_KS.exe  # 110 multi-orbit tests

# Lint + all tests in one command
bash test_all.sh                             # lint + 388 Fortran + integration
bash lint_check.sh                           # lint only (line length, precision, common blocks)
```

**Total: 508 automated checks + 5 lint rules.**

---

## 5. Input Files

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

### `input.DAT` — Simulation parameters

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

Two-block table: scale heights and densities for 291 altitude levels (60–630 km). Drag is suppressed above 500 km.

---

## 6. Output Files

### `KSROP_YYYYMMDDTHHMMSS.oem` — Trajectory (CCSDS OEM v2.0)

Uniquely named per run (UTC wall-clock timestamp). Contains `nrev × istep + 1` state vectors (position km, velocity km/s) at every integration step.

**Early termination:** Propagation stops on re-entry (alt < 80 km) or divergence (NaN). The OEM is truncated at the last valid state.

### `ksrop.opm` — Initial Keplerian elements

Written once at run start with the epoch state vector and osculating elements (a, e, i, Ω, ω, M).

---

## 7. Force Model Configuration

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

## 8. Method

The propagator uses **KS regularisation** to remove the 1/r singularity. The 3D equations of motion are lifted to a 4D harmonic oscillator in KS space, integrated with **Runge–Kutta–Gill 4th order** using the generalised eccentric anomaly E as the independent variable.

| KS Element | Rate equation |
|---|---|
| `z(j+1)`, `z(j+5)` | State — conservative (geo + third-body) + drag + SRP |
| `z(1)` | Time — geopotential + third-body + drag + SRP |
| `z(10)` | Energy — non-conservative: dω/ds = −½(u̇·(q_drag+q_srp))/ω |

Step size is scaled by Γ = ω/ω_Kep to maintain accuracy across eccentricities.

---

## 9. Subroutines Reference

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

## 10. Performance

| Configuration | Throughput | Wall time (10 rev × 360 steps) |
|---|---|---|
| Two-body | ~140,000 steps/s | ~31 ms |
| J2 geopotential | ~62,000 steps/s | ~57 ms |
| Degree-50 geopotential | ~55,000 steps/s | ~65 ms |
| With drag | ~108,000 steps/s | ~34 ms |

EGM2008 streaming read: O(n²) lines for degree n (J2 = 3 lines, not 2.4M).

---

## 11. Known Issues

- EGM2008 file (~231 MB) not included; set `ngeo_deg = 0` for point-mass gravity.
- `aLegP` hardcoded to degree 49; results incorrect beyond that.
- SRP is cannonball only; no tesseral harmonics or geometry-dependent variations.

---

## 12. Revision History

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
