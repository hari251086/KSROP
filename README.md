# KSROP — KS Regular Orbit Propagator

Orbit propagation using **Kustaanheimo–Stiefel (KS) regular elements** with a Runge–Kutta–Gill 4th-order numerical integrator for Earth-orbiting satellites.

**Author:** Harishkumar Sellamuthu · hari251086@gmail.com  
**Copyright:** 2018, Harishkumar Sellamuthu, All Rights Reserved

---

## Features

| Perturbation | Status | Notes |
|---|---|---|
| Earth oblateness (Jn) | Active | EGM2008 model, configurable up to degree 2190 |
| Luni-solar gravity | Active | Sun and Moon position vectors, configurable degree |
| Atmospheric drag | Active | Oblate, co-rotating exponential atmosphere referenced to perigee conditions; tabulated density (ATM.DAT) |
| Solar radiation pressure | Planned | Cannon-ball model (not yet implemented) |

---

## Files

### Source

| File | Description |
|---|---|
| `driver_KS.F` | Main program — reads OPM input, initialises KS elements, runs integration loop |
| `Subrouts.F` | Subroutines — coordinate transforms, integrator, OPM I/O, force models, utilities |
| `Legendre.F` | Zonal Legendre polynomial evaluation (`aLegP`, `p_polynomial_value`) |

### Input

| File | Description |
|---|---|
| `input.opm` | **CCSDS OPM v2.0** — initial state: epoch, position, velocity |
| `input.DAT` | Simulation parameters: revolutions, steps, force flags, drag coefficients |
| `const_new.dat` | Physical constants and geopotential degree settings |
| `ATM.DAT` | Tabulated atmospheric density and scale height (60–630 km, 291 entries) |
| `EGM2008_to2190_TideFree` | EGM2008 geopotential coefficients (~231 MB, required when `ngeo_deg ≥ 2`) |

### Output (generated at runtime)

| File | Format | Contents |
|---|---|---|
| `KSROP_YYYYMMDDTHHMMSS.oem` | **CCSDS OEM v2.0** | State trajectory: epoch + X Y Z Xdot Ydot Zdot at every integration step |
| `ksrop.opm` | **CCSDS OPM v2.0** | Initial-epoch state vector and Keplerian elements |
| `regular.out` | Internal | KS regular elements (debug) |

### Tests and tools

| File | Description |
|---|---|
| `test_subrouts.F` | Fortran unit tests — 20 tests covering core subroutines |
| `test_driver.py` | Python integration test — 10 checks on two-body propagation |
| `benchmark.py` | Performance profiler — timing across force model configurations |
| `Makefile` | Unix/Linux build targets: `all`, `tests`, `run`, `test`, `clean` |
| `build.bat` | Windows build script for Intel Fortran |

---

## Building

Requires Intel Fortran (`ifort`) or GNU Fortran (`gfortran`).

### Windows (Intel Fortran)

```bat
build.bat           :: build driver_KS.exe
build.bat tests     :: build test_subrouts.exe
build.bat test      :: build + run all tests
build.bat clean     :: remove build artefacts
```

### Unix / Linux / macOS

```bash
make            # build driver_KS
make tests      # build test_subrouts
make test       # build + run all tests (unit + integration)
make clean      # remove build artefacts
```

### Manual

```bash
# Propagator
ifort driver_KS.F Subrouts.F Legendre.F -o driver_KS.exe

# Unit tests
ifort test_subrouts.F Subrouts.F Legendre.F -o test_subrouts.exe
```

---

## Running

```bash
./driver_KS.exe        # Windows: driver_KS.exe
python test_driver.py driver_KS.exe   # integration test
python benchmark.py   driver_KS.exe   # performance profile
```

---

## Input Files

### `input.opm` — Initial state (CCSDS OPM v2.0)

Contains the orbital initial conditions. The driver reads `EPOCH`, `X`, `Y`, `Z`, `X_DOT`, `Y_DOT`, `Z_DOT` from the `STATE_VECTOR` block; all other keywords are ignored.

```
CCSDS_OPM_VERS = 2.0
CREATION_DATE  = 2016-09-20T00:00:00.000
ORIGINATOR     = KSROP

META_START
OBJECT_NAME    = SATELLITE
OBJECT_ID      = UNKNOWN
CENTER_NAME    = EARTH
REF_FRAME      = EME2000
TIME_SYSTEM    = UTC
META_STOP

STATE_VECTOR
EPOCH          = 2016-09-20T00:00:00.000
X              =        0.000000 [km]
Y              =    -5888.972700 [km]
Z              =    -3400.000000 [km]
X_DOT          =        9.500000 [km/s]
Y_DOT          =        0.000000 [km/s]
Z_DOT          =        0.000000 [km/s]
```

> The output `ksrop.opm` uses the same format, so it can be fed directly back as `input.opm` for a subsequent propagation (chained runs).

### `input.DAT` — Simulation parameters

Three lines only (orbital initial conditions have moved to `input.opm`):

```
nrev  istep  tole          ! Revolutions, steps/revolution, integrator tolerance
n_geo  n_sun  n_moon       ! Force model flags (non-zero = on, 0 = off)
BN  IDRAG  WE_rot  EPS_f  FR_rot   ! Drag parameters
```

**Example:**

```
10  360  1d-15
10  0  0
50.0  0  7.2921150d-5  3.35281066d-3  1.0
```

**Drag parameter line:**

| Parameter | Example | Description |
|---|---|---|
| `BN` | `50.0` | Ballistic coefficient CdA/m (kg/m²) |
| `IDRAG` | `1` | Drag: 1 = on, 0 = off |
| `WE_rot` | `7.2921150d-5` | Earth rotation rate (rad/s) |
| `EPS_f` | `3.35281066d-3` | Earth flattening (1/298.257) |
| `FR_rot` | `1.0` | Atmospheric co-rotation factor |

### `const_new.dat` — Physical constants

```
mu  R_Earth  AU  mu_Sun  mu_Moon
ngeo_deg  nsun_deg  nmoon_deg
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

### `ATM.DAT` — Atmosphere table

Two-block table in `8F10.3` / `8E10.3` format:

- **Block 1:** Density scale heights H (km) for 291 altitude levels
- **Block 2:** Atmospheric densities ρ for the same 291 levels

Altitude grid: 60–200 km in 1 km steps, 200–630 km in 2 km steps. Drag is suppressed automatically above 500 km.

---

## Output Files

### `KSROP_YYYYMMDDTHHMMSS.oem` — State trajectory (CCSDS OEM v2.0)

The filename encodes the **current UTC wall-clock time** at which the file is written (e.g. `KSROP_20260606T193024.oem`). Each run produces a uniquely named file; no previous output is overwritten.

```
CCSDS_OEM_VERS = 2.0
CREATION_DATE  = 2026-06-06T19:30:24.000   ! UTC wall-clock time of writing
ORIGINATOR     = KSROP

META_START
OBJECT_NAME    = SATELLITE
CENTER_NAME    = EARTH
REF_FRAME      = EME2000
TIME_SYSTEM    = UTC
START_TIME     = 2016-09-20T00:00:00.000   ! orbital epoch of first data point
STOP_TIME      = 2016-09-20T...            ! orbital epoch of last data point
META_STOP

DATA_START
2016-09-20T00:00:00.000        0.000000   -5888.972700   -3400.000000     9.500000000     0.000000000     0.000000000
...
DATA_STOP
```

**Time fields:**

| Field | Source |
|---|---|
| `CREATION_DATE` | Current UTC wall-clock (`utc_now_epoch`) |
| `START_TIME` | Orbital epoch of the first buffered state |
| `STOP_TIME` | Orbital epoch of the last buffered state (exact, not estimated) |
| Data line epochs | Orbital epoch at each completed integration step |

**`nrev × istep + 1`** data lines — one per completed RKG step plus the initial state — *unless the run stops early* (see **Early Termination** below), in which case the OEM is truncated at the last useful state. The entire trajectory is held in memory during the run and written atomically after propagation completes (or stops). Position in km (F16.6), velocity in km/s (F14.9).

For the production config (10 rev × 360 steps/rev) this produces **3,601 data lines**.

### Early Termination — re-entry and divergence detection

After every integration step, the propagator checks the freshly-computed altitude (`h_alt = R(1) - R_Earth`) and stops early — truncating the OEM and skipping the remaining steps/revolutions — under either of these conditions:

| Condition | Trigger | Final OEM entry | Console message |
|---|---|---|---|
| **Re-entry** | `h_alt < 80 km` | The just-computed (finite) low-altitude state | `Re-entry occurred: altitude = <h_alt> km (< 80 km) at epoch <epoch>` |
| **Divergence** | `h_alt` is non-finite (NaN), via the `x .ne. x` test | The last *valid* (finite) state — the NaN point is discarded | `Integration diverged (non-finite state); last valid epoch was <epoch>` |

Both cases print a follow-up line (`Propagation stopped; ...`) and proceed straight to writing the OEM, so the file always ends with a finite, physically meaningful state — never with NaN garbage. Divergence typically occurs during a steep orbital decay where the fixed per-revolution eccentric-anomaly stepping cannot resolve an unresolvably fast low-perigee pass.

### `ksrop.opm` — Initial Keplerian elements (CCSDS OPM v2.0)

Written once at the start of each run with the initial epoch state vector and Keplerian orbital elements (a, e, i, RAAN, AOP, M).

---

## Force Model Flags

`n_force` in `input.DAT` and `ngeo_deg` in `const_new.dat` control which perturbations are active:

| Flag | Effect |
|---|---|
| `n_force(1) = 0` | Geopotential off (point mass) |
| `n_force(1) ≠ 0` | Geopotential on — degree set by `ngeo_deg` in `const_new.dat` |
| `n_force(2) = 0` | Solar gravity off |
| `n_force(2) ≠ 0` | Solar gravity on — degree set by `nsun_deg` |
| `n_force(3) = 0` | Lunar gravity off |
| `n_force(3) ≠ 0` | Lunar gravity on — degree set by `nmoon_deg` |
| `IDRAG = 0` | Atmospheric drag off |
| `IDRAG = 1` | Atmospheric drag on |

---

## Method

The propagator uses **KS regularisation** to remove the 1/r singularity. The 3D equations of motion are lifted to a 4D harmonic oscillator in KS space and integrated with **Runge–Kutta–Gill 4th order** using the generalised eccentric anomaly E as the independent variable.

Perturbing forces are transformed to KS space via the **L(u)** matrix and added to the KS element ODEs:

| Element | Rate equation |
|---|---|
| `z(j+1)`, `z(j+5)` | State elements — conservative (geo + third-body) + drag |
| `z(1)` | Time element — geopotential + third-body + drag contributions |
| `z(10)` | Energy element — drag dissipation: dω/ds = −½(u̇·q_drag)/ω |

The step size is scaled by the frequency ratio Γ = ω/ω_Kep to maintain accuracy across eccentricities.

---

## Performance

Benchmarked on the production orbit (a ≈ 14,770 km, e ≈ 0.8):

| Configuration | Throughput | Wall time (10 rev, 360 steps/rev) |
|---|---|---|
| Two-body (no perturbations) | ~140,000 steps/s | ~31 ms |
| J2 geopotential (ngeo=2) | ~62,000 steps/s | ~57 ms |
| Degree-50 geopotential (ngeo=50) | ~55,000 steps/s | ~65 ms |
| With atmospheric drag (IDRAG=1) | ~108,000 steps/s | ~34 ms |

**EGM2008 file read:** Streaming `geo_coeff` reads only O(n²) lines for degree n — J2 reads 3 lines instead of 2,401,333. File-read cost for any degree is now <1 ms (was 2,600 ms before optimisation).

---

## Tests

### Unit tests — `test_subrouts.F`

```bash
ifort test_subrouts.F Subrouts.F Legendre.F -o test_subrouts.exe
./test_subrouts.exe
```

20 tests: `dotp3`, `dotp4`, `vmn`, `cross`, `INTPOL` (×2), `aLegP` P2/P3/P4, `car2ks`→`ks2car` roundtrip (×6), `car2oe` a/e/i.

### Integration test — `test_driver.py`

```bash
python test_driver.py driver_KS.exe
```

10 checks on a 1-revolution two-body propagation: OEM row count (361 = 1×360+1), DATA_START/DATA_STOP structure, energy conservation, angular momentum conservation, orbit closure (position and velocity), semi-major axis vs vis-viva, semi-major axis conserved, OPM output present.

### Performance profiling — `benchmark.py`

```bash
python benchmark.py driver_KS.exe
```

Times integrator, step-size sensitivity, drag overhead, and EGM2008 file-read cost across multiple configurations.

---

## Subroutines Reference

### Coordinate transforms

| Subroutine | Description |
|---|---|
| `car2ks(x,xd,u,us,w)` | Cartesian → KS |
| `ks2car(u,us,x,xd,w)` | KS → Cartesian |
| `ks2ksr(y,u,us,E,cse,sie)` | KS regular elements → u, us |
| `car2oe(x,xd,pek)` | Cartesian → Keplerian elements |
| `oe2car(pek,x,xd,tol)` | Keplerian → Cartesian (Kepler equation solver) |
| `u2uu(u,uu)` | KS index swap |
| `u2qu(u,qu,j)` | KS rearrangement for luni-solar terms |

### Force models and coefficients

| Subroutine | Description |
|---|---|
| `geo_coeff(n,c_j)` | Stream EGM2008 zonal harmonics up to degree n (reads O(n²) lines) |
| `force_models(n_for,ngeo,s,m)` | Apply force model on/off flags |
| `INTPOL(XT,YT,M1,X1,Y1)` | Linear interpolation in sorted atmosphere table |

### Ephemeris

| Subroutine | Description |
|---|---|
| `solarnpv(dj,s)` | Sun position vector (geocentric inertial, km) |
| `lunarpv(dj,tm)` | Moon position vector (geocentric inertial, km) |
| `aLegP(n,x,P)` | Zonal Legendre polynomials (computed to degree 49) |

### I/O and time

| Subroutine | Description |
|---|---|
| `read_opm(iunit,x,xd,cal)` | CCSDS OPM v2.0 parser — extracts epoch and state vector |
| `write_opm(iunit,epochstr,x,xd,pek)` | CCSDS OPM v2.0 writer — emits header, state vector and Keplerian elements |
| `read_oem(iunit,maxpts,traj_jd,traj_x,traj_xd,npts)` | CCSDS OEM v2.0 parser — reads the `DATA_START`/`DATA_STOP` ephemeris block into Julian-date/position/velocity buffers |
| `write_oem(iunit,creation_str,start_str,stop_str,traj_jd,traj_x,traj_xd,npts)` | CCSDS OEM v2.0 writer — emits header and `DATA_START`/`DATA_STOP` ephemeris block |
| `read_cdm(iunit,tca_cal,miss_dist,pc,x1,xd1,cov1,x2,xd2,cov2)` | CCSDS CDM v1.0 parser — extracts TCA, miss distance, collision probability, and each OBJECT1/OBJECT2 state vector + lower-triangular RTN covariance (expanded to a symmetric 6×6, ordered R,T,N,RDOT,TDOT,NDOT); raw numeric values, no unit conversion (cf. `read_opm`) |
| `write_cdm(iunit,creation_str,tca_str,miss_dist,rel_speed,pc,name1,x1,xd1,cov1,name2,x2,xd2,cov2)` | CCSDS CDM v1.0 writer — emits header/relative-geometry block and one OBJECTn block per object (state vector + RTN covariance) via `write_cdm_object` |
| `cdm_object_field`, `cdm_cov_index`, `cdm_rtn_index` | Internal helpers for `read_cdm`/`write_cdm`: map OBJECTn keyword/value pairs and CDM covariance keywords (`CR_R` … `CNDOT_NDOT`) to/from the symmetric 6×6 RTN matrix |
| `parse_epoch(estr,cal)` | Parse `YYYY-MM-DDTHH:MM:SS.sss` → cal(6) |
| `jd2epoch(djd,epochstr)` | Julian date → CCSDS epoch string |
| `utc_now_epoch(epochstr,compact)` | Current UTC wall-clock → CCSDS epoch string and compact filename token |
| `cal2jd(cal,djulian)` | Calendar date → Julian date |

### Vector utilities

| Function | Description |
|---|---|
| `dotp3(x,y)` | 3-vector dot product |
| `dotp4(x,y)` | 4-vector dot product |
| `vmn(x)` | 3-vector magnitude |
| `cross(x,y,z)` | 3-vector cross product |

### Integrator

| Function | Description |
|---|---|
| `rkgil(n,y,f,x,h,nt)` | Runge–Kutta–Gill 4th-order step (4-stage, fixed step) |

---

## Known Issues

- `car2oe` may produce NaN for near-circular or near-equatorial orbits (special cases partially handled).
- `aLegP` is internally hardcoded to evaluate polynomials up to degree 49; results for `ngeo_deg ≥ 50` are incorrect beyond that degree.
- EGM2008 file (~231 MB) is not included in the repository; set `ngeo_deg = 0` in `const_new.dat` to run without it (point-mass gravity).
- Solar radiation pressure is not yet implemented.

---

## Revision History

| Date | Change |
|---|---|
| 2018-06-15 | Initial program, J2 only |
| 2018-06-16 | Revolution-by-revolution output; nth-degree Legendre polynomial |
| 2018-09-13 | nth-degree geopotential up to 2190 (EGM2008 Jn) |
| 2021-07-21 | Legendre polynomial subroutine (`aLegP`) added |
| 2026-06-06 | Luni-solar third-body perturbations |
| 2026-06-06 | Atmospheric drag (oblate exponential model, ATM.DAT); fixed time-element `Tau_term` bug |
| 2026-06-06 | Unit tests (20), integration test (10), benchmark script |
| 2026-06-06 | Performance optimisation: streaming `geo_coeff` (2600 ms → <1 ms for J2); aLegP guards; removed 115 MB static array; 3× integrator speedup |
| 2026-06-06 | CCSDS OEM v2.0 output; CCSDS OPM v2.0 output (`ksrop.opm`) |
| 2026-06-06 | CCSDS OPM v2.0 input (`input.opm`); `input.DAT` reduced to simulation parameters only |
| 2026-06-06 | OEM trajectory buffered in memory; written atomically after run; `STOP_TIME` is exact final epoch |
| 2026-06-06 | OEM output at every integration step (`nrev × istep + 1` data points) |
| 2026-06-06 | OEM filename `KSROP_YYYYMMDDTHHMMSS.oem` and `CREATION_DATE` use current UTC wall-clock |
| 2026-06-07 | Atmospheric drag model replaced with an oblate, co-rotating exponential atmosphere referenced to perigee conditions (ported from `KSJLSDNP2.F`'s physical model and logic) |
| 2026-06-07 | Early-termination check added: propagation stops and truncates the OEM at the last useful state on re-entry (`altitude < 80 km`) or numerical divergence (non-finite state, `x .ne. x` NaN test) |
| 2026-06-07 | Code-structure cleanup: removed dead subroutines `car2sph`, `sph2ks`, `car2ksnew` (~165 lines) and consolidated the 4×-duplicated sun/moon auxiliary computation into a shared `third_body_aux` subroutine, fixing a latent first-step sign-convention inconsistency in the luni-solar perturbation terms along the way (~300 fewer lines overall) |
| 2026-06-07 | Fixed `cal2jd`: the time-of-day fields were read from the wrong `cal` indices, offsetting every output epoch/timestamp by ~84 days from the true input epoch (propagation itself — which runs on elapsed seconds — was unaffected; only the date *labels* were wrong) |
| 2026-06-07 | Fixed a UTF-8 byte-order-mark in the tracked `input/input.DAT` that made `driver_KS.exe` crash immediately (`list-directed I/O syntax error`) on a fresh checkout |
| 2026-06-07 | Moved CCSDS OPM/OEM file I/O out of `driver_KS.F` and into `Subrouts.F` as reusable subroutines: added `write_opm`, `read_oem`, `write_oem` alongside the existing `read_opm`; `driver_KS.F` now calls these instead of writing the records inline |
| 2026-06-07 | Added CCSDS CDM v1.0 (Conjunction Data Message) I/O to `Subrouts.F`: `read_cdm`/`write_cdm` plus internal helpers `cdm_object_field`, `cdm_cov_index`, `cdm_rtn_index`, parsing/emitting the relative-geometry summary and each OBJECT1/OBJECT2 state vector + lower-triangular RTN covariance; `test_subrouts.F` gained a read/roundtrip test against the public CCSDS 508.0-B-1 sample CDM (`input/cdm_sample.kvn`) |
