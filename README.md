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
| Atmospheric drag | Active | Oblate exponential model, tabulated density (ATM.DAT) |
| Solar radiation pressure | Planned | Cannon-ball model (not yet implemented) |

---

## Files

| File | Description |
|---|---|
| `driver_KS.F` | Main program — reads inputs, initialises KS elements, runs integration loop |
| `Subrouts.F` | Subroutines — coordinate transforms, integrator, force models, utility functions |
| `Legendre.F` | Associated Legendre polynomial evaluation (`aLegP`) |
| `const_new.dat` | Physical constants and geopotential degree settings |
| `input.DAT` | Initial conditions and simulation parameters |
| `ATM.DAT` | Tabulated atmospheric density and scale height (60–630 km, 291 entries) |
| `EGM2008_to2190_TideFree` | EGM2008 geopotential coefficients (required when `ngeo_deg > 1`) |

### Output files (generated at runtime)

| File | Contents |
|---|---|
| `state.out` | Cartesian state vector (position km, velocity km/s) per revolution |
| `regular.out` | KS regular elements |
| `kepler.out` | Keplerian orbital elements |

---

## Building

Requires an Intel Fortran compiler (`ifort`) or compatible (e.g. `gfortran`).

```bash
# Intel Fortran
ifort driver_KS.F Subrouts.F Legendre.F -o ksrop

# GNU Fortran
gfortran driver_KS.F Subrouts.F Legendre.F -o ksrop
```

Run:

```bash
./ksrop
```

---

## Input Files

### `const_new.dat`

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
| `ngeo_deg` | 2–2190 | Geopotential degree (0 or 1 = point mass) |
| `nsun_deg` | 0–2190 | Solar gravity degree |
| `nmoon_deg` | 0–2190 | Lunar gravity degree |

### `input.DAT`

Six lines, read in order:

```
x1  x2  x3                              ! Initial position (km, geocentric inertial)
xd1  xd2  xd3                           ! Initial velocity (km/s, geocentric inertial)
nrev  istep  tole                        ! Revolutions, steps/rev, integrator tolerance
yyyy  mm  dd  hh  mm  ss                 ! Epoch (calendar date and time)
n_geo  n_sun  n_moon                     ! Force model flags (1=on, 0=off)
BN  IDRAG  WE_rot  EPS_f  FR_rot        ! Drag parameters (see below)
```

**Example (`input.DAT`):**

```
0.0  -5888.9727  -3400.0
9.5  0.0  0.0
1  360  1d-15
2016  09  20  0  0  0
2  1  1
50.0  1  7.2921150d-5  3.35281066d-3  1.0
```

**Drag parameter line:**

| Parameter | Example | Description |
|---|---|---|
| `BN` | `50.0` | Ballistic coefficient CdA/m (kg/m²) |
| `IDRAG` | `1` | Drag switch: 1 = on, 0 = off |
| `WE_rot` | `7.2921150d-5` | Earth rotation rate (rad/s) |
| `EPS_f` | `3.35281066d-3` | Earth flattening (1/298.257) |
| `FR_rot` | `1.0` | Atmospheric co-rotation factor |

### `ATM.DAT`

Two-block table read in `8F10.3` / `8E10.3` format:

- **Block 1:** Density scale heights H (km) for 291 altitude levels  
- **Block 2:** Atmospheric densities ρ (kg/m³ after ×10¹⁰ scaling) for the same 291 levels  

Altitude grid: 60–200 km in 1 km steps, 200–630 km in 2 km steps.  
Drag is suppressed automatically when the satellite altitude exceeds 500 km.

---

## Force Model Flags

The three values on line 5 of `input.DAT` (and `ngeo_deg` in `const_new.dat`) control which perturbations are active:

```
n_force(1)  — Geopotential (0 = point mass, 1 = active; degree set by ngeo_deg)
n_force(2)  — Solar gravity (0 = off, 1 = on; degree set by nsun_deg)
n_force(3)  — Lunar gravity (0 = off, 1 = on; degree set by nmoon_deg)
IDRAG       — Atmospheric drag (0 = off, 1 = on; line 6 of input.DAT)
```

---

## Method

The propagator uses the **KS regularisation** to remove the 1/r singularity of the two-body problem. The 3D equations of motion are lifted to a 4D harmonic oscillator in KS space, then integrated with the **Runge–Kutta–Gill** scheme using the generalised eccentric anomaly E as the independent variable.

Perturbing forces are projected into KS space via the **L(u)** matrix transformation and added to the KS element equations:

- **z(j+1), z(j+5)** — KS state element rates (conservative + drag force)
- **z(1)** — Time element rate (geopotential + third-body + drag contributions)
- **z(10)** — Energy element rate (drag dissipation: dω/ds = −½(u̇·qdrag)/ω)

The step size scales with the frequency ratio Γ = ω/ω_Kep to maintain accuracy across different eccentricities.

---

## Subroutines Reference

| Subroutine / Function | Description |
|---|---|
| `car2ks(x,xd,u,us,w)` | Cartesian → KS coordinates |
| `ks2car(u,us,x,xd,w)` | KS → Cartesian coordinates |
| `ks2ksr(y,u,us,E,cse,sie)` | KS regular elements → u, us |
| `car2oe(x,xd,pek)` | Cartesian → Keplerian orbital elements |
| `oe2car(pek,x,xd)` | Keplerian → Cartesian (solves Kepler equation) |
| `geo_coeff(n,c_j)` | Read EGM2008 Jn coefficients from file |
| `aLegP(n,x,P)` | Zonal Legendre polynomials P2…Pn |
| `solarnpv(dj,s)` | Sun position vector (geocentric inertial, km) |
| `lunarpv(dj,tm)` | Moon position vector (geocentric inertial, km) |
| `INTPOL(XT,YT,M1,X1,Y1)` | Linear interpolation in sorted table |
| `cal2jd(cal,djulian)` | Calendar date → Julian date |
| `force_models(n_for,ngeo,s,m)` | Apply force model on/off flags |
| `rkgil(n,y,f,x,h,nt)` | Runge–Kutta–Gill 4th-order integrator step |
| `u2uu(u,uu)` | KS variable index swap (u→uu) |
| `u2qu(u,qu,j)` | KS variable rearrangement for luni-solar terms |
| `dotp3(x,y)` | 3-vector dot product |
| `dotp4(x,y)` | 4-vector dot product |
| `vmn(x)` | 3-vector magnitude |
| `cross(x,y,z)` | 3-vector cross product |

---

## Known Issues

- `car2oe` may produce NaN outputs for near-circular or near-equatorial orbits (special cases partially handled but not fully validated).
- EGM2008 coefficient file (`EGM2008_to2190_TideFree`) is large (~500 MB) and not included in the repository; set `ngeo_deg = 2` to use J2 only without the file.
- Solar radiation pressure is not yet implemented.

---

## Revision History

| Date | Change |
|---|---|
| 2018-06-15 | Initial program, J2 only |
| 2018-06-16 | Revolution-by-revolution trajectory dump; nth-degree Legendre polynomial |
| 2018-09-13 | nth-degree geopotential up to 2190×0 (EGM2008 Jn) |
| 2025 | Luni-solar third-body perturbations |
| 2025 | Atmospheric drag (oblate exponential model, ATM.DAT); fixed time-element bug |
