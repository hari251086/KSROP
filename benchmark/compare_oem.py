"""compare_oem.py -- compare two CCSDS OEM files and report 3D position error.

Usage:
    python compare_oem.py <ref.oem> <cmp.oem> [--plot] [--csv out.csv]

The reference OEM (ref) is used as the base time grid.  The comparison OEM
(cmp) is linearly interpolated onto that grid.  Only the overlapping time
span is evaluated.

Output:
    max 3D error [km and m], RMS error [m], time of max error [s from epoch]
    Optional: time-history written to CSV, optional matplotlib plot.
"""

import argparse
import math
import sys
from datetime import datetime, timezone


def parse_epoch(s):
    s = s.strip()
    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%jT%H:%M:%S.%f',
                '%Y-%m-%dT%H:%M:%S',   '%Y-%jT%H:%M:%S'):
        try:
            return datetime.strptime(s, fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            pass
    raise ValueError('Unrecognised epoch format: ' + s)


def epoch_to_sec(dt, ref_dt):
    return (dt - ref_dt).total_seconds()


def read_oem(path):
    """Return list of (t_sec_from_first, x, y, z, vx, vy, vz) tuples."""
    records = []
    in_data = False
    t0 = None
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.strip() == 'DATA_START':
                in_data = True
                continue
            if line.strip() == 'DATA_STOP':
                in_data = False
                continue
            if not in_data:
                continue
            line = line.strip()
            if not line:
                continue
            # First 23 chars = epoch, rest = 6 floats
            if len(line) < 24:
                continue
            epoch_str = line[:23]
            try:
                dt = parse_epoch(epoch_str)
            except ValueError:
                continue
            try:
                vals = list(map(float, line[23:].split()))
            except ValueError:
                continue
            if len(vals) < 6:
                continue
            if t0 is None:
                t0 = dt
            t_sec = epoch_to_sec(dt, t0)
            records.append((t_sec, vals[0], vals[1], vals[2],
                            vals[3], vals[4], vals[5]))
    return records, t0


def interp1(ts, xs, t):
    """Linear interpolation of scalar array xs at t; clamp to endpoints."""
    if t <= ts[0]:
        return xs[0]
    if t >= ts[-1]:
        return xs[-1]
    # binary search
    lo, hi = 0, len(ts) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t:
            lo = mid
        else:
            hi = mid
    frac = (t - ts[lo]) / (ts[hi] - ts[lo])
    return xs[lo] + frac * (xs[hi] - xs[lo])


def main():
    ap = argparse.ArgumentParser(description='Compare two CCSDS OEM files')
    ap.add_argument('ref',  help='Reference OEM (driver_KS output)')
    ap.add_argument('cmp',  help='Comparison OEM (KSBENCH output)')
    ap.add_argument('--plot', action='store_true', help='Show matplotlib plot')
    ap.add_argument('--csv', metavar='FILE', help='Write error time-history to CSV')
    args = ap.parse_args()

    ref_rec, ref_t0 = read_oem(args.ref)
    cmp_rec, cmp_t0 = read_oem(args.cmp)

    if not ref_rec:
        sys.exit('No data records found in ' + args.ref)
    if not cmp_rec:
        sys.exit('No data records found in ' + args.cmp)

    # Re-align cmp times to ref_t0
    if ref_t0 != cmp_t0:
        offset = epoch_to_sec(cmp_t0, ref_t0)
        cmp_rec = [(t + offset, x, y, z, vx, vy, vz)
                   for t, x, y, z, vx, vy, vz in cmp_rec]

    # Build cmp arrays for interpolation
    cmp_t  = [r[0] for r in cmp_rec]
    cmp_x  = [r[1] for r in cmp_rec]
    cmp_y  = [r[2] for r in cmp_rec]
    cmp_z  = [r[3] for r in cmp_rec]

    # Overlap window
    t_start = max(ref_rec[0][0],  cmp_t[0])
    t_end   = min(ref_rec[-1][0], cmp_t[-1])
    if t_end <= t_start:
        sys.exit('No overlapping time span between the two OEM files.')

    # Evaluate error at each ref epoch within overlap
    err_t  = []
    err_dr = []
    sum2   = 0.0
    n      = 0
    for rec in ref_rec:
        t = rec[0]
        if t < t_start or t > t_end:
            continue
        xi = interp1(cmp_t, cmp_x, t)
        yi = interp1(cmp_t, cmp_y, t)
        zi = interp1(cmp_t, cmp_z, t)
        dr = math.sqrt((rec[1]-xi)**2 + (rec[2]-yi)**2 + (rec[3]-zi)**2)
        err_t.append(t)
        err_dr.append(dr)
        sum2 += dr*dr
        n    += 1

    if n == 0:
        sys.exit('No overlapping ref points after filtering.')

    max_err = max(err_dr)
    rms_err = math.sqrt(sum2 / n)
    t_max   = err_t[err_dr.index(max_err)]
    span    = t_end - t_start

    print('Reference OEM : ' + args.ref)
    print('Comparison OEM: ' + args.cmp)
    print('Epoch (ref t0): ' + str(ref_t0))
    print('Overlap span  : %.1f s (%.3f days)' % (span, span/86400.0))
    print('Points compared: %d' % n)
    print()
    print('Max 3D position error : %.6f km  (%.1f m)  at t = %.1f s' %
          (max_err, max_err*1000.0, t_max))
    print('RMS 3D position error : %.6f km  (%.1f m)' %
          (rms_err, rms_err*1000.0))

    if args.csv:
        with open(args.csv, 'w') as fh:
            fh.write('t_sec,dr_km,dr_m\n')
            for t, dr in zip(err_t, err_dr):
                fh.write('%.3f,%.9f,%.3f\n' % (t, dr, dr*1000.0))
        print('Error time-history written to: ' + args.csv)

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            t_day = [t/86400.0 for t in err_t]
            dr_m  = [dr*1000.0 for dr in err_dr]
            plt.figure(figsize=(10, 4))
            plt.plot(t_day, dr_m, lw=0.8)
            plt.xlabel('Time [days from epoch]')
            plt.ylabel('3D position error [m]')
            plt.title('KSROP vs KSBENCH: 3D position error')
            plt.grid(True)
            plt.tight_layout()
            plt.show()
        except ImportError:
            print('matplotlib not available; skipping plot.')


if __name__ == '__main__':
    main()
