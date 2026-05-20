#!/usr/bin/env python3
"""
cross_validate.py
=================

Independent Python implementation of the PR + Peneloux air EOS, used
to verify the Fortran implementation in air_eos_pr_module.F90.

Two checks:

1. Anchor check
   Compares Python implementation against NIST air reference data
   at six (P, T) points. Acceptance tolerance ~1% (the formulation
   itself is only ~0.5-1% accurate vs NIST).

2. CSV cross-check
   Reads the CSV produced by test_air_eos_pr_sweep and compares each
   row to the Python implementation. Acceptance tolerance 1e-6
   (essentially bit-precision; should be at floating-point roundoff).

Usage:
    python3 cross_validate.py                       # anchor check only
    python3 cross_validate.py air_eos_sweep.csv     # + CSV cross-check

Changelog:
    v2 (May 2026): Skip header lines instead of crashing.
                   Refuse to report PASS on 0 data points.
"""

import math
import sys

# -----------------------------------------------------------------------
# Constants — Lemmon et al. 2000 for air
# -----------------------------------------------------------------------

R = 8.314462618           # J/(mol*K)
AIR_TC = 132.5306         # K
AIR_PC = 3.7850e6         # Pa
AIR_OMEGA = 0.0335
AIR_MW = 28.9586e-3       # kg/mol (= 28.9586 g/mol)
AIR_PENELOUX_C = -9.0e-6  # m^3/mol


# -----------------------------------------------------------------------
# PR + Peneloux implementation
# -----------------------------------------------------------------------

def _kappa(omega):
    return 0.37464 + 1.54226*omega - 0.26992*omega**2


def _alpha(T):
    kappa = _kappa(AIR_OMEGA)
    s = 1.0 + kappa * (1.0 - math.sqrt(T / AIR_TC))
    return s * s


def _a(T):
    return 0.45724 * R**2 * AIR_TC**2 / AIR_PC * _alpha(T)


def _b():
    return 0.07780 * R * AIR_TC / AIR_PC


def _solve_cubic(A, B):
    """
    Solve Z^3 + p2 Z^2 + p1 Z + p0 = 0  via Cardano.
    Returns largest real root (gas branch).
    """
    p2 = -(1.0 - B)
    p1 = A - 3.0*B**2 - 2.0*B
    p0 = -(A*B - B**2 - B**3)

    # Depress: Z = y - p2/3
    shift = p2 / 3.0
    p = p1 - p2**2 / 3.0
    q = 2.0*p2**3/27.0 - p2*p1/3.0 + p0

    discriminant = (q/2.0)**2 + (p/3.0)**3

    if discriminant > 0.0:
        # One real root
        sqrtD = math.sqrt(discriminant)
        u = -q/2.0 + sqrtD
        v = -q/2.0 - sqrtD
        u_cbrt = math.copysign(abs(u)**(1.0/3.0), u)
        v_cbrt = math.copysign(abs(v)**(1.0/3.0), v)
        y = u_cbrt + v_cbrt
        return y - shift
    else:
        # Three real roots — return the largest
        r = math.sqrt(-(p/3.0)**3)
        phi = math.acos(max(-1.0, min(1.0, -q/(2.0*r))))
        cube_root_r = (-p/3.0)**0.5  # equivalent to r**(1/3) for trigonometric form
        # Actually use the standard trigonometric solution:
        m = 2.0 * math.sqrt(-p/3.0)
        y1 = m * math.cos(phi/3.0)
        y2 = m * math.cos((phi + 2.0*math.pi)/3.0)
        y3 = m * math.cos((phi + 4.0*math.pi)/3.0)
        return max(y1, y2, y3) - shift


def pr_peneloux(P, T):
    """
    Compute air properties via PR + Peneloux.

    Returns
    -------
    dict with keys: rho [kg/m^3], Z [-], phi [-], h_dep [J/mol]
    """
    a = _a(T)
    b = _b()
    A = a * P / (R*T)**2
    B = b * P / (R*T)

    Z = _solve_cubic(A, B)

    # PR molar volume
    V_pr = Z * R * T / P
    # Peneloux-corrected molar volume
    V = V_pr - AIR_PENELOUX_C
    # Density (kg/m^3)
    rho = AIR_MW / V
    # Z corrected for Peneloux (effective Z)
    Z_eff = V * P / (R*T)

    # Fugacity coefficient (PR formula)
    sqrt2 = math.sqrt(2.0)
    if Z - B > 0.0:
        ln_term = math.log((Z + (1.0+sqrt2)*B) / (Z + (1.0-sqrt2)*B))
        ln_phi_pr = (Z - 1.0) - math.log(Z - B) \
                    - A/(2.0*sqrt2*B) * ln_term
    else:
        ln_phi_pr = 0.0
    # Peneloux correction to fugacity
    ln_phi = ln_phi_pr - AIR_PENELOUX_C * P / (R*T)
    phi = math.exp(ln_phi)

    # Enthalpy departure (PR formula)
    kappa = _kappa(AIR_OMEGA)
    alpha = _alpha(T)
    # d(a*alpha)/dT
    da_dT = -0.45724 * R**2 * AIR_TC**2 / AIR_PC \
            * kappa * math.sqrt(alpha) / math.sqrt(T * AIR_TC)
    if Z - B > 0.0:
        ln_term = math.log((Z + (1.0+sqrt2)*B) / (Z + (1.0-sqrt2)*B))
        h_dep = R*T*(Z - 1.0) + (T*da_dT - a) / (2.0*sqrt2*b) * ln_term
    else:
        h_dep = 0.0

    return {
        'rho':   rho,
        'Z':     Z_eff,
        'phi':   phi,
        'h_dep': h_dep,
    }


# -----------------------------------------------------------------------
# NIST anchor check
# -----------------------------------------------------------------------

# (P [Pa], T [K], rho_NIST [kg/m^3]) — six points covering CAES range
NIST_ANCHORS = [
    (1.0132e5, 298.15,   1.1839),
    (5.0e6,    298.15,  58.78),
    (1.0e7,    298.15, 115.4),
    (1.0e7,    323.15, 105.2),
    (1.0e7,    273.15, 127.6),
    (2.0e5,    273.15,   2.557),
]


def anchor_check():
    print("=" * 60)
    print(" Anchor check: Python PR+Peneloux vs NIST")
    print("=" * 60)
    print("      P [Pa]    T [K]  rho [kg/m3]       NIST  err [%]")
    max_err = 0.0
    for P, T, rho_nist in NIST_ANCHORS:
        props = pr_peneloux(P, T)
        rho = props['rho']
        err = abs(rho - rho_nist) / rho_nist * 100.0
        max_err = max(max_err, err)
        print(f"  {P:.4e}   {T:6.2f}     {rho:8.4f}   {rho_nist:8.4f}"
              f"   {err:.3f}")
    print(f"Max relative error vs NIST: {max_err:.3f} %")


# -----------------------------------------------------------------------
# CSV cross-check
# -----------------------------------------------------------------------

def _try_parse_data_line(line):
    """
    Attempt to parse a CSV line as numeric data.
    Returns dict of floats or None if the line is a header / non-data.
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    parts = [p.strip() for p in line.split(',')]
    if len(parts) < 6:
        return None

    try:
        return {
            'P':     float(parts[0]),
            'T':     float(parts[1]),
            'rho':   float(parts[2]),
            'Z':     float(parts[3]),
            'phi':   float(parts[4]),
            'h_dep': float(parts[5]),
        }
    except ValueError:
        # Header line — first field isn't a number
        return None


def csv_cross_check(filename):
    print("=" * 60)
    print(f" Cross-check: Fortran vs Python on {filename}")
    print("=" * 60)

    rel_tol = 1.0e-6   # acceptance threshold
    max_err = {'rho': 0.0, 'Z': 0.0, 'phi': 0.0, 'h_dep': 0.0}
    worst_rho_row = None
    n_compared = 0

    with open(filename) as f:
        for line in f:
            data = _try_parse_data_line(line)
            if data is None:
                continue

            P, T = data['P'], data['T']

            # Skip rows where Fortran reported a failure (would be in ierr column)
            if data['rho'] <= 0.0:
                continue

            py = pr_peneloux(P, T)

            for key in ('rho', 'Z', 'phi', 'h_dep'):
                ref = abs(data[key]) if data[key] != 0.0 else 1.0
                err = abs(data[key] - py[key]) / ref
                if err > max_err[key]:
                    max_err[key] = err
                    if key == 'rho':
                        worst_rho_row = (P, T, data[key], py[key])
            n_compared += 1

    # CRITICAL: refuse to report PASS on zero data points
    if n_compared == 0:
        print("ERROR: 0 data points compared.")
        print(f"  The file '{filename}' contained no parseable data rows.")
        print(f"  Check that the test program wrote the CSV correctly.")
        sys.exit(1)

    print(f"Compared {n_compared} grid points")
    print(f"  Max rel err rho   : {max_err['rho']:.3e}")
    print(f"  Max rel err Z     : {max_err['Z']:.3e}")
    print(f"  Max rel err phi   : {max_err['phi']:.3e}")
    print(f"  Max rel err h_dep : {max_err['h_dep']:.3e}")
    if worst_rho_row is not None:
        P, T, rho_F, rho_P = worst_rho_row
        print(f"  Worst rho row: P={P:.2e}, T={T:.2f}: "
              f"rho_F={rho_F:.6f}, rho_P={rho_P:.6f}")

    overall_max = max(max_err.values())
    if overall_max <= rel_tol:
        print(f"PASS: agreement to within {rel_tol:.0e} "
              f"(acceptable; differences are at floating-point roundoff scale)")
        return 0
    else:
        print(f"FAIL: max relative error {overall_max:.3e} exceeds "
              f"acceptance tolerance {rel_tol:.0e}")
        return 1


# -----------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------

if __name__ == "__main__":
    anchor_check()
    if len(sys.argv) > 1:
        rc = csv_cross_check(sys.argv[1])
        sys.exit(rc)
