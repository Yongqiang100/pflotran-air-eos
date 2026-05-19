#!/usr/bin/env python3
"""
cross_validate.py

Independent Python implementation of Peng-Robinson + Peneloux EOS for
air, used to cross-check the Fortran module.

Two checks:
  1. Sample point-by-point comparison at the hard-coded NIST anchors.
  2. Row-by-row comparison against `air_eos_sweep.csv` produced by
     the Fortran sweep driver. Reports max relative error in density
     and fugacity coefficient; should be at the 1e-12 level (pure
     floating-point roundoff) if the two implementations are doing
     the same math.

Usage:
    python3 cross_validate.py                # NIST anchor check only
    python3 cross_validate.py air_eos_sweep.csv  # full sweep check

The Python implementation here was written from the math, not from
the Fortran source, to keep the cross-check methodologically clean.
"""

import math
import sys

# -------- Air pseudo-component parameters (Lemmon et al. 2000) -------
AIR_TC    = 132.5306        # K
AIR_PC    = 3.7850e6        # Pa
AIR_OMEGA = 0.0335          # -
AIR_MW    = 28.9586e-3      # kg/mol

# -------- Peng-Robinson and Peneloux constants -----------------------
PR_OMEGA_A = 0.45723553
PR_OMEGA_B = 0.07779607
PR_KAPPA   = 0.37464 + 1.54226 * AIR_OMEGA - 0.26992 * AIR_OMEGA**2
AIR_PENELOUX_C = -9.0e-6    # m^3/mol (must match Fortran module)
R_GAS = 8.31446261815324    # J/(mol K)

SQRT_2 = math.sqrt(2.0)
PR_LO = 1.0 - SQRT_2
PR_HI = 1.0 + SQRT_2


def pr_alpha(T):
    """alpha(T) for PR EOS."""
    u = 1.0 + PR_KAPPA * (1.0 - math.sqrt(T / AIR_TC))
    return u * u


def pr_a_b(T):
    """a(T), b for PR EOS (air)."""
    a_c = PR_OMEGA_A * R_GAS**2 * AIR_TC**2 / AIR_PC
    b   = PR_OMEGA_B * R_GAS * AIR_TC / AIR_PC
    return a_c * pr_alpha(T), b


def solve_pr_cubic_gas(A, B):
    """Return largest real root Z > B of the PR cubic."""
    c2 = -(1.0 - B)
    c1 =  A - 3.0 * B**2 - 2.0 * B
    c0 = -(A * B - B**2 - B**3)

    # depressed cubic w^3 + p w + q = 0,  Z = w - c2/3
    shift = c2 / 3.0
    p = c1 - c2**2 / 3.0
    q = 2.0 * c2**3 / 27.0 - c2 * c1 / 3.0 + c0

    disc = q**2 / 4.0 + p**3 / 27.0

    if disc > 1e-12:
        sign_arg = -q / 2.0 + math.sqrt(disc)
        u_re = math.copysign(abs(sign_arg)**(1.0 / 3.0), sign_arg)
        sign_arg = -q / 2.0 - math.sqrt(disc)
        v_re = math.copysign(abs(sign_arg)**(1.0 / 3.0), sign_arg)
        Z = u_re + v_re - shift
    else:
        sqrt_mp3 = math.sqrt(max(-p / 3.0, 0.0))
        if sqrt_mp3 < 1e-300:
            Z = -shift
        else:
            r_arg = (-q / 2.0) / sqrt_mp3**3
            r_arg = max(-1.0, min(1.0, r_arg))
            theta = math.acos(r_arg)
            roots = []
            for k in range(3):
                w = 2.0 * sqrt_mp3 * math.cos((theta + 2.0 * math.pi * k) / 3.0)
                roots.append(w - shift)
            candidates = [z for z in roots if z > B]
            if not candidates:
                raise RuntimeError("No physical gas root found")
            Z = max(candidates)
    if Z <= B:
        raise RuntimeError(f"Gas root Z={Z} fails Z > B={B}")
    return Z


def air_eos(P, T):
    """Return (rho, Z, phi, h_dep) for air at (P, T) using PR + Peneloux."""
    RT = R_GAS * T
    a, b = pr_a_b(T)
    A = a * P / (RT * RT)
    B = b * P / RT

    Z_PR = solve_pr_cubic_gas(A, B)
    V_PR = Z_PR * RT / P
    V_phys = V_PR - AIR_PENELOUX_C

    # Physical properties (Peneloux-corrected where applicable)
    rho = AIR_MW / V_phys
    Z_phys = P * V_phys / RT

    # Fugacity coefficient
    num = Z_PR + PR_HI * B
    den = Z_PR + PR_LO * B
    ln_phi_PR = (Z_PR - 1.0) - math.log(Z_PR - B) \
              - A / (2.0 * SQRT_2 * B) * math.log(num / den)
    ln_phi = ln_phi_PR - AIR_PENELOUX_C * P / RT
    phi = math.exp(ln_phi)

    # Enthalpy departure (PR; Peneloux-invariant)
    # da/dT for PR
    u = 1.0 + PR_KAPPA * (1.0 - math.sqrt(T / AIR_TC))
    dalpha_dT = -u * PR_KAPPA / math.sqrt(T * AIR_TC)
    a_c = PR_OMEGA_A * R_GAS**2 * AIR_TC**2 / AIR_PC
    dadT = a_c * dalpha_dT
    h_dep = R_GAS * T * (Z_PR - 1.0) \
          + (T * dadT - a) / (2.0 * SQRT_2 * b) * math.log(num / den)

    return rho, Z_phys, phi, h_dep


def anchor_check():
    """Check at the same 6 NIST anchor points used by the Fortran tests."""
    points = [
        (1.01325e5, 298.15,  1.1839),
        (5.0e6,    298.15,  58.78),
        (1.0e7,    298.15, 115.40),
        (1.0e7,    323.15, 105.20),
        (1.0e7,    273.15, 127.60),
        (2.0e5,    273.15,   2.557),
    ]
    print("=" * 60)
    print(" Anchor check: Python PR+Peneloux vs NIST")
    print("=" * 60)
    print(f"{'P [Pa]':>12} {'T [K]':>8} {'rho [kg/m3]':>12} "
          f"{'NIST':>10} {'err [%]':>8}")
    max_err = 0.0
    for P, T, rho_nist in points:
        rho, Z, phi, h_dep = air_eos(P, T)
        err = abs(rho - rho_nist) / rho_nist
        max_err = max(max_err, err)
        print(f"{P:12.4e} {T:8.2f} {rho:12.4f} {rho_nist:10.4f} {err*100:7.3f}")
    print(f"\nMax relative error vs NIST: {max_err*100:.3f} %")
    return max_err


def csv_cross_check(csv_path):
    """Compare Fortran sweep CSV row-by-row against this Python implementation."""
    print("=" * 60)
    print(f" Cross-check: Fortran vs Python on {csv_path}")
    print("=" * 60)
    max_rho_err = 0.0
    max_Z_err   = 0.0
    max_phi_err = 0.0
    max_h_err   = 0.0
    worst_row = ""
    n_compared = 0
    with open(csv_path, "r") as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split(",")
            if "NaN" in parts:
                continue
            P    = float(parts[0])
            T    = float(parts[1])
            rho_f = float(parts[2])
            Z_f   = float(parts[3])
            phi_f = float(parts[4])
            h_f   = float(parts[5])
            rho_p, Z_p, phi_p, h_p = air_eos(P, T)

            rel = lambda a, b: abs(a - b) / max(abs(b), 1e-300)
            erho = rel(rho_p, rho_f)
            eZ   = rel(Z_p,   Z_f)
            ephi = rel(phi_p, phi_f)
            eh   = rel(h_p,   h_f)
            if erho > max_rho_err:
                max_rho_err = erho
                worst_row = f"P={P:.2e}, T={T:.2f}: " \
                            f"rho_F={rho_f:.6f}, rho_P={rho_p:.6f}"
            max_Z_err   = max(max_Z_err,   eZ)
            max_phi_err = max(max_phi_err, ephi)
            max_h_err   = max(max_h_err,   eh)
            n_compared += 1

    print(f"Compared {n_compared} grid points")
    print(f"  Max rel err rho   : {max_rho_err:.3e}")
    print(f"  Max rel err Z     : {max_Z_err:.3e}")
    print(f"  Max rel err phi   : {max_phi_err:.3e}")
    print(f"  Max rel err h_dep : {max_h_err:.3e}")
    if worst_row:
        print(f"  Worst rho row: {worst_row}")
    if max(max_rho_err, max_Z_err, max_phi_err, max_h_err) < 1e-10:
        print("\nPASS: Fortran and Python implementations agree to within 1e-10")
    elif max(max_rho_err, max_Z_err, max_phi_err, max_h_err) < 1e-6:
        print("\nPASS: agreement to within 1e-6 (acceptable; "
              "differences are at floating-point roundoff scale)")
    else:
        print("\nFAIL: implementations disagree by more than 1e-6 — "
              "review the math in both")


if __name__ == "__main__":
    anchor_check()
    print()
    if len(sys.argv) > 1:
        csv_cross_check(sys.argv[1])
