! =====================================================================
! Air_EOS_PR_module
!
! Peng-Robinson equation of state for air as a pseudo-pure component,
! with Peneloux volume translation calibrated against NIST air data
! for the supercritical regime relevant to compressed air energy
! storage in aquifers (50-150 bar, 0-80 C).
!
! Intended as a drop-in EOS module for PFLOTRAN's GENERAL mode in
! place of the default Span-Wagner CO2 EOS. Returns gas-phase density,
! compressibility, fugacity coefficient, and enthalpy departure plus
! analytic derivatives with respect to pressure and temperature
! suitable for assembly into the Newton solver Jacobian.
!
! Air pseudo-component critical properties from:
!   Lemmon, E.W., Jacobsen, R.T., Penoncello, S.G., Friend, D.G., 2000.
!   Thermodynamic properties of air and mixtures of nitrogen, argon,
!   and oxygen from 60 to 2000 K at pressures to 2000 MPa.
!   J. Phys. Chem. Ref. Data 29, 331-385.
!
! Peng-Robinson EOS formulation from:
!   Peng, D.-Y., Robinson, D.B., 1976. A new two-constant equation of
!   state. Ind. Eng. Chem. Fundam. 15, 59-64.
!
! Peneloux volume translation from:
!   Peneloux, A., Rauzy, E., Freze, R., 1982. A consistent correction
!   for Redlich-Kwong-Soave volumes. Fluid Phase Equilib. 8, 7-23.
!   (Constant c calibrated against NIST air data, not from the
!   original Rackett-based estimator, which is poorly suited to
!   supercritical light gases.)
!
! PFLOTRAN integration:
!   Compile with -DPFLOTRAN_INTEGRATION to enable PFLOTRAN/PETSc types
!   and pull IDEAL_GAS_CONSTANT from PFLOTRAN_Constants_module instead
!   of the local fallback. The standalone build (no define) compiles
!   and runs against the test driver without any PFLOTRAN dependency.
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

module Air_EOS_PR_module

#ifdef PFLOTRAN_INTEGRATION
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module, only : IDEAL_GAS_CONSTANT
#endif

  implicit none

  private

  ! ------------------------------------------------------------------
  ! Public interface
  ! ------------------------------------------------------------------
  public :: AirEOSPRProperties
  public :: AirEOSPRDensity
  public :: AirEOSPRFugacityCoeff
  public :: AirEOSPREnthalpyDeparture
  public :: AirEOSPRVerifyDerivatives

  ! ------------------------------------------------------------------
  ! Air pseudo-component parameters (Lemmon et al. 2000)
  ! ------------------------------------------------------------------
  PetscReal, parameter :: AIR_TC    = 132.5306d0    ! Critical T [K]
  PetscReal, parameter :: AIR_PC    = 3.7850d6      ! Critical P [Pa]
  PetscReal, parameter :: AIR_OMEGA = 0.0335d0      ! Acentric factor [-]
  PetscReal, parameter :: AIR_MW    = 28.9586d-3    ! Molar mass [kg/mol]

  ! ------------------------------------------------------------------
  ! Peneloux volume translation parameter, calibrated against NIST
  ! air reference data at 50-100 bar, 0-50 C (median of 4 anchor
  ! points; max residual across calibration set < 1.5 %).
  !
  ! V_physical = V_PR - AIR_PENELOUX_C
  !
  ! Negative because standard PR with single-pseudo-component air
  ! parameters UNDERPREDICTS molar volume in the supercritical regime
  ! above air's critical temperature (Tc = 132.5 K). The Peneloux
  ! shift restores agreement to better than 1.5 % below 100 bar.
  !
  ! Set AIR_PENELOUX_C = 0.0d0 to disable translation and recover
  ! pure Peng-Robinson behavior (for code-comparison studies).
  ! ------------------------------------------------------------------
  PetscReal, parameter :: AIR_PENELOUX_C = -9.0d-6  ! [m^3/mol]

  ! ------------------------------------------------------------------
  ! Peng-Robinson universal constants
  ! ------------------------------------------------------------------
  PetscReal, parameter :: PR_OMEGA_A = 0.45723553d0
  PetscReal, parameter :: PR_OMEGA_B = 0.07779607d0

  ! Air-specific kappa parameter, precomputed from acentric factor
  PetscReal, parameter :: PR_KAPPA = 0.37464d0 &
                                   + 1.54226d0 * AIR_OMEGA &
                                   - 0.26992d0 * AIR_OMEGA * AIR_OMEGA

  ! Numerical constants
  PetscReal, parameter :: SQRT_2     = 1.41421356237309504880d0
  PetscReal, parameter :: TWO_SQRT_2 = 2.82842712474619009760d0
  PetscReal, parameter :: PR_LO      = -0.41421356237309504880d0  ! 1 - sqrt(2)
  PetscReal, parameter :: PR_HI      =  2.41421356237309504880d0  ! 1 + sqrt(2)
  PetscReal, parameter :: PI_CONST   = 3.14159265358979323846d0
  PetscReal, parameter :: DISC_TOL   = 1.0d-12  ! cubic discriminant tolerance

  ! Gas constant fallback for standalone build
#ifndef PFLOTRAN_INTEGRATION
  PetscReal, parameter :: IDEAL_GAS_CONSTANT = 8.31446261815324d0  ! J/(mol K)
#endif

contains

! =====================================================================
! Subroutine: AirEOSPRProperties
!
! Main interface returning the full PR + Peneloux air gas-phase
! property set with analytic (P, T) derivatives.
!
! Input:
!   P      [Pa]   Pressure (must be > 0)
!   T      [K]    Temperature (must be > 0)
!
! Output:
!   rho      [kg/m^3]      Physical mass density (with Peneloux)
!   Z        [-]           Physical compressibility (with Peneloux)
!   phi      [-]           Fugacity coefficient (with Peneloux correction)
!   h_dep    [J/mol]       Enthalpy departure h - h^ig (PR, unaffected
!                          by linear volume translation)
!   drho_dP  [kg/m^3/Pa]   dRho/dP at constant T
!   drho_dT  [kg/m^3/K]    dRho/dT at constant P
!   dphi_dP  [1/Pa]        dPhi/dP at constant T
!   dphi_dT  [1/K]         dPhi/dT at constant P
!   ierr     [-]           Error code (0 = success)
! =====================================================================
  subroutine AirEOSPRProperties(P, T, rho, Z, phi, h_dep, &
                                drho_dP, drho_dT, dphi_dP, dphi_dT, &
                                ierr)
    PetscReal, intent(in)       :: P, T
    PetscReal, intent(out)      :: rho, Z, phi, h_dep
    PetscReal, intent(out)      :: drho_dP, drho_dT, dphi_dP, dphi_dT
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: a_T, dadT, b
    PetscReal :: A_PR, B_PR
    PetscReal :: dA_dP, dA_dT_d, dB_dP, dB_dT
    PetscReal :: Z_PR, V_PR, V_phys
    PetscReal :: dV_dP, dV_dT
    PetscReal :: ln_phi_PR, ln_phi
    PetscReal :: dlnphi_dZ, dlnphi_dA, dlnphi_dB
    PetscReal :: dZ_dP, dZ_dT
    PetscReal :: dlnphi_PR_dP, dlnphi_PR_dT
    PetscReal :: dlnphi_dP, dlnphi_dT
    PetscReal :: RT

    ierr = 0

    if (P <= 0.0d0 .or. T <= 0.0d0) then
      ierr = 1
      return
    end if

    RT = IDEAL_GAS_CONSTANT * T

    ! PR attractive (a) and repulsive (b) parameters and da/dT
    call PRComputeAB(T, a_T, dadT, b)

    ! Dimensionless A, B and their (P, T) derivatives
    !   A = a P / (RT)^2,  B = b P / (RT)
    !   dA/dT = P (T da/dT - 2 a) / (R^2 T^3)
    !   dB/dT = -b P / (R T^2)
    A_PR = a_T * P / (RT * RT)
    B_PR = b * P / RT

    dA_dP   = a_T / (RT * RT)
    dA_dT_d = P * (T * dadT - 2.0d0 * a_T) &
              / (IDEAL_GAS_CONSTANT * IDEAL_GAS_CONSTANT * T * T * T)
    dB_dP   = b / RT
    dB_dT   = -b * P / (RT * T)

    ! Solve PR cubic for gas-phase Z (internal, pre-Peneloux)
    call PRSolveCubicGas(A_PR, B_PR, Z_PR, ierr)
    if (ierr /= 0) return

    ! Peneloux volume translation:
    !   V_phys = V_PR - c
    !   Z_phys = Z_PR - c P / (RT)
    !   rho    = MW / V_phys
    V_PR   = Z_PR * RT / P
    V_phys = V_PR - AIR_PENELOUX_C
    Z      = P * V_phys / RT
    rho    = AIR_MW / V_phys

    ! Density derivatives: dV_phys/d(P,T) = dV_PR/d(P,T) since c is
    ! constant; use V_phys in the prefactor for rho.
    call PRMolarVolumeDerivatives(P, T, V_PR, a_T, dadT, b, dV_dP, dV_dT)
    drho_dP = -AIR_MW / (V_phys * V_phys) * dV_dP
    drho_dT = -AIR_MW / (V_phys * V_phys) * dV_dT

    ! Fugacity coefficient with Peneloux correction:
    !   ln(phi_phys) = ln(phi_PR) - c P / (RT)
    ! Standard result for constant volume translation (Peneloux et al.
    ! 1982 Eq. 6; Privat et al. 2012 review).
    ln_phi_PR = PRLnFugacityCoeff(Z_PR, A_PR, B_PR)
    ln_phi    = ln_phi_PR - AIR_PENELOUX_C * P / RT
    phi       = exp(ln_phi)

    ! Fugacity-coefficient derivatives. Compute d(ln phi_PR)/d(P, T)
    ! via the chain rule, then add the Peneloux correction's own
    ! P- and T-derivatives.
    call PRFugacityCoeffPartials(Z_PR, A_PR, B_PR, &
                                 dlnphi_dZ, dlnphi_dA, dlnphi_dB)
    call PRCompressibilityPartialsFromAB(Z_PR, A_PR, B_PR, dZ_dP, dZ_dT, &
                                         dA_dP, dA_dT_d, dB_dP, dB_dT)

    dlnphi_PR_dP = dlnphi_dZ * dZ_dP + dlnphi_dA * dA_dP + dlnphi_dB * dB_dP
    dlnphi_PR_dT = dlnphi_dZ * dZ_dT + dlnphi_dA * dA_dT_d + dlnphi_dB * dB_dT

    ! d/dP [-c P / (RT)] = -c / (RT)
    ! d/dT [-c P / (RT)] =  c P / (R T^2)
    dlnphi_dP = dlnphi_PR_dP - AIR_PENELOUX_C / RT
    dlnphi_dT = dlnphi_PR_dT &
              + AIR_PENELOUX_C * P / (IDEAL_GAS_CONSTANT * T * T)

    dphi_dP = phi * dlnphi_dP
    dphi_dT = phi * dlnphi_dT

    ! Enthalpy departure (PR formulation, uses internal Z_PR).
    ! Linear volume translation does not change h - h^ig because
    ! (dP/dT)_V is invariant under the shift.
    h_dep = PREnthalpyDepartureCore(T, Z_PR, A_PR, B_PR, a_T, dadT, b)

  end subroutine AirEOSPRProperties

! =====================================================================
! Lightweight density-only interface (with Peneloux)
! =====================================================================
  subroutine AirEOSPRDensity(P, T, rho, drho_dP, drho_dT, ierr)
    PetscReal, intent(in)       :: P, T
    PetscReal, intent(out)      :: rho, drho_dP, drho_dT
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: a_T, dadT, b
    PetscReal :: A_PR, B_PR, Z_PR, V_PR, V_phys, RT
    PetscReal :: dV_dP, dV_dT

    ierr = 0
    if (P <= 0.0d0 .or. T <= 0.0d0) then
      ierr = 1
      return
    end if

    RT = IDEAL_GAS_CONSTANT * T
    call PRComputeAB(T, a_T, dadT, b)
    A_PR = a_T * P / (RT * RT)
    B_PR = b * P / RT

    call PRSolveCubicGas(A_PR, B_PR, Z_PR, ierr)
    if (ierr /= 0) return

    V_PR   = Z_PR * RT / P
    V_phys = V_PR - AIR_PENELOUX_C
    rho    = AIR_MW / V_phys

    call PRMolarVolumeDerivatives(P, T, V_PR, a_T, dadT, b, dV_dP, dV_dT)
    drho_dP = -AIR_MW / (V_phys * V_phys) * dV_dP
    drho_dT = -AIR_MW / (V_phys * V_phys) * dV_dT
  end subroutine AirEOSPRDensity

! =====================================================================
! Fugacity coefficient (with Peneloux); returns -1 on failure
! =====================================================================
  function AirEOSPRFugacityCoeff(P, T, ierr) result(phi)
    PetscReal, intent(in)       :: P, T
    PetscErrorCode, intent(out) :: ierr
    PetscReal                   :: phi
    PetscReal :: a_T, dadT, b, A_PR, B_PR, Z_PR, RT, ln_phi

    ierr = 0
    phi  = -1.0d0
    if (P <= 0.0d0 .or. T <= 0.0d0) then
      ierr = 1
      return
    end if

    RT = IDEAL_GAS_CONSTANT * T
    call PRComputeAB(T, a_T, dadT, b)
    A_PR = a_T * P / (RT * RT)
    B_PR = b * P / RT

    call PRSolveCubicGas(A_PR, B_PR, Z_PR, ierr)
    if (ierr /= 0) return

    ln_phi = PRLnFugacityCoeff(Z_PR, A_PR, B_PR) &
           - AIR_PENELOUX_C * P / RT
    phi    = exp(ln_phi)
  end function AirEOSPRFugacityCoeff

! =====================================================================
! Enthalpy departure h - h^ig in J/mol (PR; Peneloux-invariant)
! =====================================================================
  function AirEOSPREnthalpyDeparture(P, T, ierr) result(h_dep)
    PetscReal, intent(in)       :: P, T
    PetscErrorCode, intent(out) :: ierr
    PetscReal                   :: h_dep
    PetscReal :: a_T, dadT, b, A_PR, B_PR, Z_PR, RT

    ierr  = 0
    h_dep = 0.0d0
    if (P <= 0.0d0 .or. T <= 0.0d0) then
      ierr = 1
      return
    end if

    RT = IDEAL_GAS_CONSTANT * T
    call PRComputeAB(T, a_T, dadT, b)
    A_PR = a_T * P / (RT * RT)
    B_PR = b * P / RT
    call PRSolveCubicGas(A_PR, B_PR, Z_PR, ierr)
    if (ierr /= 0) return
    h_dep = PREnthalpyDepartureCore(T, Z_PR, A_PR, B_PR, a_T, dadT, b)
  end function AirEOSPREnthalpyDeparture

! =====================================================================
! Diagnostic: analytic vs finite-difference derivative check
! =====================================================================
  subroutine AirEOSPRVerifyDerivatives(P_ref, T_ref, max_rel_err, ierr)
    PetscReal, intent(in)       :: P_ref, T_ref
    PetscReal, intent(out)      :: max_rel_err
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: rho0, Z0, phi0, h0
    PetscReal :: drho_dP_a, drho_dT_a, dphi_dP_a, dphi_dT_a
    PetscReal :: rho_p, rho_m, phi_p, phi_m
    PetscReal :: Z_dum, h_dum, d1, d2, d3, d4
    PetscReal :: drho_dP_fd, drho_dT_fd, dphi_dP_fd, dphi_dT_fd
    PetscReal :: dP, dT, err
    PetscErrorCode :: ie

    ierr = 0
    max_rel_err = 0.0d0

    call AirEOSPRProperties(P_ref, T_ref, rho0, Z0, phi0, h0, &
                            drho_dP_a, drho_dT_a, dphi_dP_a, dphi_dT_a, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    dP = max(1.0d0, 1.0d-6 * P_ref)
    dT = max(1.0d-3, 1.0d-6 * T_ref)

    call AirEOSPRProperties(P_ref + dP, T_ref, rho_p, Z_dum, phi_p, h_dum, &
                            d1, d2, d3, d4, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    call AirEOSPRProperties(P_ref - dP, T_ref, rho_m, Z_dum, phi_m, h_dum, &
                            d1, d2, d3, d4, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    drho_dP_fd = (rho_p - rho_m) / (2.0d0 * dP)
    dphi_dP_fd = (phi_p - phi_m) / (2.0d0 * dP)

    call AirEOSPRProperties(P_ref, T_ref + dT, rho_p, Z_dum, phi_p, h_dum, &
                            d1, d2, d3, d4, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    call AirEOSPRProperties(P_ref, T_ref - dT, rho_m, Z_dum, phi_m, h_dum, &
                            d1, d2, d3, d4, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    drho_dT_fd = (rho_p - rho_m) / (2.0d0 * dT)
    dphi_dT_fd = (phi_p - phi_m) / (2.0d0 * dT)

    err = abs(drho_dP_a - drho_dP_fd) / max(abs(drho_dP_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(drho_dT_a - drho_dT_fd) / max(abs(drho_dT_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(dphi_dP_a - dphi_dP_fd) / max(abs(dphi_dP_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(dphi_dT_a - dphi_dT_fd) / max(abs(dphi_dT_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
  end subroutine AirEOSPRVerifyDerivatives

! =====================================================================
! ====               PRIVATE HELPER ROUTINES                        ===
! =====================================================================

  subroutine PRComputeAB(T, a_T, dadT, b)
    PetscReal, intent(in)  :: T
    PetscReal, intent(out) :: a_T, dadT, b
    PetscReal :: a_c, sqrt_Tr, u, alpha_T, dalpha_dT

    a_c = PR_OMEGA_A * IDEAL_GAS_CONSTANT**2 * AIR_TC**2 / AIR_PC
    b   = PR_OMEGA_B * IDEAL_GAS_CONSTANT * AIR_TC / AIR_PC

    sqrt_Tr   = sqrt(T / AIR_TC)
    u         = 1.0d0 + PR_KAPPA * (1.0d0 - sqrt_Tr)
    alpha_T   = u * u
    dalpha_dT = -u * PR_KAPPA / sqrt(T * AIR_TC)

    a_T  = a_c * alpha_T
    dadT = a_c * dalpha_dT
  end subroutine PRComputeAB

! ---------------------------------------------------------------------
! PRSolveCubicGas: solve the dimensionless PR cubic for the gas root.
! Hardened against the disc ~ 0 boundary via the trig branch and
! clamped acos argument.
! ---------------------------------------------------------------------
  subroutine PRSolveCubicGas(A_PR, B_PR, Z_gas, ierr)
    PetscReal, intent(in)       :: A_PR, B_PR
    PetscReal, intent(out)      :: Z_gas
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: c2, c1, c0
    PetscReal :: p, q, disc, shift
    PetscReal :: w1, w2, w3, Z_candidates(3)
    PetscReal :: r_arg, theta, sqrt_mp3
    PetscReal :: sign_arg, u_re, v_re
    PetscInt  :: i
    logical   :: found

    ierr = 0

    c2 = -(1.0d0 - B_PR)
    c1 =  A_PR - 3.0d0 * B_PR * B_PR - 2.0d0 * B_PR
    c0 = -(A_PR * B_PR - B_PR * B_PR - B_PR * B_PR * B_PR)

    shift = c2 / 3.0d0
    p     = c1 - c2 * c2 / 3.0d0
    q     = 2.0d0 * c2**3 / 27.0d0 - c2 * c1 / 3.0d0 + c0

    disc = (q * q) / 4.0d0 + (p * p * p) / 27.0d0

    if (disc > DISC_TOL) then
      sign_arg = -q / 2.0d0 + sqrt(disc)
      u_re     = sign(abs(sign_arg)**(1.0d0/3.0d0), sign_arg)
      sign_arg = -q / 2.0d0 - sqrt(disc)
      v_re     = sign(abs(sign_arg)**(1.0d0/3.0d0), sign_arg)
      Z_gas    = u_re + v_re - shift
    else
      sqrt_mp3 = sqrt(max(-p / 3.0d0, 0.0d0))
      if (sqrt_mp3 < tiny(1.0d0)) then
        Z_gas = -shift
      else
        r_arg = (-q / 2.0d0) / (sqrt_mp3**3)
        if (r_arg >  1.0d0) r_arg =  1.0d0
        if (r_arg < -1.0d0) r_arg = -1.0d0
        theta = acos(r_arg)
        w1 = 2.0d0 * sqrt_mp3 * cos( theta / 3.0d0 )
        w2 = 2.0d0 * sqrt_mp3 * cos((theta + 2.0d0 * PI_CONST) / 3.0d0 )
        w3 = 2.0d0 * sqrt_mp3 * cos((theta + 4.0d0 * PI_CONST) / 3.0d0 )
        Z_candidates(1) = w1 - shift
        Z_candidates(2) = w2 - shift
        Z_candidates(3) = w3 - shift
        found = .false.
        Z_gas = 0.0d0
        do i = 1, 3
          if (Z_candidates(i) > B_PR) then
            if (.not. found .or. Z_candidates(i) > Z_gas) then
              Z_gas = Z_candidates(i)
              found = .true.
            end if
          end if
        end do
        if (.not. found) then
          ierr = 2
          return
        end if
      end if
    end if

    if (Z_gas <= B_PR) then
      ierr = 3
      return
    end if
  end subroutine PRSolveCubicGas

  subroutine PRMolarVolumeDerivatives(P, T, V_mol, a_T, dadT, b, &
                                       dV_dP, dV_dT)
    PetscReal, intent(in)  :: P, T, V_mol, a_T, dadT, b
    PetscReal, intent(out) :: dV_dP, dV_dT
    PetscReal :: V_minus_b, denom_attr, df_dV, df_dT

    V_minus_b  = V_mol - b
    denom_attr = V_mol * V_mol + 2.0d0 * b * V_mol - b * b

    df_dV = -IDEAL_GAS_CONSTANT * T / (V_minus_b * V_minus_b) &
          + a_T * (2.0d0 * V_mol + 2.0d0 * b) / (denom_attr * denom_attr)
    df_dT =  IDEAL_GAS_CONSTANT / V_minus_b - dadT / denom_attr

    dV_dP =  1.0d0 / df_dV
    dV_dT = -df_dT / df_dV
  end subroutine PRMolarVolumeDerivatives

  function PRLnFugacityCoeff(Z, A_PR, B_PR) result(ln_phi)
    PetscReal, intent(in) :: Z, A_PR, B_PR
    PetscReal             :: ln_phi
    PetscReal             :: num, den

    num    = Z + PR_HI * B_PR
    den    = Z + PR_LO * B_PR
    ln_phi = (Z - 1.0d0) - log(Z - B_PR) &
           - A_PR / (TWO_SQRT_2 * B_PR) * log(num / den)
  end function PRLnFugacityCoeff

  subroutine PRFugacityCoeffPartials(Z, A_PR, B_PR, &
                                     dlnphi_dZ, dlnphi_dA, dlnphi_dB)
    PetscReal, intent(in)  :: Z, A_PR, B_PR
    PetscReal, intent(out) :: dlnphi_dZ, dlnphi_dA, dlnphi_dB
    PetscReal :: num, den, L, dL_dZ, dL_dB, inv_2s2B

    num   = Z + PR_HI * B_PR
    den   = Z + PR_LO * B_PR
    L     = log(num / den)
    dL_dZ = 1.0d0 / num - 1.0d0 / den
    dL_dB = PR_HI / num - PR_LO / den

    inv_2s2B = 1.0d0 / (TWO_SQRT_2 * B_PR)

    dlnphi_dZ = 1.0d0 - 1.0d0 / (Z - B_PR) - A_PR * inv_2s2B * dL_dZ
    dlnphi_dA = -inv_2s2B * L
    dlnphi_dB =  1.0d0 / (Z - B_PR) &
              +  A_PR / (TWO_SQRT_2 * B_PR * B_PR) * L &
              -  A_PR * inv_2s2B * dL_dB
  end subroutine PRFugacityCoeffPartials

  subroutine PRCompressibilityPartialsFromAB(Z, A_PR, B_PR, dZ_dP, dZ_dT, &
                                             dA_dP, dA_dT, dB_dP, dB_dT)
    PetscReal, intent(in)  :: Z, A_PR, B_PR
    PetscReal, intent(in)  :: dA_dP, dA_dT, dB_dP, dB_dT
    PetscReal, intent(out) :: dZ_dP, dZ_dT
    PetscReal :: c2_loc, dF_dZ, dF_dA, dF_dB, dZ_dA, dZ_dB

    c2_loc = -(1.0d0 - B_PR)
    dF_dZ  = 3.0d0 * Z * Z + 2.0d0 * c2_loc * Z &
           + (A_PR - 3.0d0 * B_PR * B_PR - 2.0d0 * B_PR)
    dF_dA  = Z - B_PR
    dF_dB  = Z * Z + ( -6.0d0 * B_PR - 2.0d0 ) * Z &
                   + ( -A_PR + 2.0d0 * B_PR + 3.0d0 * B_PR * B_PR )

    dZ_dA = -dF_dA / dF_dZ
    dZ_dB = -dF_dB / dF_dZ

    dZ_dP = dZ_dA * dA_dP + dZ_dB * dB_dP
    dZ_dT = dZ_dA * dA_dT + dZ_dB * dB_dT
  end subroutine PRCompressibilityPartialsFromAB

  function PREnthalpyDepartureCore(T, Z, A_PR, B_PR, a_T, dadT, b) result(h_R)
    PetscReal, intent(in) :: T, Z, A_PR, B_PR, a_T, dadT, b
    PetscReal             :: h_R
    PetscReal             :: num, den

    num = Z + PR_HI * B_PR
    den = Z + PR_LO * B_PR
    h_R = IDEAL_GAS_CONSTANT * T * (Z - 1.0d0) &
        + (T * dadT - a_T) / (TWO_SQRT_2 * b) * log(num / den)
  end function PREnthalpyDepartureCore

end module Air_EOS_PR_module
