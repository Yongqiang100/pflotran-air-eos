! =====================================================================
! Air_WaterVapour_module
!
! Water vapour partitioning between the gas and aqueous phases for
! compressed-air energy storage in saline aquifers. Computes the
! water mole fraction y_w in the gas phase given total pressure,
! temperature, brine salinity, and the water-vapour fugacity
! coefficient.
!
! Working equation (Modified Raoult's law with Poynting correction):
!   y_w * phi_w * P = a_w * P_sat(T) * Poynting(P, T)
!
! Components:
!   1. P_sat(T) - IAPWS-IF97 Region 4 saturation pressure equation
!      (Wagner et al. 2000). Valid 273.15 K to 647.096 K (critical).
!   2. a_w(m_NaCl) - water activity from ideal mixing on the molar
!      scale: a_w = 55.51 / (55.51 + nu*m), nu = 2 for NaCl. Accurate
!      to ~ 1 % below m_NaCl = 2 mol/kg vs Pitzer; adequate for CAES
!      brines. Can be replaced with Pitzer when extending to higher
!      salinity.
!   3. Poynting(P, T) = exp[V_w (P - P_sat) / (R T)], where V_w is the
!      molar volume of liquid water (~ 18 cm^3/mol).
!
! Why this matters for your science:
!   When dry compressed air is injected into a brine-saturated aquifer,
!   water evaporates into the gas phase, concentrating the residual
!   brine. If the brine reaches halite saturation (~ 6.1 mol/kg at
!   25 C), halite precipitates near the wellbore, altering porosity
!   and permeability. This effect was documented in the PG&E pilot
!   (Medeiros et al. 2018) and is a known concern for CAESA performance.
!   This module provides the y_w that drives that calculation.
!
! References:
!   Wagner, W., Cooper, J.R., Dittmann, A., et al., 2000. The IAPWS
!     Industrial Formulation 1997 for the Thermodynamic Properties of
!     Water and Steam. J. Eng. Gas Turbines Power 122, 150-184.
!   Pitzer, K.S., 1973. Thermodynamics of electrolytes. I. Theoretical
!     basis and general equations. J. Phys. Chem. 77, 268-277.
!   Spycher, N., Pruess, K., 2005. CO2-H2O mixtures in the geological
!     sequestration of CO2. II. Partitioning in chloride brines at
!     12-100 C and up to 600 bar. Geochim. Cosmochim. Acta 69, 3309-3320.
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

module Air_WaterVapour_module

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
  public :: WaterSatPressure
  public :: WaterActivityNaCl
  public :: WaterPoyntingFactor
  public :: WaterVapourMoleFraction
  public :: WaterMolarVolumeLiquid
  public :: WaterVapourVerifyDerivatives

  ! ------------------------------------------------------------------
  ! IAPWS-IF97 Region 4 coefficients (10 numerical constants n1..n10).
  ! From IAPWS-IF97 (2007 release), Equation 30.
  ! Saturation pressure: A*beta^2 + B*beta + C = 0,
  ! beta = (P_sat / P*)^(1/4), where the coefficients depend on
  ! theta = T/T* + n9/(T/T* - n10).
  ! ------------------------------------------------------------------
  PetscReal, parameter :: IAPWS_N1  =  0.11670521452767d4
  PetscReal, parameter :: IAPWS_N2  = -0.72421316703206d6
  PetscReal, parameter :: IAPWS_N3  = -0.17073846940092d2
  PetscReal, parameter :: IAPWS_N4  =  0.12020824702470d5
  PetscReal, parameter :: IAPWS_N5  = -0.32325550322333d7
  PetscReal, parameter :: IAPWS_N6  =  0.14915108613530d2
  PetscReal, parameter :: IAPWS_N7  = -0.48232657361591d4
  PetscReal, parameter :: IAPWS_N8  =  0.40511340542057d6
  PetscReal, parameter :: IAPWS_N9  = -0.23855557567849d0
  PetscReal, parameter :: IAPWS_N10 =  0.65017534844798d3

  PetscReal, parameter :: IAPWS_T_STAR = 1.0d0         ! K
  PetscReal, parameter :: IAPWS_P_STAR = 1.0d6         ! Pa (1 MPa)

  ! Valid range for IAPWS-IF97 Region 4
  PetscReal, parameter :: T_MIN_IAPWS = 273.15d0       ! K
  PetscReal, parameter :: T_MAX_IAPWS = 647.096d0      ! K (critical)

  ! ------------------------------------------------------------------
  ! Water properties (held constant in this module)
  ! ------------------------------------------------------------------
  PetscReal, parameter :: WATER_MW          = 18.01528d-3   ! kg/mol
  PetscReal, parameter :: WATER_MOLALITY_REF = 55.50843d0   ! mol per kg H2O
  PetscReal, parameter :: WATER_V_LIQUID    = 1.807d-5      ! m^3/mol (~ at 25 C, 1 atm)

  ! NaCl dissociation
  PetscReal, parameter :: NU_NACL = 2.0d0   ! Na+ and Cl-

  ! Gas constant fallback
#ifndef PFLOTRAN_INTEGRATION
  PetscReal, parameter :: IDEAL_GAS_CONSTANT = 8.31446261815324d0  ! J/(mol K)
#endif

contains

! =====================================================================
! Function: WaterSatPressure
!
! Returns P_sat(T) [Pa] and its T derivative via IAPWS-IF97 Region 4.
!
! Input:
!   T          [K]    Temperature (must lie in [273.15, 647.096])
!
! Output:
!   P_sat      [Pa]
!   dPsat_dT   [Pa/K]
!   ierr       [-]    0 success; 1 out of range
! =====================================================================
  subroutine WaterSatPressure(T, P_sat, dPsat_dT, ierr)
    PetscReal, intent(in)       :: T
    PetscReal, intent(out)      :: P_sat, dPsat_dT
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: theta, A, B, C, disc, beta
    PetscReal :: dtheta_dT
    PetscReal :: dA_dtheta, dB_dtheta, dC_dtheta
    PetscReal :: dA_dT, dB_dT, dC_dT
    PetscReal :: ddisc_dT, dbeta_dT
    PetscReal :: denom

    ierr     = 0
    P_sat    = 0.0d0
    dPsat_dT = 0.0d0

    if (T < T_MIN_IAPWS .or. T > T_MAX_IAPWS) then
      ierr = 1
      return
    end if

    ! theta = T/T* + n9/(T/T* - n10)
    !       = T + n9/(T - n10) at T_STAR = 1
    theta = T + IAPWS_N9 / (T - IAPWS_N10)

    A = theta * theta + IAPWS_N1 * theta + IAPWS_N2
    B = IAPWS_N3 * theta * theta + IAPWS_N4 * theta + IAPWS_N5
    C = IAPWS_N6 * theta * theta + IAPWS_N7 * theta + IAPWS_N8

    disc = B * B - 4.0d0 * A * C
    if (disc < 0.0d0) then
      ierr = 2
      return
    end if

    ! IAPWS recommended form: beta = 2C / (-B + sqrt(B^2 - 4AC))
    denom = -B + sqrt(disc)
    beta = 2.0d0 * C / denom
    P_sat = IAPWS_P_STAR * beta**4

    ! ----- Derivative w.r.t. T -----
    !
    ! d(theta)/dT = 1 - n9/(T - n10)^2
    dtheta_dT = 1.0d0 - IAPWS_N9 / ((T - IAPWS_N10) * (T - IAPWS_N10))

    dA_dtheta = 2.0d0 * theta + IAPWS_N1
    dB_dtheta = 2.0d0 * IAPWS_N3 * theta + IAPWS_N4
    dC_dtheta = 2.0d0 * IAPWS_N6 * theta + IAPWS_N7

    dA_dT = dA_dtheta * dtheta_dT
    dB_dT = dB_dtheta * dtheta_dT
    dC_dT = dC_dtheta * dtheta_dT

    ddisc_dT = 2.0d0 * B * dB_dT - 4.0d0 * (dA_dT * C + A * dC_dT)

    ! d(beta)/dT via quotient rule applied to beta = 2C / (-B + sqrt(disc))
    ! d(-B + sqrt(disc))/dT = -dB_dT + ddisc_dT / (2 sqrt(disc))
    dbeta_dT = (2.0d0 * dC_dT * denom &
              - 2.0d0 * C * (-dB_dT + ddisc_dT / (2.0d0 * sqrt(disc)))) &
              / (denom * denom)

    ! dP_sat/dT = P* * 4 * beta^3 * d(beta)/dT
    dPsat_dT = IAPWS_P_STAR * 4.0d0 * beta**3 * dbeta_dT

  end subroutine WaterSatPressure

! =====================================================================
! Function: WaterActivityNaCl
!
! Returns water activity a_w in NaCl brine via ideal-mixing on the
! molar scale: a_w = 55.51 / (55.51 + nu*m), nu = 2 for NaCl.
!
! Accuracy: within 1 % of Pitzer model for m < 2 mol/kg, within
! 5 % for m < 5 mol/kg. For higher salinities, replace with a
! proper Pitzer implementation.
!
! Note: This formulation has no explicit T dependence beyond the
! near-temperature-independence of pure-water reference molality.
! Acceptable for CAES T range (0-80 C).
! =====================================================================
  subroutine WaterActivityNaCl(m_NaCl, a_w, da_dm, ierr)
    PetscReal, intent(in)       :: m_NaCl
    PetscReal, intent(out)      :: a_w, da_dm
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: denom

    ierr = 0
    a_w  = 1.0d0
    da_dm = 0.0d0

    if (m_NaCl < 0.0d0) then
      ierr = 1
      return
    end if

    denom = WATER_MOLALITY_REF + NU_NACL * m_NaCl
    a_w   = WATER_MOLALITY_REF / denom
    da_dm = -WATER_MOLALITY_REF * NU_NACL / (denom * denom)

  end subroutine WaterActivityNaCl

! =====================================================================
! Function: WaterPoyntingFactor
!
! Returns the Poynting correction factor for water vapour pressure
! at elevated total pressure:
!   Poynting(P, T) = exp[V_w (P - P_sat) / (R T)]
!
! At 100 bar and 40 C, Poynting ~ 1.07, a 7 % effect — non-trivial
! for CAES storage conditions.
! =====================================================================
  subroutine WaterPoyntingFactor(P, T, P_sat, V_w, &
                                  poynting, dpy_dP, dpy_dT, ierr)
    PetscReal, intent(in)       :: P, T, P_sat, V_w
    PetscReal, intent(out)      :: poynting, dpy_dP, dpy_dT
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: arg, RT

    ierr = 0
    poynting = 1.0d0
    dpy_dP   = 0.0d0
    dpy_dT   = 0.0d0

    if (T <= 0.0d0 .or. P <= 0.0d0) then
      ierr = 1
      return
    end if

    RT = IDEAL_GAS_CONSTANT * T
    arg = V_w * (P - P_sat) / RT
    poynting = exp(arg)

    ! d(arg)/dP = V_w / RT  (P_sat treated as constant here; chain rule
    ! through P_sat handled in WaterVapourMoleFraction)
    dpy_dP = poynting * V_w / RT

    ! d(arg)/dT = -V_w (P - P_sat) / (R T^2)  (again, P_sat held;
    ! chain rule via P_sat applied externally)
    dpy_dT = poynting * (-V_w * (P - P_sat) / (RT * T))

  end subroutine WaterPoyntingFactor

! =====================================================================
! Function: WaterMolarVolumeLiquid
!
! Returns molar volume of liquid water [m^3/mol]. Held constant in
! this module at the 25 C, 1 atm value. Variation across CAES T
! range is ~ 3 %, which propagates to ~ 0.2 % in y_w; negligible
! relative to the IAPWS-IF97 precision.
! =====================================================================
  function WaterMolarVolumeLiquid(T) result(V_w)
    PetscReal, intent(in) :: T
    PetscReal             :: V_w
    PetscReal             :: T_unused
    T_unused = T   ! avoid unused-arg warning; placeholder for future T dependence
    V_w = WATER_V_LIQUID
  end function WaterMolarVolumeLiquid

! =====================================================================
! Subroutine: WaterVapourMoleFraction
!
! Returns the equilibrium water vapour mole fraction y_w in the gas
! phase, given total pressure P, temperature T, brine salinity, and
! the water-vapour fugacity coefficient phi_w.
!
! Working equation:
!   y_w = a_w(m) * P_sat(T) * Poynting(P, T) / (phi_w * P)
!
! Input:
!   P            [Pa]
!   T            [K]
!   m_NaCl       [mol/kg]   Brine NaCl molality
!   phi_w        [-]        Water-vapour fugacity coefficient
!                           (pass 1.0 for ideal gas; in practice use
!                            the bulk-gas fugacity coefficient from
!                            Air_EOS_PR_module as an approximation)
!
! Output:
!   y_w          [-]        Water mole fraction in gas phase
!   dyw_dP       [1/Pa]
!   dyw_dT       [1/K]
!   dyw_dm       [kg/mol]
!   ierr         [-]
! =====================================================================
  subroutine WaterVapourMoleFraction(P, T, m_NaCl, phi_w, &
                                      y_w, dyw_dP, dyw_dT, dyw_dm, ierr)
    PetscReal, intent(in)       :: P, T, m_NaCl, phi_w
    PetscReal, intent(out)      :: y_w, dyw_dP, dyw_dT, dyw_dm
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: P_sat, dPsat_dT
    PetscReal :: a_w, da_dm
    PetscReal :: poynting, dpy_dP_partial, dpy_dT_partial
    PetscReal :: V_w, RT
    PetscReal :: dyw_dPsat, dyw_dpy

    ierr  = 0
    y_w   = 0.0d0
    dyw_dP = 0.0d0; dyw_dT = 0.0d0; dyw_dm = 0.0d0

    if (P <= 0.0d0 .or. T <= 0.0d0 .or. m_NaCl < 0.0d0 .or. phi_w <= 0.0d0) then
      ierr = 1
      return
    end if

    call WaterSatPressure(T, P_sat, dPsat_dT, ierr)
    if (ierr /= 0) return

    call WaterActivityNaCl(m_NaCl, a_w, da_dm, ierr)
    if (ierr /= 0) return

    V_w = WaterMolarVolumeLiquid(T)
    call WaterPoyntingFactor(P, T, P_sat, V_w, poynting, &
                              dpy_dP_partial, dpy_dT_partial, ierr)
    if (ierr /= 0) return

    ! y_w = a_w * P_sat * Poynting / (phi_w * P)
    y_w = a_w * P_sat * poynting / (phi_w * P)

    ! Derivative w.r.t. P:
    !   y_w depends on P explicitly (through 1/P) and through Poynting.
    !   d(y_w)/dP|_T = a_w*P_sat*[d(Poynting)/dP / (phi_w*P)
    !                  - Poynting / (phi_w*P^2)]
    dyw_dP = a_w * P_sat / phi_w * &
             ( dpy_dP_partial / P - poynting / (P * P) )

    ! Derivative w.r.t. T:
    !   y_w depends on T through P_sat(T) and Poynting(P, T, P_sat).
    !   d(Poynting)/dT|_full = dpy_dT_partial - V_w/(RT)*dPsat_dT
    !   (the partial holds P_sat fixed; the chain through P_sat adds
    !    the second term: d/dT[V_w (P - P_sat)/RT] includes -V_w/(RT)*dPsat_dT)
    !
    !   y_w = a_w/(phi_w*P) * P_sat * poynting
    !   dyw_dT = a_w/(phi_w*P) * [dPsat_dT*poynting + P_sat*d(poynting)/dT_full]
    RT = IDEAL_GAS_CONSTANT * T
    dyw_dT = a_w / (phi_w * P) * &
             ( dPsat_dT * poynting + &
               P_sat * (dpy_dT_partial - V_w * dPsat_dT / RT * poynting) )

    ! Derivative w.r.t. m_NaCl:
    !   y_w depends on m only through a_w.
    dyw_dm = da_dm * P_sat * poynting / (phi_w * P)

  end subroutine WaterVapourMoleFraction

! =====================================================================
! Subroutine: WaterVapourVerifyDerivatives
!
! Compares analytic vs centered-FD derivatives for y_w(P, T, m).
! =====================================================================
  subroutine WaterVapourVerifyDerivatives(P, T, m_NaCl, phi_w, &
                                           max_rel_err, ierr)
    PetscReal, intent(in)       :: P, T, m_NaCl, phi_w
    PetscReal, intent(out)      :: max_rel_err
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: y0, dP_a, dT_a, dm_a
    PetscReal :: yp, ym, dPdum, dTdum, dmdum
    PetscReal :: dP_fd, dT_fd, dm_fd
    PetscReal :: dP, dT, dm, err
    PetscErrorCode :: ie

    ierr = 0
    max_rel_err = 0.0d0

    call WaterVapourMoleFraction(P, T, m_NaCl, phi_w, &
                                  y0, dP_a, dT_a, dm_a, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    dP = max(1.0d0, 1.0d-6 * P)
    dT = max(1.0d-4, 1.0d-6 * T)
    dm = max(1.0d-6, 1.0d-6 * max(m_NaCl, 1.0d0))

    call WaterVapourMoleFraction(P + dP, T, m_NaCl, phi_w, &
                                  yp, dPdum, dTdum, dmdum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    call WaterVapourMoleFraction(P - dP, T, m_NaCl, phi_w, &
                                  ym, dPdum, dTdum, dmdum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    dP_fd = (yp - ym) / (2.0d0 * dP)

    call WaterVapourMoleFraction(P, T + dT, m_NaCl, phi_w, &
                                  yp, dPdum, dTdum, dmdum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    call WaterVapourMoleFraction(P, T - dT, m_NaCl, phi_w, &
                                  ym, dPdum, dTdum, dmdum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    dT_fd = (yp - ym) / (2.0d0 * dT)

    if (m_NaCl >= dm) then
      call WaterVapourMoleFraction(P, T, m_NaCl + dm, phi_w, &
                                    yp, dPdum, dTdum, dmdum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      call WaterVapourMoleFraction(P, T, m_NaCl - dm, phi_w, &
                                    ym, dPdum, dTdum, dmdum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      dm_fd = (yp - ym) / (2.0d0 * dm)
    else
      call WaterVapourMoleFraction(P, T, m_NaCl + dm, phi_w, &
                                    yp, dPdum, dTdum, dmdum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      dm_fd = (yp - y0) / dm
    end if

    err = abs(dP_a - dP_fd) / max(abs(dP_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(dT_a - dT_fd) / max(abs(dT_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(dm_a - dm_fd) / max(abs(dm_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)

  end subroutine WaterVapourVerifyDerivatives

end module Air_WaterVapour_module
