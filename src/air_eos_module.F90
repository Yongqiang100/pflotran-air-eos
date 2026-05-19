! =====================================================================
! Air_EOS_module
!
! PFLOTRAN-facing wrapper that consolidates the four verified
! building-block modules into a single interface matching PFLOTRAN's
! GENERAL-mode auxiliary-structure conventions.
!
! This module is the *integration surface* for PFLOTRAN's GENERAL mode.
! It does no new physics — every call dispatches to one of the
! underlying modules and returns the result in the shape PFLOTRAN's
! Newton solver expects (property + analytic Jacobian terms in a
! single subroutine call).
!
! Dispatches:
!   AirEOSGetGasProperties      -> Air_EOS_PR_module + Air_Properties_module
!   AirEOSGetWaterVapourFraction-> Air_WaterVapour_module
!   AirEOSGetAirSolubility      -> Air_Henry_module
!   AirEOSGetTotalEnthalpy      -> Air_EOS_PR_module + Air_Properties_module
!
! PFLOTRAN GENERAL-mode integration points (see INTEGRATION.md for
! line-level patch instructions):
!
!   1. src/pflotran/general_aux.F90
!      Replace co2_span_wagner dispatcher calls with calls into this
!      module when a new "air_general" fluid option is set.
!
!   2. src/pflotran/general.F90
!      Add primary-variable switching cases for air-water (analogous
!      to existing CO2-water but with simpler phase-transition logic
!      since air does not undergo near-CAES-conditions phase change).
!
!   3. src/pflotran/input_aux.F90
!      Register "FLUID_TYPE AIR_GENERAL" keyword in the deck parser.
!
! Conventions matching PFLOTRAN GENERAL mode:
!   - Pressure in Pa, temperature in K
!   - Density in kg/m^3, enthalpy in J/mol (PFLOTRAN converts to J/kg
!     internally where needed)
!   - Mass fractions, not mole fractions, are PFLOTRAN's typical
!     primary variables; this wrapper accepts/returns mole fractions
!     and PFLOTRAN-side code converts via molar masses
!   - All derivatives are with respect to PFLOTRAN's primary variables
!     P and T at constant composition; chain rule through composition
!     handled by the calling GENERAL-mode auxiliary routines
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

module Air_EOS_module

#ifdef PFLOTRAN_INTEGRATION
#include "petsc/finclude/petscsys.h"
  use petscsys
#endif

  use Air_EOS_PR_module
  use Air_Henry_module
  use Air_WaterVapour_module
  use Air_Properties_module

  implicit none

  private

  ! ------------------------------------------------------------------
  ! Public interface — five top-level routines exposing all underlying
  ! functionality through a PFLOTRAN-style dispatcher pattern.
  ! ------------------------------------------------------------------
  public :: AirEOSGetGasProperties
  public :: AirEOSGetTotalEnthalpy
  public :: AirEOSGetWaterVapourFraction
  public :: AirEOSGetAirSolubility
  public :: AirEOSGetDissolvedO2
  public :: AirEOSVerify

  ! ------------------------------------------------------------------
  ! Reference state for total enthalpy. PFLOTRAN typically uses
  ! H(298.15 K, 1 atm, ideal gas) = 0 as the reference; we follow that.
  ! ------------------------------------------------------------------
  PetscReal, parameter :: T_REF_ENTH = 298.15d0   ! K
  PetscReal, parameter :: AIR_MW_KG = 28.9586d-3  ! kg/mol (PFLOTRAN uses SI)

contains

! =====================================================================
! Subroutine: AirEOSGetGasProperties
!
! Returns the full set of gas-phase thermodynamic and transport
! properties at (P, T), with derivatives needed for Jacobian assembly.
! This is the single most-frequently-called routine by PFLOTRAN's
! GENERAL-mode auxiliary update.
!
! Output convention matches PFLOTRAN: all derivatives w.r.t. P and T
! at constant composition (PFLOTRAN-side code handles chain rule
! through composition variables).
! =====================================================================
  subroutine AirEOSGetGasProperties(P, T, &
                                     rho, mu, fug_coeff, h_dep, &
                                     drho_dP, drho_dT, &
                                     dmu_dP, dmu_dT, &
                                     dfug_dP, dfug_dT, &
                                     dhdep_dP, dhdep_dT, &
                                     ierr)
    PetscReal, intent(in)       :: P, T
    PetscReal, intent(out)      :: rho, mu, fug_coeff, h_dep
    PetscReal, intent(out)      :: drho_dP, drho_dT
    PetscReal, intent(out)      :: dmu_dP, dmu_dT
    PetscReal, intent(out)      :: dfug_dP, dfug_dT
    PetscReal, intent(out)      :: dhdep_dP, dhdep_dT
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: Z_pr
    PetscReal :: dmu_drho

    ierr = 0
    rho = 0.0d0; mu = 0.0d0; fug_coeff = 0.0d0; h_dep = 0.0d0
    drho_dP = 0.0d0; drho_dT = 0.0d0
    dmu_dP  = 0.0d0; dmu_dT  = 0.0d0
    dfug_dP = 0.0d0; dfug_dT = 0.0d0
    dhdep_dP = 0.0d0; dhdep_dT = 0.0d0

    if (P <= 0.0d0 .or. T <= 0.0d0) then
      ierr = 1
      return
    end if

    ! ----- PR + Peneloux: density, fugacity coefficient, enthalpy departure -----
    call AirEOSPRProperties(P, T, rho, Z_pr, fug_coeff, h_dep, &
                            drho_dP, drho_dT, dfug_dP, dfug_dT, ierr)
    if (ierr /= 0) return

    ! Enthalpy departure derivatives are not directly returned by
    ! AirEOSPRProperties (only h_dep itself). Compute via FD on the
    ! existing fugacity-departure routine — these are small numbers
    ! and PFLOTRAN's Newton solver tolerates ~10^-6 FD accuracy here.
    ! TODO: add analytic h_dep derivatives to AirEOSPRProperties
    !       when integrating in PFLOTRAN production runs (deliverable 1b).
    call ComputeEnthDepDerivativesFD(P, T, dhdep_dP, dhdep_dT, ierr)
    if (ierr /= 0) return

    ! ----- L-J 2004: viscosity -----
    call AirViscosity(T, rho, mu, dmu_dT, dmu_drho, ierr)
    if (ierr /= 0) return

    ! Chain rule: dmu/dP|_T = dmu/drho|_T * drho/dP|_T
    !             dmu/dT|_P = dmu/dT|_rho + dmu/drho|_T * drho/dT|_P
    dmu_dP = dmu_drho * drho_dP
    dmu_dT = dmu_dT + dmu_drho * drho_dT

  end subroutine AirEOSGetGasProperties

! =====================================================================
! Subroutine: AirEOSGetTotalEnthalpy
!
! Returns the total molar enthalpy h(P, T) - h_ref(T_REF) [J/mol]
! by combining the ideal-gas reference enthalpy and the PR enthalpy
! departure. The reference state is (T_REF, ideal-gas) so
! h(T_REF, 1 atm, ideal) = 0.
!
! For PFLOTRAN: this is the gas-phase total enthalpy needed in the
! energy balance equation.
! =====================================================================
  subroutine AirEOSGetTotalEnthalpy(P, T, h_total, dh_dP, dh_dT, ierr)
    PetscReal, intent(in)       :: P, T
    PetscReal, intent(out)      :: h_total, dh_dP, dh_dT
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: h_ig, dhig_dT, h_dep, dhdep_dP, dhdep_dT
    PetscReal :: rho_dum, mu_dum, fug_dum
    PetscReal :: drho_dP_dum, drho_dT_dum, dmu_dP_dum, dmu_dT_dum
    PetscReal :: dfug_dP_dum, dfug_dT_dum

    ierr = 0
    h_total = 0.0d0; dh_dP = 0.0d0; dh_dT = 0.0d0

    ! Ideal-gas contribution
    call AirIdealGasEnthalpy(T, h_ig, dhig_dT, ierr)
    if (ierr /= 0) return

    ! Departure contribution
    call AirEOSGetGasProperties(P, T, &
                                 rho_dum, mu_dum, fug_dum, h_dep, &
                                 drho_dP_dum, drho_dT_dum, &
                                 dmu_dP_dum, dmu_dT_dum, &
                                 dfug_dP_dum, dfug_dT_dum, &
                                 dhdep_dP, dhdep_dT, ierr)
    if (ierr /= 0) return

    h_total = h_ig + h_dep
    dh_dP   = dhdep_dP
    dh_dT   = dhig_dT + dhdep_dT

  end subroutine AirEOSGetTotalEnthalpy

! =====================================================================
! Subroutine: AirEOSGetWaterVapourFraction
!
! Returns the equilibrium water vapour mole fraction y_w in the gas
! phase, given total pressure, temperature, and brine salinity.
!
! Uses the bulk-gas fugacity coefficient from PR+Peneloux as a proxy
! for phi_w (the dilute-water-vapour approximation valid for CAES
! conditions where y_w < 1 % at T < 60 C).
! =====================================================================
  subroutine AirEOSGetWaterVapourFraction(P, T, m_NaCl, &
                                           y_w, dyw_dP, dyw_dT, dyw_dm, &
                                           ierr)
    PetscReal, intent(in)       :: P, T, m_NaCl
    PetscReal, intent(out)      :: y_w, dyw_dP, dyw_dT, dyw_dm
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: phi_air
    PetscErrorCode :: ie

    ierr = 0
    y_w = 0.0d0; dyw_dP = 0.0d0; dyw_dT = 0.0d0; dyw_dm = 0.0d0

    ! Bulk gas fugacity coefficient
    phi_air = AirEOSPRFugacityCoeff(P, T, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    ! Water vapour mole fraction with full derivatives
    call WaterVapourMoleFraction(P, T, m_NaCl, phi_air, &
                                  y_w, dyw_dP, dyw_dT, dyw_dm, ierr)

  end subroutine AirEOSGetWaterVapourFraction

! =====================================================================
! Subroutine: AirEOSGetAirSolubility
!
! Returns the equilibrium aqueous mole fractions of O2, N2, CO2, Ar
! given gas-phase composition and brine salinity. Uses the PR+Peneloux
! fugacity coefficient for the high-pressure correction.
! =====================================================================
  subroutine AirEOSGetAirSolubility(P, T, m_NaCl, y_gas, &
                                     x_aq, dxaq_dP, dxaq_dT, dxaq_dm, &
                                     ierr)
    PetscReal, intent(in)       :: P, T, m_NaCl
    PetscReal, intent(in)       :: y_gas(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: x_aq(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: dxaq_dP(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: dxaq_dT(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: dxaq_dm(N_HENRY_SPECIES)
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: phi_air
    PetscReal :: phi_gas(N_HENRY_SPECIES)
    PetscErrorCode :: ie

    ierr = 0
    x_aq = 0.0d0
    dxaq_dP = 0.0d0; dxaq_dT = 0.0d0; dxaq_dm = 0.0d0

    ! Use bulk fugacity coefficient for all species. A rigorous
    ! multicomponent treatment would compute species-specific
    ! fugacity coefficients from a multicomponent EOS; this is
    ! deferred to a future enhancement (deliverable 1b).
    phi_air = AirEOSPRFugacityCoeff(P, T, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    phi_gas = phi_air

    call AirHenryEquilibrium(P, T, m_NaCl, y_gas, phi_gas, &
                              x_aq, dxaq_dP, dxaq_dT, dxaq_dm, ierr)

  end subroutine AirEOSGetAirSolubility

! =====================================================================
! Subroutine: AirEOSGetDissolvedO2
!
! Convenience wrapper returning only dissolved O2 concentration with
! derivatives. The most-frequently-called solubility routine in the
! chemistry coupling (O2 drives the pyrite oxidation kinetics).
! =====================================================================
  subroutine AirEOSGetDissolvedO2(P, T, m_NaCl, y_O2, &
                                   c_O2_aq, dc_dP, dc_dT, dc_dm, ierr)
    PetscReal, intent(in)       :: P, T, m_NaCl, y_O2
    PetscReal, intent(out)      :: c_O2_aq, dc_dP, dc_dT, dc_dm
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: phi_air
    PetscErrorCode :: ie

    ierr = 0

    phi_air = AirEOSPRFugacityCoeff(P, T, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    call AirHenryDissolvedO2(P, T, m_NaCl, y_O2, phi_air, &
                              c_O2_aq, dc_dP, dc_dT, dc_dm, ierr)

  end subroutine AirEOSGetDissolvedO2

! =====================================================================
! Subroutine: AirEOSVerify
!
! Self-test routine called by PFLOTRAN at startup to confirm the
! wrapper produces sensible numbers. Returns 0 on success.
! =====================================================================
  subroutine AirEOSVerify(verbose, ierr)
    logical, intent(in)         :: verbose
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: P, T, m
    PetscReal :: rho, mu, fug, h_dep, h_tot
    PetscReal :: drho_dP, drho_dT, dmu_dP, dmu_dT
    PetscReal :: dfug_dP, dfug_dT, dhdep_dP, dhdep_dT
    PetscReal :: dh_dP, dh_dT
    PetscReal :: y_w, dyw_dP, dyw_dT, dyw_dm
    PetscReal :: c_O2, dc_dP, dc_dT, dc_dm

    ierr = 0

    ! Representative CAES state
    P = 8.0d6
    T = 313.15d0
    m = 0.5d0

    call AirEOSGetGasProperties(P, T, &
                                 rho, mu, fug, h_dep, &
                                 drho_dP, drho_dT, dmu_dP, dmu_dT, &
                                 dfug_dP, dfug_dT, dhdep_dP, dhdep_dT, ierr)
    if (ierr /= 0) return

    call AirEOSGetTotalEnthalpy(P, T, h_tot, dh_dP, dh_dT, ierr)
    if (ierr /= 0) return

    call AirEOSGetWaterVapourFraction(P, T, m, y_w, dyw_dP, dyw_dT, dyw_dm, ierr)
    if (ierr /= 0) return

    call AirEOSGetDissolvedO2(P, T, m, 0.20946d0, c_O2, dc_dP, dc_dT, dc_dm, ierr)
    if (ierr /= 0) return

    ! Sanity-check the bounds of all returned values
    if (rho < 1.0d0 .or. rho > 1000.0d0) ierr = ierr + 10
    if (mu  < 1.0d-6 .or. mu > 1.0d-3)  ierr = ierr + 20
    if (fug < 0.5d0 .or. fug > 2.0d0)   ierr = ierr + 30
    if (y_w < 0.0d0 .or. y_w > 0.1d0)   ierr = ierr + 40
    if (c_O2 < 0.0d0)                    ierr = ierr + 50

    if (verbose) then
      print '(a)',     '  Air_EOS_module self-check at 80 bar, 40 C, 0.5 M NaCl:'
      print '(a,f7.2,a)', '    rho   = ', rho,    ' kg/m^3'
      print '(a,es10.3,a)', '    mu    = ', mu,    ' Pa.s'
      print '(a,f7.4)',   '    fug   = ', fug
      print '(a,f7.2,a)', '    h_tot = ', h_tot,  ' J/mol (relative to ideal-gas at 298.15 K)'
      print '(a,es10.3)', '    y_w   = ', y_w
      print '(a,es10.3,a)', '    c_O2  = ', c_O2, ' mol/m^3'
      print '(a,i0)',     '    ierr  = ', ierr
    end if

  end subroutine AirEOSVerify

! =====================================================================
! Private helper: finite-difference enthalpy departure derivatives
! Placeholder until analytic dh_dep/d{P,T} added to PR module.
! =====================================================================
  subroutine ComputeEnthDepDerivativesFD(P, T, dhdep_dP, dhdep_dT, ierr)
    PetscReal, intent(in)       :: P, T
    PetscReal, intent(out)      :: dhdep_dP, dhdep_dT
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: dP, dT, hp, hm
    PetscErrorCode :: ie

    ierr = 0
    dhdep_dP = 0.0d0; dhdep_dT = 0.0d0

    dP = max(1.0d0, 1.0d-6 * P)
    dT = max(1.0d-3, 1.0d-6 * T)

    hp = AirEOSPREnthalpyDeparture(P + dP, T, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    hm = AirEOSPREnthalpyDeparture(P - dP, T, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    dhdep_dP = (hp - hm) / (2.0d0 * dP)

    hp = AirEOSPREnthalpyDeparture(P, T + dT, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    hm = AirEOSPREnthalpyDeparture(P, T - dT, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    dhdep_dT = (hp - hm) / (2.0d0 * dT)

  end subroutine ComputeEnthDepDerivativesFD

end module Air_EOS_module
