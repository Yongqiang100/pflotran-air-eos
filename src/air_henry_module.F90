! =====================================================================
! Air_Henry_module
!
! Henry's-law solubility module for the four atmospheric components
! O2, N2, CO2, Ar in saline brines. Provides per-species Henry's
! constants H^cp(T, salinity) and the equilibrium aqueous mole
! fractions x_aq given a gas-phase composition y_gas and the gas-
! phase fugacity coefficients phi_i from Air_EOS_PR_module.
!
! Temperature dependence: Van't Hoff form
!   H(T) = H_ref * exp[ B * (1/T - 1/T_ref) ]
! where B = d(ln H)/d(1/T) is tabulated per species (Sander 2015).
!
! Salinity dependence: Schumpe (1993) / Weisenberger & Schumpe (1996)
! reduced to an effective Setschenow constant for NaCl brines:
!   log10(H(T, m) / H(T, 0)) = -K_S * m_NaCl
! where m_NaCl is molality of NaCl in mol/kg-H2O and K_S is the
! species-specific Setschenow coefficient. This formulation is
! adequate for CAES brine compositions; for more complex multi-ion
! brines, the full Schumpe sum should be used instead.
!
! Coupling to gas-phase fugacity: at CAES pressures (50-100 bar),
! fugacity corrections matter, so the dissolved mole fraction is
!   x_aq_i = (y_i * phi_i * P * H^cp_i) / rho_w_molar
! where rho_w_molar ~ 55509 mol/m^3 is the molar density of liquid
! water (approximate; held constant here, refined when coupled to
! the IAPWS water EOS in deliverable 3).
!
! References:
!   Sander, R., 2015. Compilation of Henry's law constants (v4.0)
!     for water as solvent. Atmos. Chem. Phys. 15, 4399-4981.
!   Schumpe, A., 1993. The estimation of gas solubilities in salt
!     solutions. Chem. Eng. Sci. 48, 153-158.
!   Weisenberger, S., Schumpe, A., 1996. Estimation of gas
!     solubilities in salt solutions at temperatures from 273 K
!     to 363 K. AIChE J. 42, 298-300.
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

module Air_Henry_module

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
  public :: AirHenryConstant
  public :: AirHenryEquilibrium
  public :: AirHenryDissolvedO2
  public :: AirHenryVerifyDerivatives

  ! ------------------------------------------------------------------
  ! Species enumeration (public, used as indices into parameter arrays)
  ! ------------------------------------------------------------------
  PetscInt, parameter, public :: HENRY_O2  = 1
  PetscInt, parameter, public :: HENRY_N2  = 2
  PetscInt, parameter, public :: HENRY_CO2 = 3
  PetscInt, parameter, public :: HENRY_AR  = 4
  PetscInt, parameter, public :: N_HENRY_SPECIES = 4

  ! ------------------------------------------------------------------
  ! Reference temperature for Van't Hoff
  ! ------------------------------------------------------------------
  PetscReal, parameter :: HENRY_T_REF = 298.15d0   ! K

  ! ------------------------------------------------------------------
  ! Henry's-law constants at T_ref = 298.15 K in pure water [mol/(m^3 Pa)]
  ! From Sander (2015) compilation, geometric mean of reported values.
  !
  ! Index order: O2, N2, CO2, Ar (HENRY_* constants above)
  ! ------------------------------------------------------------------
  PetscReal, parameter, dimension(N_HENRY_SPECIES) :: HENRY_H_REF = &
    [ 1.3d-5,   &   ! O2:  ~ 1.3e-5 mol/(m^3 Pa)
      6.4d-6,   &   ! N2:  ~ 6.4e-6 mol/(m^3 Pa)
      3.3d-4,   &   ! CO2: ~ 3.3e-4 mol/(m^3 Pa)
      1.4d-5 ]      ! Ar:  ~ 1.4e-5 mol/(m^3 Pa)

  ! ------------------------------------------------------------------
  ! Van't Hoff coefficients B = d(ln H)/d(1/T) [K]
  ! Sander (2015) tabulated values; H(T) = H_ref * exp[B (1/T - 1/T_ref)]
  ! ------------------------------------------------------------------
  PetscReal, parameter, dimension(N_HENRY_SPECIES) :: HENRY_VANT_HOFF = &
    [ 1500.0d0,  &  ! O2
      1300.0d0,  &  ! N2
      2400.0d0,  &  ! CO2
      1500.0d0 ]    ! Ar

  ! ------------------------------------------------------------------
  ! Effective Setschenow coefficients for NaCl brines [L/mol]
  ! log10(H(T,m)/H(T,0)) = -K_S * m_NaCl  (salting-out: H decreases,
  ! i.e. less gas dissolves per unit partial pressure)
  ! From Weisenberger & Schumpe (1996) and compiled literature.
  ! ------------------------------------------------------------------
  PetscReal, parameter, dimension(N_HENRY_SPECIES) :: HENRY_SETSCHENOW = &
    [ 0.143d0,   &  ! O2
      0.135d0,   &  ! N2
      0.119d0,   &  ! CO2
      0.142d0 ]     ! Ar

  ! ------------------------------------------------------------------
  ! Molar density of liquid water at standard conditions [mol/m^3]
  ! Used for converting concentration to mole fraction. Will be made
  ! variable in deliverable 3 when coupled to a water EOS; for now
  ! held constant since the variation is small over the CAES T range.
  ! ------------------------------------------------------------------
  PetscReal, parameter :: WATER_MOLAR_DENSITY = 55509.0d0  ! mol/m^3

  ! ------------------------------------------------------------------
  ! Natural logarithm of 10, for Setschenow base-10 -> base-e conversion
  ! ------------------------------------------------------------------
  PetscReal, parameter :: LN_10 = 2.30258509299404568402d0

  ! ------------------------------------------------------------------
  ! Gas constant fallback for standalone build
  ! ------------------------------------------------------------------
#ifndef PFLOTRAN_INTEGRATION
  PetscReal, parameter :: IDEAL_GAS_CONSTANT = 8.31446261815324d0  ! J/(mol K)
#endif

contains

! =====================================================================
! Subroutine: AirHenryConstant
!
! Returns Henry's-law constant H^cp(T, m_NaCl) for one species in
! mol/(m^3 Pa), with analytic derivatives w.r.t. T and salinity.
!
! Input:
!   species_id  [-]            One of HENRY_O2, HENRY_N2, HENRY_CO2, HENRY_AR
!   T           [K]            Temperature
!   m_NaCl      [mol/kg-H2O]   NaCl molality (use 0 for pure water)
!
! Output:
!   H_cp        [mol/(m^3 Pa)] Henry's constant
!   dH_dT       [1/(m^3 Pa K)] dH^cp/dT at constant salinity
!   dH_dm       [kg/(m^3 Pa)]  dH^cp/dm at constant T
!   ierr        [-]            0 success, 1 bad species id, 2 bad T or m
! =====================================================================
  subroutine AirHenryConstant(species_id, T, m_NaCl, H_cp, dH_dT, dH_dm, ierr)
    PetscInt,  intent(in)       :: species_id
    PetscReal, intent(in)       :: T, m_NaCl
    PetscReal, intent(out)      :: H_cp, dH_dT, dH_dm
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: H_ref, B, K_S
    PetscReal :: vh_factor, salt_factor

    ierr  = 0
    H_cp  = 0.0d0
    dH_dT = 0.0d0
    dH_dm = 0.0d0

    if (species_id < 1 .or. species_id > N_HENRY_SPECIES) then
      ierr = 1
      return
    end if
    if (T <= 0.0d0 .or. m_NaCl < 0.0d0) then
      ierr = 2
      return
    end if

    H_ref = HENRY_H_REF(species_id)
    B     = HENRY_VANT_HOFF(species_id)
    K_S   = HENRY_SETSCHENOW(species_id)

    ! Van't Hoff temperature factor: exp[B(1/T - 1/T_ref)]
    vh_factor = exp(B * (1.0d0/T - 1.0d0/HENRY_T_REF))

    ! Setschenow salinity factor: 10^(-K_S * m) = exp(-ln10 * K_S * m)
    salt_factor = exp(-LN_10 * K_S * m_NaCl)

    H_cp = H_ref * vh_factor * salt_factor

    ! d(H)/dT |_m :
    !   d/dT [exp(B(1/T - 1/T_ref))] = exp(...) * (-B/T^2)
    !   d/dT [salt_factor] = 0
    !   so dH/dT = H * (-B/T^2)
    dH_dT = H_cp * (-B / (T * T))

    ! d(H)/dm |_T :
    !   d/dm [exp(-ln10 * K_S * m)] = exp(...) * (-ln10 * K_S)
    !   so dH/dm = H * (-ln10 * K_S)
    dH_dm = H_cp * (-LN_10 * K_S)

  end subroutine AirHenryConstant

! =====================================================================
! Subroutine: AirHenryEquilibrium
!
! Returns the equilibrium aqueous mole fractions for all four
! components given gas-phase composition y_gas, total pressure P,
! gas-phase fugacity coefficients phi_gas (typically from
! Air_EOS_PR_module), temperature T, and brine salinity m_NaCl.
!
! Working equation:
!   c_aq_i = y_i * phi_i * P * H^cp_i   [mol/m^3]
!   x_aq_i = c_aq_i / rho_w_molar
!
! Input:
!   P             [Pa]           Total gas-phase pressure
!   T             [K]            Temperature
!   m_NaCl        [mol/kg-H2O]   Brine salinity (NaCl)
!   y_gas(:)      [-]            Gas-phase mole fractions
!                                (size N_HENRY_SPECIES; should sum to 1)
!   phi_gas(:)    [-]            Gas-phase fugacity coefficients
!                                (size N_HENRY_SPECIES; pass 1.0 for ideal gas)
!
! Output:
!   x_aq(:)       [-]            Equilibrium aqueous mole fractions
!   dxaq_dP(:)    [1/Pa]
!   dxaq_dT(:)    [1/K]
!   dxaq_dm(:)    [kg/mol]
!   ierr          [-]            Error code
!
! Notes:
!   * Fugacity coefficient is treated as independent input (not
!     differentiated through here); the caller is responsible for
!     using phi(P, T) from the EOS module and applying the full
!     chain rule if needed.
!   * Assumes dilute aqueous concentrations (Henry's-law regime).
!     Not valid above ~ 1 % aqueous mole fraction; for CAES O2
!     conditions this assumption holds with wide margin.
! =====================================================================
  subroutine AirHenryEquilibrium(P, T, m_NaCl, y_gas, phi_gas, x_aq, &
                                  dxaq_dP, dxaq_dT, dxaq_dm, ierr)
    PetscReal, intent(in)       :: P, T, m_NaCl
    PetscReal, intent(in)       :: y_gas(N_HENRY_SPECIES)
    PetscReal, intent(in)       :: phi_gas(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: x_aq(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: dxaq_dP(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: dxaq_dT(N_HENRY_SPECIES)
    PetscReal, intent(out)      :: dxaq_dm(N_HENRY_SPECIES)
    PetscErrorCode, intent(out) :: ierr

    PetscReal      :: H_cp, dH_dT, dH_dm
    PetscReal      :: c_aq, fugacity
    PetscInt       :: i
    PetscErrorCode :: ie

    ierr = 0
    x_aq    = 0.0d0
    dxaq_dP = 0.0d0
    dxaq_dT = 0.0d0
    dxaq_dm = 0.0d0

    if (P <= 0.0d0 .or. T <= 0.0d0 .or. m_NaCl < 0.0d0) then
      ierr = 1
      return
    end if

    do i = 1, N_HENRY_SPECIES
      call AirHenryConstant(i, T, m_NaCl, H_cp, dH_dT, dH_dm, ie)
      if (ie /= 0) then
        ierr = ie
        return
      end if

      fugacity = y_gas(i) * phi_gas(i) * P
      c_aq     = fugacity * H_cp
      x_aq(i)  = c_aq / WATER_MOLAR_DENSITY

      ! Derivatives. f = y*phi*P. We treat y and phi as fixed inputs
      ! to this routine; H depends on T and m only.
      !   dxaq/dP  = y * phi * H / rho_w
      !   dxaq/dT  = y * phi * P * (dH/dT) / rho_w
      !   dxaq/dm  = y * phi * P * (dH/dm) / rho_w
      dxaq_dP(i) = y_gas(i) * phi_gas(i) * H_cp / WATER_MOLAR_DENSITY
      dxaq_dT(i) = y_gas(i) * phi_gas(i) * P * dH_dT / WATER_MOLAR_DENSITY
      dxaq_dm(i) = y_gas(i) * phi_gas(i) * P * dH_dm / WATER_MOLAR_DENSITY
    end do

  end subroutine AirHenryEquilibrium

! =====================================================================
! Subroutine: AirHenryDissolvedO2
!
! Convenience wrapper returning only dissolved O2 properties. Useful
! for chemistry coupling, since O2 is the primary reactive species
! driving pyrite oxidation in the CAES context.
!
! Input:
!   P       [Pa]           Total gas-phase pressure
!   T       [K]            Temperature
!   m_NaCl  [mol/kg-H2O]   Salinity
!   y_O2    [-]            Gas-phase O2 mole fraction (use 0.21 for air)
!   phi_O2  [-]            Gas-phase O2 fugacity coefficient
!
! Output:
!   c_O2_aq  [mol/m^3]     Dissolved O2 concentration
!   dc_dP    [mol/(m^3 Pa)]
!   dc_dT    [mol/(m^3 K)]
!   dc_dm    [kg/m^3]
!   ierr     [-]
! =====================================================================
  subroutine AirHenryDissolvedO2(P, T, m_NaCl, y_O2, phi_O2, &
                                  c_O2_aq, dc_dP, dc_dT, dc_dm, ierr)
    PetscReal, intent(in)       :: P, T, m_NaCl, y_O2, phi_O2
    PetscReal, intent(out)      :: c_O2_aq, dc_dP, dc_dT, dc_dm
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: H_cp, dH_dT, dH_dm
    PetscReal :: fug

    ierr = 0
    c_O2_aq = 0.0d0; dc_dP = 0.0d0; dc_dT = 0.0d0; dc_dm = 0.0d0

    if (P <= 0.0d0 .or. T <= 0.0d0 .or. m_NaCl < 0.0d0) then
      ierr = 1
      return
    end if

    call AirHenryConstant(HENRY_O2, T, m_NaCl, H_cp, dH_dT, dH_dm, ierr)
    if (ierr /= 0) return

    fug      = y_O2 * phi_O2 * P
    c_O2_aq  = fug * H_cp
    dc_dP    = y_O2 * phi_O2 * H_cp
    dc_dT    = fug * dH_dT
    dc_dm    = fug * dH_dm

  end subroutine AirHenryDissolvedO2

! =====================================================================
! Subroutine: AirHenryVerifyDerivatives
!
! Compares analytic vs centered-FD derivatives of H^cp w.r.t. T and m
! at a chosen test point. Diagnostic only.
! =====================================================================
  subroutine AirHenryVerifyDerivatives(species_id, T, m_NaCl, &
                                        max_rel_err, ierr)
    PetscInt,  intent(in)       :: species_id
    PetscReal, intent(in)       :: T, m_NaCl
    PetscReal, intent(out)      :: max_rel_err
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: H0, dH_dT_a, dH_dm_a
    PetscReal :: Hp, Hm, dHd, dHd_dum
    PetscReal :: dH_dT_fd, dH_dm_fd
    PetscReal :: dT, dm, err
    PetscErrorCode :: ie

    ierr = 0
    max_rel_err = 0.0d0

    call AirHenryConstant(species_id, T, m_NaCl, H0, dH_dT_a, dH_dm_a, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    dT = max(1.0d-4, 1.0d-6 * T)
    dm = max(1.0d-6, 1.0d-6 * max(m_NaCl, 1.0d0))

    call AirHenryConstant(species_id, T + dT, m_NaCl, Hp, dHd, dHd_dum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    call AirHenryConstant(species_id, T - dT, m_NaCl, Hm, dHd, dHd_dum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    dH_dT_fd = (Hp - Hm) / (2.0d0 * dT)

    ! For salinity FD: use centered when away from boundary, forward
    ! when m_NaCl is too close to zero (m - dm would go negative).
    if (m_NaCl >= dm) then
      call AirHenryConstant(species_id, T, m_NaCl + dm, Hp, dHd, dHd_dum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      call AirHenryConstant(species_id, T, m_NaCl - dm, Hm, dHd, dHd_dum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      dH_dm_fd = (Hp - Hm) / (2.0d0 * dm)
    else
      ! Forward FD at the boundary m_NaCl ~ 0
      call AirHenryConstant(species_id, T, m_NaCl + dm, Hp, dHd, dHd_dum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      call AirHenryConstant(species_id, T, m_NaCl, H0, dHd, dHd_dum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      dH_dm_fd = (Hp - H0) / dm
    end if

    err = abs(dH_dT_a - dH_dT_fd) / max(abs(dH_dT_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(dH_dm_a - dH_dm_fd) / max(abs(dH_dm_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)

  end subroutine AirHenryVerifyDerivatives

end module Air_Henry_module
