! =====================================================================
! Air_Properties_module
!
! Auxiliary thermodynamic and transport properties for air, completing
! the set needed for the energy balance and the Darcy flux:
!
!   - Air viscosity mu(T, rho) via the Lemmon-Jacobsen (2004) reference
!     correlation: Chapman-Enskog dilute-gas viscosity plus a five-term
!     residual polynomial in (T, rho). Accuracy < 1% across the CAES
!     P-T regime (0-200 bar, 0-100 C).
!
!   - Air ideal-gas heat capacity cp_ig(T) and reference enthalpy
!     h_ig(T) - h_ig(T_ref) via the NIST Shomate polynomial (nitrogen
!     parameters, which match air cp to ~ 1 % across our T range).
!     Combined with the Peng-Robinson enthalpy departure from
!     Air_EOS_PR_module, this gives the total specific enthalpy
!     h(P, T) - h_ref needed for the energy conservation equation.
!
! References:
!   Lemmon, E.W., Jacobsen, R.T., 2004. Viscosity and thermal
!     conductivity equations for nitrogen, oxygen, argon, and air.
!     Int. J. Thermophysics 25, 21-69.
!   NIST Chemistry WebBook, Shomate coefficients for N2 (100-500 K
!     and 500-2000 K ranges).
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

module Air_Properties_module

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
  public :: AirViscosity
  public :: AirIdealGasEnthalpy
  public :: AirIdealGasHeatCapacity
  public :: AirViscosityVerifyDerivatives
  public :: AirEnthalpyVerifyDerivatives

  ! ------------------------------------------------------------------
  ! Air parameters used by both viscosity and ideal-gas EOS
  ! ------------------------------------------------------------------
  PetscReal, parameter :: AIR_MW_GMOL = 28.9586d0       ! g/mol
  PetscReal, parameter :: AIR_TC_VISC = 132.6312d0      ! K (L-J 2004 value)
  PetscReal, parameter :: AIR_RHOC_VISC = 302.5507d0    ! kg/m^3 (= 10.4477 mol/L * MW)
  PetscReal, parameter :: AIR_T_REF_ENTH = 298.15d0     ! K, enthalpy reference

  ! ------------------------------------------------------------------
  ! Chapman-Enskog parameters for air (Lemmon-Jacobsen 2004)
  ! sigma in nm; epsilon/k in K
  ! ------------------------------------------------------------------
  PetscReal, parameter :: AIR_SIGMA_LJ = 0.360d0        ! nm (collision diameter)
  PetscReal, parameter :: AIR_EPSILON_K = 103.3d0       ! K (energy parameter)

  ! Pre-factor in eta^0(T) = K * sqrt(M * T) / (sigma^2 * Omega(T*))
  ! returning eta^0 in micro-Pa.s
  PetscReal, parameter :: ETA0_PREFACTOR = 0.0266958d0

  ! Bich-Buchholz collision-integral fit coefficients for nitrogen/air
  ! ln(Omega) = b0 + b1*s + b2*s^2 + b3*s^3 + b4*s^4, where s = ln(T*)
  PetscReal, parameter :: B_OMEGA(0:4) =                                  &
    [  0.431d0,                                                           &
      -0.4623d0,                                                          &
       0.08406d0,                                                         &
       0.005341d0,                                                        &
      -0.00331d0 ]

  ! ------------------------------------------------------------------
  ! Residual viscosity polynomial coefficients (L-J 2004 Table 2 for air)
  ! eta^r(tau, delta) = sum_i N_i * tau^t_i * delta^d_i * exp(-gamma_i * delta^l_i)
  !   tau = Tc/T,  delta = rho/rho_c
  ! Result in micro-Pa.s
  ! ------------------------------------------------------------------
  PetscInt, parameter :: N_VISC_TERMS = 5
  PetscReal, parameter, dimension(N_VISC_TERMS) :: N_VISC =               &
    [ 10.72d0,    1.122d0,    0.002019d0, -8.876d0,    -0.02916d0 ]
  PetscReal, parameter, dimension(N_VISC_TERMS) :: T_VISC =               &
    [  0.2d0,     0.05d0,     2.4d0,       0.6d0,       3.6d0 ]
  PetscReal, parameter, dimension(N_VISC_TERMS) :: D_VISC =               &
    [  1.0d0,     4.0d0,      9.0d0,       1.0d0,       8.0d0 ]
  PetscReal, parameter, dimension(N_VISC_TERMS) :: L_VISC =               &
    [  0.0d0,     0.0d0,      0.0d0,       1.0d0,       1.0d0 ]
  PetscReal, parameter, dimension(N_VISC_TERMS) :: G_VISC =               &
    [  0.0d0,     0.0d0,      0.0d0,       1.0d0,       1.0d0 ]

  ! ------------------------------------------------------------------
  ! Shomate polynomial coefficients for N2, used as proxy for air.
  ! cp/[J/(mol*K)] = A + B*t + C*t^2 + D*t^3 + E/t^2,   t = T/1000
  ! H(T) - H(T_ref)/[kJ/mol] = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F - H_offset
  ! Values from NIST WebBook (100-500 K range).
  ! ------------------------------------------------------------------
  PetscReal, parameter :: SHOMATE_A      =  28.98641d0
  PetscReal, parameter :: SHOMATE_B      =   1.853978d0
  PetscReal, parameter :: SHOMATE_C      =  -9.647459d0
  PetscReal, parameter :: SHOMATE_D      =  16.63537d0
  PetscReal, parameter :: SHOMATE_E      =   0.000117d0
  PetscReal, parameter :: SHOMATE_F      =  -8.671914d0
  PetscReal, parameter :: SHOMATE_H_OFF  =   0.0d0

  ! Gas constant fallback
#ifndef PFLOTRAN_INTEGRATION
  PetscReal, parameter :: IDEAL_GAS_CONSTANT = 8.31446261815324d0  ! J/(mol K)
#endif

contains

! =====================================================================
! Subroutine: AirViscosity
!
! Returns the dynamic viscosity of air mu(T, rho) [Pa.s] via the
! Lemmon-Jacobsen (2004) correlation: eta = eta^0(T) + eta^r(T, rho).
!
! Input:
!   T          [K]       Temperature (recommended range 70-2000 K)
!   rho        [kg/m^3]  Air density
!
! Output:
!   mu         [Pa.s]    Viscosity
!   dmu_dT     [Pa.s/K]
!   dmu_drho   [Pa.s.m^3/kg]
!   ierr       [-]       Error code (0 = success)
! =====================================================================
  subroutine AirViscosity(T, rho, mu, dmu_dT, dmu_drho, ierr)
    PetscReal, intent(in)       :: T, rho
    PetscReal, intent(out)      :: mu, dmu_dT, dmu_drho
    PetscErrorCode, intent(out) :: ierr

    PetscReal :: T_star, lnTs, omega, dlnOmega_dlnTs
    PetscReal :: eta0_uPa, deta0_dT
    PetscReal :: tau, delta, eta_r_uPa, deta_r_dT, deta_r_drho
    PetscReal :: term, exp_arg, dterm_dT_factor, dterm_drho_factor
    PetscReal :: M_T_sqrt, M_T_sqrt_dT
    PetscInt  :: i

    ierr = 0
    mu = 0.0d0; dmu_dT = 0.0d0; dmu_drho = 0.0d0

    if (T <= 0.0d0 .or. rho < 0.0d0) then
      ierr = 1
      return
    end if

    ! ----- Dilute-gas viscosity eta^0(T) -----
    T_star = T / AIR_EPSILON_K
    lnTs   = log(T_star)

    ! Collision integral Omega(T*) = exp(sum b_i (ln T*)^i)
    omega = exp( B_OMEGA(0)                                    &
               + B_OMEGA(1) * lnTs                             &
               + B_OMEGA(2) * lnTs * lnTs                      &
               + B_OMEGA(3) * lnTs**3                          &
               + B_OMEGA(4) * lnTs**4 )

    ! d(ln Omega)/d(ln T*) = sum_{i>=1} i * b_i * (ln T*)^(i-1)
    dlnOmega_dlnTs =       B_OMEGA(1)                          &
                    + 2.0d0 * B_OMEGA(2) * lnTs                &
                    + 3.0d0 * B_OMEGA(3) * lnTs * lnTs         &
                    + 4.0d0 * B_OMEGA(4) * lnTs**3

    ! eta^0 in micro-Pa.s (Chapman-Enskog)
    M_T_sqrt = sqrt(AIR_MW_GMOL * T)
    eta0_uPa = ETA0_PREFACTOR * M_T_sqrt / (AIR_SIGMA_LJ**2 * omega)

    ! d eta^0/dT analytically:
    !   eta^0 = K * sqrt(M*T) / (sigma^2 * Omega)
    !   ln(eta^0) = ln(K) + 0.5*ln(M*T) - 2*ln(sigma) - ln(Omega)
    !   d(ln eta^0)/dT = 0.5/T - d(ln Omega)/dT
    !   d(ln Omega)/dT = d(ln Omega)/d(ln T*) * (1/T)
    !   so d(ln eta^0)/dT = (0.5 - dlnOmega_dlnTs) / T
    deta0_dT = eta0_uPa * (0.5d0 - dlnOmega_dlnTs) / T

    ! ----- Residual viscosity eta^r(tau, delta) -----
    tau   = AIR_TC_VISC / T
    delta = rho / AIR_RHOC_VISC

    eta_r_uPa   = 0.0d0
    deta_r_dT   = 0.0d0
    deta_r_drho = 0.0d0

    if (delta > 0.0d0) then
      do i = 1, N_VISC_TERMS
        ! exp(-gamma * delta^l). For terms with gamma=0 or l=0, exp(0)=1.
        if (G_VISC(i) > 0.0d0 .and. L_VISC(i) > 0.0d0) then
          exp_arg = -G_VISC(i) * delta**L_VISC(i)
        else
          exp_arg = 0.0d0
        end if
        term = N_VISC(i) * tau**T_VISC(i) * delta**D_VISC(i) * exp(exp_arg)
        eta_r_uPa = eta_r_uPa + term

        ! d(term)/dT at constant rho:
        !   d/dT [tau^t] = t * tau^t * d(ln tau)/dT = t * tau^t * (-1/T)
        !   (delta and exp(-gamma*delta^l) don't depend on T)
        dterm_dT_factor = -T_VISC(i) / T

        ! d(term)/drho at constant T:
        !   d/drho [delta^d] = d * delta^(d-1) * (1/rho_c) = d * delta^d / rho
        !   d/drho [exp(-gamma*delta^l)] = exp(...) * (-gamma * l * delta^(l-1) / rho_c)
        !                                = exp(...) * (-gamma * l * delta^l / rho)
        !   product rule: term * (d/rho - gamma * l * delta^l / rho)
        if (G_VISC(i) > 0.0d0 .and. L_VISC(i) > 0.0d0) then
          dterm_drho_factor = ( D_VISC(i)                                   &
                              - G_VISC(i) * L_VISC(i) * delta**L_VISC(i) ) &
                              / rho
        else
          dterm_drho_factor = D_VISC(i) / rho
        end if

        deta_r_dT   = deta_r_dT   + term * dterm_dT_factor
        deta_r_drho = deta_r_drho + term * dterm_drho_factor
      end do
    end if

    ! Total viscosity (convert from micro-Pa.s to Pa.s)
    mu       = (eta0_uPa + eta_r_uPa) * 1.0d-6
    dmu_dT   = (deta0_dT + deta_r_dT) * 1.0d-6
    dmu_drho = deta_r_drho           * 1.0d-6

  end subroutine AirViscosity

! =====================================================================
! Subroutine: AirIdealGasHeatCapacity
!
! Returns cp_ig(T) [J/(mol*K)] and its T derivative via Shomate.
! =====================================================================
  subroutine AirIdealGasHeatCapacity(T, cp, dcp_dT, ierr)
    PetscReal, intent(in)       :: T
    PetscReal, intent(out)      :: cp, dcp_dT
    PetscErrorCode, intent(out) :: ierr
    PetscReal :: t_red

    ierr = 0
    cp   = 0.0d0; dcp_dT = 0.0d0
    if (T <= 0.0d0) then
      ierr = 1
      return
    end if

    t_red = T / 1000.0d0

    ! cp = A + B*t + C*t^2 + D*t^3 + E/t^2
    cp = SHOMATE_A                                                          &
       + SHOMATE_B * t_red                                                  &
       + SHOMATE_C * t_red * t_red                                          &
       + SHOMATE_D * t_red**3                                               &
       + SHOMATE_E / (t_red * t_red)

    ! dcp/dT = (1/1000) * [B + 2*C*t + 3*D*t^2 - 2*E/t^3]
    dcp_dT = ( SHOMATE_B                                                    &
             + 2.0d0 * SHOMATE_C * t_red                                    &
             + 3.0d0 * SHOMATE_D * t_red * t_red                            &
             - 2.0d0 * SHOMATE_E / (t_red**3) ) / 1000.0d0

  end subroutine AirIdealGasHeatCapacity

! =====================================================================
! Subroutine: AirIdealGasEnthalpy
!
! Returns h_ig(T) - h_ig(T_REF) [J/mol] and its T derivative.
! By Shomate construction h_ig(T_REF) = 0; the user can offset.
! Note: T derivative of h equals cp by definition.
! =====================================================================
  subroutine AirIdealGasEnthalpy(T, h_minus_href, dh_dT, ierr)
    PetscReal, intent(in)       :: T
    PetscReal, intent(out)      :: h_minus_href, dh_dT
    PetscErrorCode, intent(out) :: ierr
    PetscReal :: t_red, h_kJmol, cp_dum

    ierr = 0
    h_minus_href = 0.0d0; dh_dT = 0.0d0

    if (T <= 0.0d0) then
      ierr = 1
      return
    end if

    t_red = T / 1000.0d0

    ! H(T) - H(T_REF) in kJ/mol from Shomate
    h_kJmol = SHOMATE_A * t_red                                             &
            + SHOMATE_B * t_red * t_red / 2.0d0                             &
            + SHOMATE_C * t_red**3 / 3.0d0                                  &
            + SHOMATE_D * t_red**4 / 4.0d0                                  &
            - SHOMATE_E / t_red                                             &
            + SHOMATE_F                                                     &
            - SHOMATE_H_OFF

    h_minus_href = h_kJmol * 1.0d3  ! convert kJ/mol -> J/mol

    ! dh/dT = cp by definition
    call AirIdealGasHeatCapacity(T, dh_dT, cp_dum, ierr)

  end subroutine AirIdealGasEnthalpy

! =====================================================================
! Diagnostic: viscosity analytic vs FD derivatives
! =====================================================================
  subroutine AirViscosityVerifyDerivatives(T, rho, max_rel_err, ierr)
    PetscReal, intent(in)       :: T, rho
    PetscReal, intent(out)      :: max_rel_err
    PetscErrorCode, intent(out) :: ierr
    PetscReal :: mu0, dT_a, drho_a
    PetscReal :: mu_p, mu_m, dT_dum, drho_dum
    PetscReal :: dT_fd, drho_fd, dT, drho, err
    PetscErrorCode :: ie

    ierr = 0; max_rel_err = 0.0d0

    call AirViscosity(T, rho, mu0, dT_a, drho_a, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    dT   = max(1.0d-4, 1.0d-6 * T)
    drho = max(1.0d-6, 1.0d-6 * max(rho, 1.0d0))

    call AirViscosity(T + dT, rho, mu_p, dT_dum, drho_dum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    call AirViscosity(T - dT, rho, mu_m, dT_dum, drho_dum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    dT_fd = (mu_p - mu_m) / (2.0d0 * dT)

    if (rho >= drho) then
      call AirViscosity(T, rho + drho, mu_p, dT_dum, drho_dum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      call AirViscosity(T, rho - drho, mu_m, dT_dum, drho_dum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      drho_fd = (mu_p - mu_m) / (2.0d0 * drho)
    else
      call AirViscosity(T, rho + drho, mu_p, dT_dum, drho_dum, ie)
      if (ie /= 0) then; ierr = ie; return; end if
      drho_fd = (mu_p - mu0) / drho
    end if

    err = abs(dT_a - dT_fd) / max(abs(dT_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(drho_a - drho_fd) / max(abs(drho_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
  end subroutine AirViscosityVerifyDerivatives

! =====================================================================
! Diagnostic: enthalpy analytic vs FD derivatives (should match cp)
! =====================================================================
  subroutine AirEnthalpyVerifyDerivatives(T, max_rel_err, ierr)
    PetscReal, intent(in)       :: T
    PetscReal, intent(out)      :: max_rel_err
    PetscErrorCode, intent(out) :: ierr
    PetscReal :: h0, dh_a
    PetscReal :: cp_a, dcp_dum
    PetscReal :: h_p, h_m, dh_dum
    PetscReal :: dh_fd, dT, err
    PetscErrorCode :: ie

    ierr = 0; max_rel_err = 0.0d0

    call AirIdealGasEnthalpy(T, h0, dh_a, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    ! Verify dh/dT = cp by computing cp directly
    call AirIdealGasHeatCapacity(T, cp_a, dcp_dum, ie)
    if (ie /= 0) then; ierr = ie; return; end if

    dT = max(1.0d-4, 1.0d-6 * T)
    call AirIdealGasEnthalpy(T + dT, h_p, dh_dum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    call AirIdealGasEnthalpy(T - dT, h_m, dh_dum, ie)
    if (ie /= 0) then; ierr = ie; return; end if
    dh_fd = (h_p - h_m) / (2.0d0 * dT)

    err = abs(dh_a - dh_fd) / max(abs(dh_fd), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
    err = abs(cp_a - dh_a) / max(abs(dh_a), tiny(1.0d0))
    max_rel_err = max(max_rel_err, err)
  end subroutine AirEnthalpyVerifyDerivatives

end module Air_Properties_module
