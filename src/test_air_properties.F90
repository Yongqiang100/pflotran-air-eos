! =====================================================================
! test_air_properties
!
! Verification driver for Air_Properties_module.
!
! Tests:
!   1. Viscosity anchor checks against NIST WebBook reference values
!      at atmospheric and high-pressure CAES conditions.
!   2. cp_ig at 298.15 K vs textbook value (29.12 J/(mol*K)).
!   3. Ideal-gas enthalpy increments h(T) - h(298.15) vs NIST.
!   4. Analytic vs FD derivatives (viscosity and enthalpy).
!   5. Coupled state: use PR+Peneloux rho(T, P) then feed into
!      viscosity, illustrating end-to-end usage at a CAES injection
!      state.
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

program test_air_properties

  use Air_EOS_PR_module
  use Air_Properties_module

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! Tolerances
  real(kind=dp), parameter :: TOL_VISC   = 0.030_dp   ! 3 % (L-J 2004 stated ~ 5 % near critical;
                                                      !  our regime is well above, ~ 1-2 % typical)
  real(kind=dp), parameter :: TOL_CP_REF = 0.010_dp   ! 1 % at reference state
  real(kind=dp), parameter :: TOL_H      = 0.020_dp   ! 2 % on enthalpy increments
  real(kind=dp), parameter :: TOL_H_ABS  = 1.0_dp     ! 1 J/mol absolute floor (handles 0-reference)
  real(kind=dp), parameter :: TOL_DERIV  = 1.0d-5

  integer :: ierr, n_pass, n_fail
  real(kind=dp) :: mu, dmu_dT, dmu_drho, err
  real(kind=dp) :: cp, dcp_dT, h, dh_dT
  real(kind=dp) :: max_err
  real(kind=dp) :: rho_g, Z_g, phi_g, h_dep_g
  real(kind=dp) :: drhoP, drhoT, dphiP, dphiT

  write(*, '(a)') '======================================================='
  write(*, '(a)') ' Air properties module verification'
  write(*, '(a)') ' (viscosity + ideal-gas enthalpy)'
  write(*, '(a)') '======================================================='
  write(*, *)

  ! ==================================================================
  ! TEST 1: viscosity anchors vs L-J 2004 reference
  ! Reference values from Lemmon & Jacobsen (2004) Tables for air,
  ! which is the published NIST standard. At atmospheric pressure
  ! the values match Kestin et al. measurements; at elevated P they
  ! reflect the L-J formulation our code implements (self-consistent
  ! check, not an external accuracy check at high P).
  ! ==================================================================
  write(*, '(a)') '--- Test 1: viscosity vs L-J 2004 reference (tol 3 %) ---'
  write(*, '(2x,a)')                                                       &
    '   T [K]  rho [kg/m^3]   mu [Pa.s]    L-J ref       err [%]  status'

  n_pass = 0
  n_fail = 0
  ! Atmospheric: dilute-gas regime, well-established values
  call CheckMu(298.15_dp,   1.184_dp,   1.846d-5)
  call CheckMu(273.15_dp,   1.293_dp,   1.722d-5)
  call CheckMu(350.0_dp,    1.008_dp,   2.076d-5)
  ! High pressure: L-J 2004 tabulated values for air
  call CheckMu(298.15_dp,   58.78_dp,   1.910d-5)
  call CheckMu(298.15_dp,  115.4_dp,    2.040d-5)
  call CheckMu(323.15_dp,  105.2_dp,    2.100d-5)

  write(*, '(2x,a,i0,a,i0)') 'Viscosity anchors: ', n_pass, '/', n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 2: cp_ig at reference temperature
  ! ==================================================================
  write(*, '(a)') '--- Test 2: cp_ig at 298.15 K vs textbook ---'
  call AirIdealGasHeatCapacity(298.15_dp, cp, dcp_dT, ierr)
  err = abs(cp - 29.12_dp) / 29.12_dp
  if (ierr == 0 .and. err < TOL_CP_REF) then
    write(*, '(2x,a,f7.3,a,f6.3,a)') 'cp_ig(298.15 K) = ', cp,             &
      ' J/(mol K)   err = ', err*100.0_dp, ' %   PASS'
  else
    write(*, '(2x,a,f7.3,a,f6.3,a)') 'cp_ig(298.15 K) = ', cp,             &
      ' J/(mol K)   err = ', err*100.0_dp, ' %   FAIL'
  end if
  write(*, *)

  ! ==================================================================
  ! TEST 3: ideal-gas enthalpy increments
  ! ==================================================================
  write(*, '(a)') '--- Test 3: h_ig(T) - h_ig(298.15 K) vs NIST (tol 2 %) ---'
  write(*, '(2x,a)') '   T [K]      h [J/mol]    NIST [J/mol]   err [%]  status'

  n_pass = 0
  n_fail = 0
  ! NIST Shomate-tabulated enthalpy increments for N2 / air (J/mol)
  call CheckH(298.15_dp,     0.0_dp)
  call CheckH(323.15_dp,   728.0_dp)
  call CheckH(373.15_dp,  2189.0_dp)
  call CheckH(473.15_dp,  5121.0_dp)
  call CheckH(573.15_dp,  8131.0_dp)

  write(*, '(2x,a,i0,a,i0)') 'Enthalpy anchors: ', n_pass, '/', n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 4: derivative checks
  ! ==================================================================
  write(*, '(a)') '--- Test 4: analytic vs FD derivatives (tol 1e-5) ---'
  write(*, '(2x,a)') '   property   T [K]   rho [kg/m^3]   max_err     status'

  n_pass = 0
  n_fail = 0
  call DerivCheckVisc(298.15_dp,   1.184_dp)
  call DerivCheckVisc(298.15_dp, 115.4_dp)
  call DerivCheckVisc(323.15_dp, 105.2_dp)
  call DerivCheckVisc(273.15_dp, 127.6_dp)
  call DerivCheckEnth(298.15_dp)
  call DerivCheckEnth(373.15_dp)

  write(*, '(2x,a,i0,a,i0)') 'Derivative test: ', n_pass, '/', n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 5: end-to-end coupled state at CAES conditions
  ! ==================================================================
  write(*, '(a)') '--- Test 5: coupled state at 80 bar, 40 C ---'
  call AirEOSPRProperties(8.0d6, 313.15_dp, rho_g, Z_g, phi_g, h_dep_g,  &
                          drhoP, drhoT, dphiP, dphiT, ierr)
  if (ierr == 0) then
    write(*, '(2x,a,f8.3,a)') 'From PR+Peneloux: rho = ', rho_g, ' kg/m^3'
    write(*, '(2x,a,es12.5,a)') 'h_dep (residual)        = ', h_dep_g,   &
      ' J/mol'

    call AirViscosity(313.15_dp, rho_g, mu, dmu_dT, dmu_drho, ierr)
    write(*, '(2x,a,es12.5,a)') 'Viscosity               = ', mu,         &
      ' Pa.s'
    write(*, '(2x,a,es12.5,a)') 'dmu/dT                  = ', dmu_dT,    &
      ' Pa.s/K'
    write(*, '(2x,a,es12.5,a)') 'dmu/drho                = ', dmu_drho,  &
      ' Pa.s.m^3/kg'

    call AirIdealGasEnthalpy(313.15_dp, h, dh_dT, ierr)
    write(*, '(2x,a,f8.2,a,f7.2,a)') 'h_ig(40 C) - h_ig(25 C) = ', h,    &
      ' J/mol   (cp = ', dh_dT, ' J/(mol K))'
    write(*, '(2x,a,f8.2,a)') 'Total h - h_ref         = ',                &
      h + h_dep_g, ' J/mol'
  end if
  write(*, *)

  write(*, '(a)') '======================================================='

contains

  subroutine CheckMu(T_in, rho_in, mu_ref)
    real(kind=dp), intent(in) :: T_in, rho_in, mu_ref
    real(kind=dp) :: mu_calc, dT_dum, drho_dum
    integer       :: ie
    call AirViscosity(T_in, rho_in, mu_calc, dT_dum, drho_dum, ie)
    err = abs(mu_calc - mu_ref) / mu_ref
    if (ie == 0 .and. err < TOL_VISC) then
      write(*, '(2x,f8.2,f10.3,2es14.5,f9.3,a)')                          &
        T_in, rho_in, mu_calc, mu_ref, err*100.0_dp, '  PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,f8.2,f10.3,2es14.5,f9.3,a)')                          &
        T_in, rho_in, mu_calc, mu_ref, err*100.0_dp, '  FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine CheckMu

  subroutine CheckH(T_in, h_ref)
    real(kind=dp), intent(in) :: T_in, h_ref
    real(kind=dp) :: h_calc, dh_dum, abs_err
    integer       :: ie
    call AirIdealGasEnthalpy(T_in, h_calc, dh_dum, ie)
    abs_err = abs(h_calc - h_ref)
    ! For near-zero reference (T = T_ref), use absolute tolerance;
    ! otherwise use relative tolerance.
    if (abs(h_ref) < 10.0_dp) then
      err = abs_err  ! show as J/mol when reference is near zero
      if (ie == 0 .and. abs_err < TOL_H_ABS) then
        write(*, '(2x,f8.2,2f15.3,f10.3,a)')                              &
          T_in, h_calc, h_ref, abs_err, ' J/mol PASS'
        n_pass = n_pass + 1
      else
        write(*, '(2x,f8.2,2f15.3,f10.3,a)')                              &
          T_in, h_calc, h_ref, abs_err, ' J/mol FAIL'
        n_fail = n_fail + 1
      end if
    else
      err = abs_err / abs(h_ref)
      if (ie == 0 .and. err < TOL_H) then
        write(*, '(2x,f8.2,2f15.3,f10.3,a)')                              &
          T_in, h_calc, h_ref, err*100.0_dp, '   %  PASS'
        n_pass = n_pass + 1
      else
        write(*, '(2x,f8.2,2f15.3,f10.3,a)')                              &
          T_in, h_calc, h_ref, err*100.0_dp, '   %  FAIL'
        n_fail = n_fail + 1
      end if
    end if
  end subroutine CheckH

  subroutine DerivCheckVisc(T_in, rho_in)
    real(kind=dp), intent(in) :: T_in, rho_in
    real(kind=dp) :: mre
    integer       :: ie
    call AirViscosityVerifyDerivatives(T_in, rho_in, mre, ie)
    if (ie == 0 .and. mre < TOL_DERIV) then
      write(*, '(2x,a,f8.2,f12.3,es12.3,a)') 'visc      ',                 &
        T_in, rho_in, mre, '   PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,a,f8.2,f12.3,es12.3,a)') 'visc      ',                 &
        T_in, rho_in, mre, '   FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine DerivCheckVisc

  subroutine DerivCheckEnth(T_in)
    real(kind=dp), intent(in) :: T_in
    real(kind=dp) :: mre
    integer       :: ie
    call AirEnthalpyVerifyDerivatives(T_in, mre, ie)
    if (ie == 0 .and. mre < TOL_DERIV) then
      write(*, '(2x,a,f8.2,a,es12.3,a)') 'enthalpy  ',                     &
        T_in, '       --  ', mre, '   PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,a,f8.2,a,es12.3,a)') 'enthalpy  ',                     &
        T_in, '       --  ', mre, '   FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine DerivCheckEnth

end program test_air_properties
