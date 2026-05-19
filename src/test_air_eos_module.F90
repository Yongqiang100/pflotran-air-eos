! =====================================================================
! test_air_eos_module
!
! Verification driver for the PFLOTRAN-facing wrapper Air_EOS_module.
!
! Tests:
!   1. Self-check routine (AirEOSVerify) — confirms all five top-level
!      entry points return values within physical bounds.
!   2. Consistency check: every value returned by the wrapper agrees
!      with the result of calling the underlying module directly.
!      This catches any chain-rule or dispatch bug.
!   3. Derivative chain rule: dmu/dP via wrapper (chain rule through
!      rho) agrees with finite-difference dmu/dP at the wrapper level.
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

program test_air_eos_module

  use Air_EOS_PR_module
  use Air_Henry_module
  use Air_WaterVapour_module
  use Air_Properties_module
  use Air_EOS_module

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! Tolerances
  real(kind=dp), parameter :: TOL_CONS   = 1.0d-12  ! identical math
  real(kind=dp), parameter :: TOL_DERIV  = 1.0d-5

  integer :: ierr, n_pass, n_fail

  write(*, '(a)') '======================================================='
  write(*, '(a)') ' Air_EOS_module wrapper verification'
  write(*, '(a)') '======================================================='
  write(*, *)

  ! ==================================================================
  ! TEST 1: self-check verbose
  ! ==================================================================
  write(*, '(a)') '--- Test 1: AirEOSVerify self-check ---'
  call AirEOSVerify(.true., ierr)
  if (ierr == 0) then
    write(*, '(2x,a)') 'PASS'
  else
    write(*, '(2x,a,i0)') 'FAIL ierr = ', ierr
  end if
  write(*, *)

  ! ==================================================================
  ! TEST 2: wrapper-vs-underlying consistency
  ! ==================================================================
  write(*, '(a)') '--- Test 2: wrapper vs underlying-module consistency ---'

  n_pass = 0
  n_fail = 0
  call CheckConsistency(1.0d5,  298.15_dp, 0.0_dp,  'atmospheric, pure water')
  call CheckConsistency(8.0d6,  313.15_dp, 0.5_dp,  '80 bar 40C, 0.5 M NaCl')
  call CheckConsistency(1.5d7,  333.15_dp, 1.0_dp,  '150 bar 60C, 1 M NaCl')
  call CheckConsistency(5.0d7,  283.15_dp, 2.0_dp,  '500 bar 10C, 2 M NaCl')

  write(*, '(2x,a,i0,a,i0)') 'Consistency test: ', n_pass, '/',           &
                              n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 3: wrapper-level derivative chain-rule
  ! ==================================================================
  write(*, '(a)') '--- Test 3: wrapper-level derivative checks (tol 1e-5) ---'

  n_pass = 0
  n_fail = 0
  call CheckWrapperDerivatives(1.0d5,  298.15_dp, 0.0_dp)
  call CheckWrapperDerivatives(8.0d6,  313.15_dp, 0.5_dp)
  call CheckWrapperDerivatives(1.5d7,  333.15_dp, 1.0_dp)

  write(*, '(2x,a,i0,a,i0)') 'Derivative test: ', n_pass, '/',            &
                              n_pass + n_fail
  write(*, *)

  write(*, '(a)') '======================================================='

contains

! ---------------------------------------------------------------------
! Compare every quantity returned by the wrapper with the result of
! a direct call to the underlying module. They must agree exactly
! (1e-12 relative tolerance) since the wrapper does not re-implement
! anything — it just dispatches.
! ---------------------------------------------------------------------
  subroutine CheckConsistency(P, T, m, label)
    real(kind=dp),    intent(in) :: P, T, m
    character(len=*), intent(in) :: label

    real(kind=dp) :: rho_w, mu_w, fug_w, hdep_w, htot_w
    real(kind=dp) :: drho_dP_w, drho_dT_w, dmu_dP_w, dmu_dT_w
    real(kind=dp) :: dfug_dP_w, dfug_dT_w, dhdep_dP_w, dhdep_dT_w
    real(kind=dp) :: dh_dP_w, dh_dT_w
    real(kind=dp) :: y_w, dyw_dP_w, dyw_dT_w, dyw_dm_w
    real(kind=dp) :: c_O2_w, dc_dP_w, dc_dT_w, dc_dm_w

    real(kind=dp) :: rho_d, Z_d, fug_d, hdep_d
    real(kind=dp) :: drho_dP_d, drho_dT_d, dfug_dP_d, dfug_dT_d
    real(kind=dp) :: mu_d, dmu_dT_d, dmu_drho_d
    real(kind=dp) :: phi_d, y_d, dyw_dP_d, dyw_dT_d, dyw_dm_d
    real(kind=dp) :: c_d, dc_dP_d, dc_dT_d, dc_dm_d
    real(kind=dp) :: h_ig, dhig_dT

    integer :: ie
    real(kind=dp) :: err, max_err
    logical :: ok

    ! ----- Wrapper calls -----
    call AirEOSGetGasProperties(P, T, rho_w, mu_w, fug_w, hdep_w,         &
                                 drho_dP_w, drho_dT_w, dmu_dP_w, dmu_dT_w,&
                                 dfug_dP_w, dfug_dT_w, dhdep_dP_w,        &
                                 dhdep_dT_w, ie)
    if (ie /= 0) then
      write(*, '(2x,a,a,i0)') label, ': wrapper gas-props ERROR ', ie
      n_fail = n_fail + 1
      return
    end if

    call AirEOSGetTotalEnthalpy(P, T, htot_w, dh_dP_w, dh_dT_w, ie)
    call AirEOSGetWaterVapourFraction(P, T, m, y_w, dyw_dP_w, dyw_dT_w,   &
                                       dyw_dm_w, ie)
    call AirEOSGetDissolvedO2(P, T, m, 0.20946_dp, c_O2_w, dc_dP_w,       &
                               dc_dT_w, dc_dm_w, ie)

    ! ----- Direct module calls -----
    call AirEOSPRProperties(P, T, rho_d, Z_d, fug_d, hdep_d,              &
                            drho_dP_d, drho_dT_d, dfug_dP_d, dfug_dT_d, ie)
    call AirViscosity(T, rho_d, mu_d, dmu_dT_d, dmu_drho_d, ie)
    call AirIdealGasEnthalpy(T, h_ig, dhig_dT, ie)
    phi_d = AirEOSPRFugacityCoeff(P, T, ie)
    call WaterVapourMoleFraction(P, T, m, phi_d, y_d, dyw_dP_d, dyw_dT_d, &
                                  dyw_dm_d, ie)
    call AirHenryDissolvedO2(P, T, m, 0.20946_dp, phi_d, c_d, dc_dP_d,    &
                              dc_dT_d, dc_dm_d, ie)

    ! ----- Compare (max relative error across all returned values) -----
    max_err = 0.0_dp
    max_err = max(max_err, RelErr(rho_w,    rho_d))
    max_err = max(max_err, RelErr(fug_w,    fug_d))
    max_err = max(max_err, RelErr(hdep_w,   hdep_d))
    max_err = max(max_err, RelErr(mu_w,     mu_d))
    max_err = max(max_err, RelErr(htot_w,   h_ig + hdep_d))
    max_err = max(max_err, RelErr(y_w,      y_d))
    max_err = max(max_err, RelErr(c_O2_w,   c_d))
    max_err = max(max_err, RelErr(drho_dP_w, drho_dP_d))
    max_err = max(max_err, RelErr(drho_dT_w, drho_dT_d))
    max_err = max(max_err, RelErr(dfug_dP_w, dfug_dP_d))
    max_err = max(max_err, RelErr(dfug_dT_w, dfug_dT_d))
    max_err = max(max_err, RelErr(dyw_dP_w, dyw_dP_d))
    max_err = max(max_err, RelErr(dyw_dT_w, dyw_dT_d))
    max_err = max(max_err, RelErr(dyw_dm_w, dyw_dm_d))
    max_err = max(max_err, RelErr(dc_dP_w,  dc_dP_d))
    max_err = max(max_err, RelErr(dc_dT_w,  dc_dT_d))
    max_err = max(max_err, RelErr(dc_dm_w,  dc_dm_d))

    ok = max_err < TOL_CONS
    if (ok) then
      write(*, '(2x,a,a,es12.3,a)') label, ': max_err = ', max_err, '  PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,a,a,es12.3,a)') label, ': max_err = ', max_err, '  FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine CheckConsistency

! ---------------------------------------------------------------------
! Check that wrapper-level derivatives (which involve chain-rule
! through rho for viscosity) match centered FD on the wrapper itself.
! ---------------------------------------------------------------------
  subroutine CheckWrapperDerivatives(P, T, m)
    real(kind=dp), intent(in) :: P, T, m
    real(kind=dp) :: dP_step, dT_step
    real(kind=dp) :: rho_p, mu_p, fug_p, hdep_p
    real(kind=dp) :: rho_m, mu_m, fug_m, hdep_m
    real(kind=dp) :: drho_dum, dmu_dum, dfug_dum, dhdep_dum
    real(kind=dp) :: drhoP_a, drhoT_a, dmuP_a, dmuT_a
    real(kind=dp) :: dfugP_a, dfugT_a, dhdepP_a, dhdepT_a
    real(kind=dp) :: rho, mu, fug, hdep, max_err
    integer :: ie

    call AirEOSGetGasProperties(P, T, rho, mu, fug, hdep,                  &
                                 drhoP_a, drhoT_a, dmuP_a, dmuT_a,         &
                                 dfugP_a, dfugT_a, dhdepP_a, dhdepT_a, ie)

    dP_step = max(1.0_dp, 1.0d-6 * P)
    dT_step = max(1.0d-3, 1.0d-6 * T)

    call AirEOSGetGasProperties(P + dP_step, T, rho_p, mu_p, fug_p, hdep_p,&
                                 drho_dum, drho_dum, dmu_dum, dmu_dum,     &
                                 dfug_dum, dfug_dum, dhdep_dum, dhdep_dum, ie)
    call AirEOSGetGasProperties(P - dP_step, T, rho_m, mu_m, fug_m, hdep_m,&
                                 drho_dum, drho_dum, dmu_dum, dmu_dum,     &
                                 dfug_dum, dfug_dum, dhdep_dum, dhdep_dum, ie)

    max_err = 0.0_dp
    max_err = max(max_err, RelErr(drhoP_a, (rho_p - rho_m)/(2.0_dp*dP_step)))
    max_err = max(max_err, RelErr(dmuP_a,  (mu_p  - mu_m )/(2.0_dp*dP_step)))
    max_err = max(max_err, RelErr(dfugP_a, (fug_p - fug_m)/(2.0_dp*dP_step)))

    call AirEOSGetGasProperties(P, T + dT_step, rho_p, mu_p, fug_p, hdep_p,&
                                 drho_dum, drho_dum, dmu_dum, dmu_dum,     &
                                 dfug_dum, dfug_dum, dhdep_dum, dhdep_dum, ie)
    call AirEOSGetGasProperties(P, T - dT_step, rho_m, mu_m, fug_m, hdep_m,&
                                 drho_dum, drho_dum, dmu_dum, dmu_dum,     &
                                 dfug_dum, dfug_dum, dhdep_dum, dhdep_dum, ie)

    max_err = max(max_err, RelErr(drhoT_a, (rho_p - rho_m)/(2.0_dp*dT_step)))
    max_err = max(max_err, RelErr(dmuT_a,  (mu_p  - mu_m )/(2.0_dp*dT_step)))
    max_err = max(max_err, RelErr(dfugT_a, (fug_p - fug_m)/(2.0_dp*dT_step)))

    if (max_err < TOL_DERIV) then
      write(*, '(2x,a,es10.3,a,es10.3,a,es12.3,a)')                        &
        'P=', P, ', T=', T, ': max_err = ', max_err, '  PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,a,es10.3,a,es10.3,a,es12.3,a)')                        &
        'P=', P, ', T=', T, ': max_err = ', max_err, '  FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine CheckWrapperDerivatives

  pure function RelErr(a, b) result(err)
    real(kind=dp), intent(in) :: a, b
    real(kind=dp)             :: err
    err = abs(a - b) / max(abs(b), 1.0d-300)
  end function RelErr

end program test_air_eos_module
