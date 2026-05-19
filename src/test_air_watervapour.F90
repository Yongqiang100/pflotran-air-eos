! =====================================================================
! test_air_watervapour
!
! Verification driver for Air_WaterVapour_module.
!
! Tests:
!   1. IAPWS-IF97 saturation pressure anchor checks against published
!      reference values (Wagner et al. 2000) at standard temperatures.
!   2. Water activity for NaCl brine against Pitzer reference values.
!   3. Poynting correction at high pressure for sanity.
!   4. Analytic vs finite-difference derivatives.
!   5. Halite-precipitation potential at a CAES injection state.
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

program test_air_watervapour

  use Air_EOS_PR_module
  use Air_WaterVapour_module

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! Tolerances
  real(kind=dp), parameter :: TOL_PSAT_ABS  = 0.005_dp   ! 0.5 % vs IAPWS reference
  real(kind=dp), parameter :: TOL_AW        = 0.020_dp   ! 2 % vs Pitzer
  real(kind=dp), parameter :: TOL_DERIV     = 1.0d-5

  integer :: ierr, n_pass, n_fail
  real(kind=dp) :: P_sat, dPsat_dT, err
  real(kind=dp) :: a_w, da_dm
  real(kind=dp) :: poynting, dpy_dP, dpy_dT
  real(kind=dp) :: y_w, dyw_dP, dyw_dT, dyw_dm
  real(kind=dp) :: max_err
  real(kind=dp) :: rho_g, Z_g, phi_g, h_dep_g
  real(kind=dp) :: drhoP, drhoT, dphiP, dphiT

  write(*, '(a)') '======================================================='
  write(*, '(a)') ' Air water-vapour partitioning module verification'
  write(*, '(a)') '======================================================='
  write(*, *)

  ! ==================================================================
  ! TEST 1: IAPWS-IF97 saturation pressure anchors
  ! Reference values from Wagner et al. 2000 / NIST WebBook
  ! ==================================================================
  write(*, '(a)') '--- Test 1: IAPWS-IF97 saturation pressure (tol 0.5 %) ---'
  write(*, '(2x,a)') '   T [K]   P_sat [Pa]   reference     err [%]  status'

  n_pass = 0
  n_fail = 0
  call CheckPsat(273.16_dp,    611.657_dp)        ! triple point
  call CheckPsat(298.15_dp,    3169.93_dp)        ! 25 C
  call CheckPsat(300.0_dp,     3536.81_dp)        ! IAPWS reference point
  call CheckPsat(323.15_dp,    12351.6_dp)        ! 50 C
  call CheckPsat(373.124_dp,   101325.0_dp)       ! normal boiling point
  call CheckPsat(423.15_dp,    476137.0_dp)       ! 150 C
  call CheckPsat(500.0_dp,     2638970.0_dp)      ! IAPWS reference point
  call CheckPsat(623.15_dp,    16529165.0_dp)     ! near-critical

  write(*, '(2x,a,i0,a,i0)') 'P_sat anchors: ', n_pass, '/', n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 2: Water activity for NaCl brine
  ! Pitzer reference values from Robinson & Stokes (1959) / Pitzer (1991)
  ! Note: our ideal-mixing model is an approximation; tolerance reflects
  ! its known accuracy versus the rigorous Pitzer reference.
  ! ==================================================================
  write(*, '(a)') '--- Test 2: water activity in NaCl brine (vs Pitzer, tol 2 %) ---'
  write(*, '(2x,a)') '   m [mol/kg]  a_w(model)  a_w(Pitzer)   err [%]  status'

  n_pass = 0
  n_fail = 0
  call CheckActivity(0.5_dp,   0.9836_dp)
  call CheckActivity(1.0_dp,   0.9670_dp)
  call CheckActivity(2.0_dp,   0.9316_dp)
  call CheckActivity(3.0_dp,   0.8932_dp)

  write(*, '(2x,a,i0,a,i0)') 'Activity anchors: ', n_pass, '/', n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 3: Poynting at high pressure (sanity)
  ! ==================================================================
  write(*, '(a)') '--- Test 3: Poynting correction at CAES conditions ---'
  call WaterSatPressure(313.15_dp, P_sat, dPsat_dT, ierr)
  call WaterPoyntingFactor(8.0d6, 313.15_dp, P_sat, 1.807d-5, &
                            poynting, dpy_dP, dpy_dT, ierr)
  write(*, '(2x,a,f6.4,a)') 'Poynting(80 bar, 40 C): ', poynting,        &
    '   (~ 1.06 expected)'
  call WaterPoyntingFactor(1.0d7, 313.15_dp, P_sat, 1.807d-5, &
                            poynting, dpy_dP, dpy_dT, ierr)
  write(*, '(2x,a,f6.4,a)') 'Poynting(100 bar, 40 C): ', poynting,       &
    '   (~ 1.08 expected)'
  write(*, *)

  ! ==================================================================
  ! TEST 4: analytic vs FD derivatives of y_w
  ! ==================================================================
  write(*, '(a)') '--- Test 4: analytic vs FD derivatives of y_w (tol 1e-5) ---'
  write(*, '(2x,a)') '    P [Pa]      T [K]    m [mol/kg]    max_err    status'

  n_pass = 0
  n_fail = 0
  call DerivCheck(1.0d5,  298.15_dp, 0.0_dp)
  call DerivCheck(1.0d5,  298.15_dp, 1.0_dp)
  call DerivCheck(8.0d6,  313.15_dp, 0.5_dp)
  call DerivCheck(1.0d7,  333.15_dp, 1.0_dp)
  call DerivCheck(5.0d7,  353.15_dp, 2.0_dp)

  write(*, '(2x,a,i0,a,i0)') 'Derivative test: ', n_pass, '/', n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 5: Coupled water vapour at a CAES injection state with the
  !         PR+Peneloux gas-phase fugacity coefficient.
  ! ==================================================================
  write(*, '(a)') '--- Test 5: water vapour at 80 bar, 40 C, 0.5 M NaCl ---'

  call AirEOSPRProperties(8.0d6, 313.15_dp, rho_g, Z_g, phi_g, h_dep_g,  &
                          drhoP, drhoT, dphiP, dphiT, ierr)
  write(*, '(2x,a,f7.4)') 'Bulk gas phi from PR+Peneloux: ', phi_g
  call WaterVapourMoleFraction(8.0d6, 313.15_dp, 0.5_dp, phi_g,          &
                                y_w, dyw_dP, dyw_dT, dyw_dm, ierr)
  if (ierr == 0) then
    write(*, '(2x,a,es12.5)')         'y_w [-]               : ', y_w
    write(*, '(2x,a,f6.3,a)')          'y_w as %              : ', y_w*100.0_dp, ' %'
    write(*, '(2x,a,es12.5,a)')        'Equivalent partial P  : ', y_w*8.0d6, ' Pa'
    write(*, '(2x,a,es12.5)')          'dy_w/dP [1/Pa]        : ', dyw_dP
    write(*, '(2x,a,es12.5)')          'dy_w/dT [1/K]         : ', dyw_dT
    write(*, '(2x,a,es12.5)')          'dy_w/dm [kg/mol]      : ', dyw_dm
  end if
  write(*, *)

  ! Show how water vapour scales across CAES conditions
  write(*, '(a)') '  y_w sensitivity to T at P = 80 bar, m = 0.5 mol/kg NaCl:'
  call ShowYwTrend(8.0d6,  283.15_dp, 0.5_dp)  ! 10 C
  call ShowYwTrend(8.0d6,  298.15_dp, 0.5_dp)  ! 25 C
  call ShowYwTrend(8.0d6,  313.15_dp, 0.5_dp)  ! 40 C
  call ShowYwTrend(8.0d6,  333.15_dp, 0.5_dp)  ! 60 C
  call ShowYwTrend(8.0d6,  353.15_dp, 0.5_dp)  ! 80 C
  write(*, *)

  write(*, '(a)') '======================================================='

contains

  subroutine CheckPsat(T_in, P_ref)
    real(kind=dp), intent(in) :: T_in, P_ref
    real(kind=dp) :: P_calc, dPdT
    integer       :: ie
    call WaterSatPressure(T_in, P_calc, dPdT, ie)
    err = abs(P_calc - P_ref) / P_ref
    if (ie == 0 .and. err < TOL_PSAT_ABS) then
      write(*, '(2x,f8.2,2es14.5,f9.3,a)') T_in, P_calc, P_ref,             &
                                             err*100.0_dp, '  PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,f8.2,2es14.5,f9.3,a)') T_in, P_calc, P_ref,             &
                                             err*100.0_dp, '  FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine CheckPsat

  subroutine CheckActivity(m_in, a_ref)
    real(kind=dp), intent(in) :: m_in, a_ref
    real(kind=dp) :: a_calc, da
    integer       :: ie
    call WaterActivityNaCl(m_in, a_calc, da, ie)
    err = abs(a_calc - a_ref) / a_ref
    if (ie == 0 .and. err < TOL_AW) then
      write(*, '(2x,f10.3,2f12.4,f10.3,a)') m_in, a_calc, a_ref,            &
                                              err*100.0_dp, '  PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,f10.3,2f12.4,f10.3,a)') m_in, a_calc, a_ref,            &
                                              err*100.0_dp, '  FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine CheckActivity

  subroutine DerivCheck(P_in, T_in, m_in)
    real(kind=dp), intent(in) :: P_in, T_in, m_in
    real(kind=dp) :: mre
    integer       :: ie
    call WaterVapourVerifyDerivatives(P_in, T_in, m_in, 1.0_dp, mre, ie)
    if (ie == 0 .and. mre < TOL_DERIV) then
      write(*, '(2x,2es12.4,f12.3,es12.3,a)') P_in, T_in, m_in, mre,        &
                                                ' PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,2es12.4,f12.3,es12.3,a)') P_in, T_in, m_in, mre,        &
                                                ' FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine DerivCheck

  subroutine ShowYwTrend(P_in, T_in, m_in)
    real(kind=dp), intent(in) :: P_in, T_in, m_in
    real(kind=dp) :: y_local, dPdum, dTdum, dmdum
    integer       :: ie
    call WaterVapourMoleFraction(P_in, T_in, m_in, 1.0_dp,                  &
                                  y_local, dPdum, dTdum, dmdum, ie)
    write(*, '(2x,a,f6.2,a,es10.3,a,f5.3,a)') 'T = ', T_in - 273.15_dp,     &
      ' C: y_w = ', y_local, ' (', y_local*100.0_dp, ' %)'
  end subroutine ShowYwTrend

end program test_air_watervapour
