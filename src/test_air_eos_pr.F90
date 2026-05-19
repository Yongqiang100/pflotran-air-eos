! =====================================================================
! test_air_eos_pr
!
! Verification driver for Air_EOS_PR_module with Peneloux correction.
!
! Three checks:
!   1. Compare computed density against NIST air reference values
!      across the CAES P-T regime (1 atm to 100 bar, 0-50 C).
!      Expected accuracy with Peneloux: < 1.5 % relative error.
!   2. Cross-check analytic vs centered-finite-difference derivatives.
!      Expected accuracy: < 1e-5 relative.
!   3. Print full property set at a representative CAES state.
!
! Reference data: NIST Chemistry WebBook (Lemmon & Jacobsen 2000
! reference EOS for air); values manually transcribed at the
! tabulated P, T grid points used here.
! =====================================================================
program test_air_eos_pr

  use Air_EOS_PR_module

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! ------------------------------------------------------------------
  ! NIST air reference points
  !   P [Pa], T [K], rho_NIST [kg/m^3]
  ! ------------------------------------------------------------------
  integer, parameter :: N_REF = 6
  real(kind=dp), dimension(N_REF), parameter :: P_ref =                 &
    [ 1.01325d5, 5.0d6,    1.0d7,    1.0d7,    1.0d7,    2.0d5    ]
  real(kind=dp), dimension(N_REF), parameter :: T_ref =                 &
    [ 298.15d0,  298.15d0, 298.15d0, 323.15d0, 273.15d0, 273.15d0 ]
  real(kind=dp), dimension(N_REF), parameter :: rho_ref =               &
    [ 1.1839d0,  58.78d0,  115.40d0, 105.20d0, 127.60d0, 2.557d0  ]

  ! ------------------------------------------------------------------
  ! Derivative cross-check points
  ! ------------------------------------------------------------------
  integer, parameter :: N_DERIV = 5
  real(kind=dp), dimension(N_DERIV), parameter :: P_deriv =             &
    [ 1.0d5, 1.0d6, 1.0d7, 5.0d7, 1.0d7 ]
  real(kind=dp), dimension(N_DERIV), parameter :: T_deriv =             &
    [ 298.15d0, 298.15d0, 298.15d0, 298.15d0, 350.0d0 ]

  ! Tolerances
  real(kind=dp), parameter :: TOL_DENSITY = 1.5d-2     ! 1.5 % vs NIST (with Peneloux)
  real(kind=dp), parameter :: TOL_DERIV   = 1.0d-5     ! 1e-5 analytic vs FD

  ! Working variables
  real(kind=dp) :: rho, Z, phi, h_dep
  real(kind=dp) :: drho_dP, drho_dT, dphi_dP, dphi_dT
  real(kind=dp) :: max_rel_err
  real(kind=dp) :: rho_err
  integer       :: ierr, i, n_pass, n_fail

  write(*, '(a)') '======================================================='
  write(*, '(a)') ' Air PR+Peneloux EOS module verification'
  write(*, '(a)') ' Peneloux c = -9.0e-6 m^3/mol (NIST-calibrated)'
  write(*, '(a)') '======================================================='
  write(*, *)

  ! ==================================================================
  ! TEST 1: density vs NIST reference
  ! ==================================================================
  write(*, '(a)') '--- Test 1: density vs NIST reference (tol 1.5 %) ---'
  write(*, '(2x,a)')                                                       &
    '   i    P [Pa]      T [K]    rho_PR     rho_NIST   err [%]  status'

  n_pass = 0
  n_fail = 0
  do i = 1, N_REF
    call AirEOSPRProperties(P_ref(i), T_ref(i), rho, Z, phi, h_dep,     &
                            drho_dP, drho_dT, dphi_dP, dphi_dT, ierr)
    if (ierr /= 0) then
      write(*, '(2x,i4,a,i0)') i, '  ERROR ierr = ', ierr
      n_fail = n_fail + 1
      cycle
    end if

    rho_err = abs(rho - rho_ref(i)) / rho_ref(i)

    if (rho_err < TOL_DENSITY) then
      write(*, '(2x,i4,2es12.4,2f10.3,f8.3,a)') i, P_ref(i), T_ref(i),   &
        rho, rho_ref(i), rho_err*100.0_dp, '  PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,i4,2es12.4,2f10.3,f8.3,a)') i, P_ref(i), T_ref(i),   &
        rho, rho_ref(i), rho_err*100.0_dp, '  FAIL'
      n_fail = n_fail + 1
    end if
  end do

  write(*, '(2x,a,i0,a,i0)') 'Density test: ', n_pass, '/', N_REF
  write(*, *)

  ! ==================================================================
  ! TEST 2: analytic vs finite-difference derivative check
  ! ==================================================================
  write(*, '(a)') '--- Test 2: analytic vs FD derivatives (tol 1e-5) ---'
  write(*, '(2x,a)') '   i    P [Pa]      T [K]      max_rel_err   status'

  n_pass = 0
  n_fail = 0
  do i = 1, N_DERIV
    call AirEOSPRVerifyDerivatives(P_deriv(i), T_deriv(i), max_rel_err, ierr)
    if (ierr /= 0) then
      write(*, '(2x,i4,a,i0)') i, '  ERROR ierr = ', ierr
      n_fail = n_fail + 1
      cycle
    end if
    if (max_rel_err < TOL_DERIV) then
      write(*, '(2x,i4,2es12.4,es15.4,a)')                               &
        i, P_deriv(i), T_deriv(i), max_rel_err, '   PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,i4,2es12.4,es15.4,a)')                               &
        i, P_deriv(i), T_deriv(i), max_rel_err, '   FAIL'
      n_fail = n_fail + 1
    end if
  end do

  write(*, '(2x,a,i0,a,i0)') 'Derivative test: ', n_pass, '/', N_DERIV
  write(*, *)

  ! ==================================================================
  ! TEST 3: full property set at a representative CAES state
  ! ==================================================================
  write(*, '(a)') '--- Test 3: full property set at P=80 bar, T=40 C ---'
  call AirEOSPRProperties(8.0d6, 313.15d0, rho, Z, phi, h_dep,           &
                          drho_dP, drho_dT, dphi_dP, dphi_dT, ierr)
  if (ierr == 0) then
    write(*, '(2x,a,es14.6,a)') 'rho       = ', rho,     ' kg/m^3'
    write(*, '(2x,a,es14.6,a)') 'Z         = ', Z,       ' [-]'
    write(*, '(2x,a,es14.6,a)') 'phi       = ', phi,     ' [-]'
    write(*, '(2x,a,es14.6,a)') 'h_dep     = ', h_dep,   ' J/mol'
    write(*, '(2x,a,es14.6,a)') 'drho/dP   = ', drho_dP, ' kg/m^3/Pa'
    write(*, '(2x,a,es14.6,a)') 'drho/dT   = ', drho_dT, ' kg/m^3/K'
    write(*, '(2x,a,es14.6,a)') 'dphi/dP   = ', dphi_dP, ' 1/Pa'
    write(*, '(2x,a,es14.6,a)') 'dphi/dT   = ', dphi_dT, ' 1/K'
  else
    write(*, '(2x,a,i0)') 'ERROR ierr = ', ierr
  end if
  write(*, *)
  write(*, '(a)') '======================================================='

end program test_air_eos_pr
