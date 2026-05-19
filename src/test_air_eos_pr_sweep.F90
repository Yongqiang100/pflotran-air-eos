! =====================================================================
! test_air_eos_pr_sweep
!
! Layer-2 verification: P-T grid sweep producing CSV output suitable
! for plotting and cross-checking against independent EOS sources.
!
! Sweeps:
!   P = 1 - 200 bar in 5-bar steps (40 points)
!   T = 0 - 100 C in 10 C steps (11 points)
!
! Outputs file `air_eos_sweep.csv` with columns:
!   P_Pa, T_K, rho_kg_m3, Z, phi, h_dep_J_mol,
!   drho_dP, drho_dT, dphi_dP, dphi_dT, ierr
!
! Use cases:
!   1. Plot rho(P, T) heatmap and check for monotonic-with-P,
!      smooth-with-T behaviour (no spurious oscillations).
!   2. Plot Z(P, T) and check supercritical Z > 1 trend.
!   3. Compare row-by-row against an independent PR+Peneloux
!      implementation (Python script provided).
!   4. Compare against NIST WebBook isotherms at fixed T.
!
! Also reports max/min/mean values across the grid as a sanity check
! against thermodynamic plausibility.
! =====================================================================
program test_air_eos_pr_sweep

  use Air_EOS_PR_module

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! Grid definition
  integer, parameter :: NP = 40   ! P points
  integer, parameter :: NT = 11   ! T points
  real(kind=dp), parameter :: P_MIN_BAR = 1.0_dp
  real(kind=dp), parameter :: P_MAX_BAR = 200.0_dp
  real(kind=dp), parameter :: T_MIN_C   = 0.0_dp
  real(kind=dp), parameter :: T_MAX_C   = 100.0_dp

  ! Working
  integer       :: ip, it, ierr, unit_csv, n_success, n_fail
  real(kind=dp) :: P, T, dP_bar, dT_C
  real(kind=dp) :: rho, Z, phi, h_dep
  real(kind=dp) :: drho_dP, drho_dT, dphi_dP, dphi_dT
  real(kind=dp) :: rho_min, rho_max, rho_sum
  real(kind=dp) :: Z_min, Z_max
  real(kind=dp) :: phi_min, phi_max
  integer       :: n_total

  write(*, '(a)') '======================================================='
  write(*, '(a)') ' Air PR+Peneloux EOS  P-T grid sweep'
  write(*, '(a,i0,a,i0,a)') ' Grid: ', NP, ' P points x ', NT, ' T points'
  write(*, '(a)') '======================================================='

  dP_bar = (P_MAX_BAR - P_MIN_BAR) / real(NP - 1, dp)
  dT_C   = (T_MAX_C   - T_MIN_C)   / real(NT - 1, dp)

  open(newunit=unit_csv, file='air_eos_sweep.csv', status='replace', action='write')
  write(unit_csv, '(a)')                                                   &
    'P_Pa,T_K,rho_kg_m3,Z,phi,h_dep_J_mol,drho_dP,drho_dT,dphi_dP,dphi_dT,ierr'

  rho_min =  huge(1.0_dp)
  rho_max = -huge(1.0_dp)
  rho_sum =  0.0_dp
  Z_min   =  huge(1.0_dp);  Z_max   = -huge(1.0_dp)
  phi_min =  huge(1.0_dp);  phi_max = -huge(1.0_dp)
  n_success = 0
  n_fail    = 0
  n_total   = 0

  do it = 1, NT
    T = (T_MIN_C + dT_C * real(it - 1, dp)) + 273.15_dp
    do ip = 1, NP
      P = (P_MIN_BAR + dP_bar * real(ip - 1, dp)) * 1.0e5_dp
      n_total = n_total + 1

      call AirEOSPRProperties(P, T, rho, Z, phi, h_dep,                  &
                              drho_dP, drho_dT, dphi_dP, dphi_dT, ierr)

      if (ierr == 0) then
        n_success = n_success + 1
        rho_min = min(rho_min, rho);   rho_max = max(rho_max, rho)
        rho_sum = rho_sum + rho
        Z_min   = min(Z_min,   Z);     Z_max   = max(Z_max,   Z)
        phi_min = min(phi_min, phi);   phi_max = max(phi_max, phi)
        write(unit_csv, '(es15.7,",",es15.7,",",es15.7,",",es15.7,",",   &
                          &es15.7,",",es15.7,",",es15.7,",",es15.7,",",  &
                          &es15.7,",",es15.7,",",i0)')                   &
          P, T, rho, Z, phi, h_dep, drho_dP, drho_dT, dphi_dP, dphi_dT, ierr
      else
        n_fail = n_fail + 1
        write(unit_csv, '(es15.7,",",es15.7,9(",",a),",",i0)')           &
          P, T, 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', ierr
      end if
    end do
  end do

  close(unit_csv)

  ! Summary statistics
  write(*, *)
  write(*, '(a)') '--- Grid summary ---'
  write(*, '(2x,a,i0,a,i0)') 'Success: ', n_success, ' / ', n_total
  write(*, '(2x,a,i0)')      'Failure: ', n_fail
  write(*, *)
  write(*, '(2x,a,2es12.4)') 'rho range [kg/m^3]: ', rho_min, rho_max
  write(*, '(2x,a,es12.4)')  'rho mean  [kg/m^3]: ', rho_sum / real(n_success, dp)
  write(*, '(2x,a,2f10.6)')  'Z   range  [-]    : ', Z_min, Z_max
  write(*, '(2x,a,2f10.6)')  'phi range  [-]    : ', phi_min, phi_max
  write(*, *)
  write(*, '(a)') 'CSV output written to: air_eos_sweep.csv'
  write(*, *)

  ! Quick monotonicity sanity check: rho should be monotonically
  ! increasing with P at every fixed T. Detect any decreases.
  call CheckMonotonicityWithP()

  write(*, '(a)') '======================================================='

contains

  subroutine CheckMonotonicityWithP()
    ! Re-read CSV and check rho(P) is monotone increasing at each T.
    real(kind=dp) :: P_prev, T_curr, T_prev
    real(kind=dp) :: rho_prev, rho_curr
    real(kind=dp) :: P_v, T_v, Z_v, phi_v, h_v
    real(kind=dp) :: drho_dP_v, drho_dT_v, dphi_dP_v, dphi_dT_v
    integer :: ierr_v, n_bad, unit_in, ios

    open(newunit=unit_in, file='air_eos_sweep.csv', status='old', action='read')
    read(unit_in, '(a)') ! header

    P_prev   = -1.0_dp
    T_prev   = -1.0_dp
    rho_prev = -1.0_dp
    n_bad    = 0

    do
      read(unit_in, *, iostat=ios) P_v, T_v, rho_curr, Z_v, phi_v, h_v,   &
        drho_dP_v, drho_dT_v, dphi_dP_v, dphi_dT_v, ierr_v
      if (ios /= 0) exit
      if (ierr_v /= 0) cycle

      T_curr = T_v
      if (abs(T_curr - T_prev) < 1.0e-6_dp) then
        ! Same isotherm: check rho monotone in P
        if (rho_curr <= rho_prev .or. P_v <= P_prev) then
          n_bad = n_bad + 1
          if (n_bad <= 3) then
            write(*, '(2x,a,es12.4,a,f7.2,a)')                            &
              'WARN non-monotone at P=', P_v, ', T=', T_v - 273.15_dp, ' C'
          end if
        end if
      end if
      P_prev   = P_v
      T_prev   = T_curr
      rho_prev = rho_curr
    end do
    close(unit_in)

    write(*, '(2x,a,i0,a)') 'Monotonicity check: ', n_bad,                &
                            ' non-monotone points along isotherms'
  end subroutine CheckMonotonicityWithP

end program test_air_eos_pr_sweep
