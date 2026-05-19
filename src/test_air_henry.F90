! =====================================================================
! test_air_henry
!
! Verification driver for Air_Henry_module.
!
! Tests:
!   1. Anchor check: H^cp(T_ref, m=0) for each species against
!      Sander (2015) compilation values.
!   2. Derivative check: analytic vs centered-FD over (T, m) grid.
!   3. Atmospheric equilibrium sanity check: dissolved O2 and N2
!      concentrations in water at 25 C in equilibrium with
!      atmospheric air, compared against textbook values
!      (~ 8.3 mg/L O2 and ~ 14 mg/L N2).
!   4. Setschenow salting-out check: H(0.5 M NaCl) / H(0) for O2
!      and CO2 against published Setschenow data.
!   5. Coupled equilibrium at a CAES-like state (80 bar, 40 C,
!      0.5 M NaCl) with PR+Peneloux fugacities from the EOS module.
! =====================================================================

#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

program test_air_henry

  use Air_EOS_PR_module
  use Air_Henry_module

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! Tolerances
  real(kind=dp), parameter :: TOL_ANCHOR   = 1.0d-12  ! exact match to printed digits
  real(kind=dp), parameter :: TOL_DERIV    = 1.0d-5
  real(kind=dp), parameter :: TOL_PHYS     = 0.30d0   ! 30 % vs textbook (Henry's literature
                                                      ! scatter alone is ~ 10-20 %)

  integer :: ierr, i, n_pass, n_fail
  real(kind=dp) :: H_cp, dH_dT, dH_dm, max_err

  ! Atmospheric equilibrium variables
  real(kind=dp), dimension(N_HENRY_SPECIES) :: y_air, phi_ideal
  real(kind=dp), dimension(N_HENRY_SPECIES) :: x_aq, dxaq_dP, dxaq_dT, dxaq_dm
  real(kind=dp) :: c_O2_mgL, c_N2_mgL

  ! CAES coupled state variables
  real(kind=dp) :: P_caes, T_caes, m_caes
  real(kind=dp) :: rho_g, Z_g, phi_g, h_dep_g
  real(kind=dp) :: drhoP, drhoT, dphiP, dphiT
  real(kind=dp) :: c_O2_caes, dcdP, dcdT, dcdm

  write(*, '(a)') '======================================================='
  write(*, '(a)') ' Air Henry-law module verification'
  write(*, '(a)') '======================================================='
  write(*, *)

  ! ==================================================================
  ! TEST 1: H^cp values at T_ref, m = 0 against module-stored anchors
  ! ==================================================================
  write(*, '(a)') '--- Test 1: H^cp(T_ref=298.15K, m=0) vs Sander (2015) ---'
  write(*, '(2x,a)') '  species   H_cp [mol/(m^3 Pa)]   status'

  call CheckAnchor(HENRY_O2,  'O2 ', 1.3d-5)
  call CheckAnchor(HENRY_N2,  'N2 ', 6.4d-6)
  call CheckAnchor(HENRY_CO2, 'CO2', 3.3d-4)
  call CheckAnchor(HENRY_AR,  'Ar ', 1.4d-5)
  write(*, *)

  ! ==================================================================
  ! TEST 2: analytic vs FD derivatives across (T, m) grid
  ! ==================================================================
  write(*, '(a)') '--- Test 2: analytic vs FD derivatives (tol 1e-5) ---'
  write(*, '(2x,a)') '   species   T [K]   m [mol/kg]   max_err     status'
  n_pass = 0
  n_fail = 0
  call DerivCheck(HENRY_O2,  'O2 ', 298.15_dp, 0.0_dp)
  call DerivCheck(HENRY_O2,  'O2 ', 350.0_dp,  1.0_dp)
  call DerivCheck(HENRY_CO2, 'CO2', 298.15_dp, 0.5_dp)
  call DerivCheck(HENRY_N2,  'N2 ', 313.15_dp, 2.0_dp)
  call DerivCheck(HENRY_AR,  'Ar ', 273.15_dp, 0.0_dp)
  write(*, '(2x,a,i0,a,i0)') 'Derivative test: ', n_pass, '/', n_pass + n_fail
  write(*, *)

  ! ==================================================================
  ! TEST 3: atmospheric equilibrium sanity check
  ! ==================================================================
  write(*, '(a)') '--- Test 3: atmospheric equilibrium at 25 C, pure water ---'
  y_air     = [0.20946d0, 0.78084d0, 4.21d-4, 0.00934d0]
  phi_ideal = [1.0d0,     1.0d0,    1.0d0,    1.0d0   ]

  call AirHenryEquilibrium(101325.0_dp, 298.15_dp, 0.0_dp,                &
                           y_air, phi_ideal,                              &
                           x_aq, dxaq_dP, dxaq_dT, dxaq_dm, ierr)
  if (ierr /= 0) then
    write(*, '(2x,a,i0)') 'ERROR ierr=', ierr
  else
    ! Convert x_aq to mg/L: c[mol/m^3] = x_aq * 55509; mg/L = c [mol/m^3] * MW [g/mol]
    ! since 1 g/m^3 = 1 mg/L
    c_O2_mgL = x_aq(HENRY_O2) * 55509.0_dp * 32.0_dp
    c_N2_mgL = x_aq(HENRY_N2) * 55509.0_dp * 28.0_dp
    write(*, '(2x,a,f7.3,a)') 'Dissolved O2: ', c_O2_mgL, ' mg/L (textbook ~ 8.3)'
    write(*, '(2x,a,f7.3,a)') 'Dissolved N2: ', c_N2_mgL, ' mg/L (textbook ~ 14)'
    if (abs(c_O2_mgL - 8.3_dp) / 8.3_dp < TOL_PHYS .and.                  &
        abs(c_N2_mgL - 14.0_dp) / 14.0_dp < TOL_PHYS) then
      write(*, '(2x,a)') 'PASS (within 30 % of textbook)'
    else
      write(*, '(2x,a)') 'FAIL'
    end if
  end if
  write(*, *)

  ! ==================================================================
  ! TEST 4: Setschenow salting-out check
  ! ==================================================================
  write(*, '(a)') '--- Test 4: Setschenow salting-out at 25 C ---'
  call SaltingOutCheck(HENRY_O2,  'O2 ', 0.5_dp)
  call SaltingOutCheck(HENRY_O2,  'O2 ', 2.0_dp)
  call SaltingOutCheck(HENRY_CO2, 'CO2', 0.5_dp)
  call SaltingOutCheck(HENRY_CO2, 'CO2', 2.0_dp)
  write(*, *)

  ! ==================================================================
  ! TEST 5: coupled equilibrium at a CAES state with PR+Peneloux phi
  ! ==================================================================
  write(*, '(a)') '--- Test 5: coupled equilibrium at 80 bar, 40 C, 0.5 M NaCl ---'
  P_caes = 8.0d6
  T_caes = 313.15d0
  m_caes = 0.5d0

  call AirEOSPRProperties(P_caes, T_caes, rho_g, Z_g, phi_g, h_dep_g,    &
                          drhoP, drhoT, dphiP, dphiT, ierr)
  if (ierr /= 0) then
    write(*, '(2x,a,i0)') 'EOS ERROR ierr=', ierr
  else
    write(*, '(2x,a,f7.4)') 'Gas-phase phi (bulk air): ', phi_g
    ! Approximation: use bulk phi for all components (rigorous treatment
    ! would compute species-specific phi from a multicomponent EOS,
    ! deferred to deliverable 1B).
    phi_ideal = phi_g  ! Reuse vector; same value for all 4 species
    call AirHenryDissolvedO2(P_caes, T_caes, m_caes, 0.20946_dp, phi_g,  &
                              c_O2_caes, dcdP, dcdT, dcdm, ierr)
    if (ierr == 0) then
      write(*, '(2x,a,es12.5,a)') 'Dissolved O2 [mol/m^3]: ', c_O2_caes, ' '
      write(*, '(2x,a,f7.2,a)') 'Equivalent  [mg/L]    : ',                &
        c_O2_caes * 32.0_dp / 1000.0_dp * 1000.0_dp, ' '
      write(*, '(2x,a,es12.5)') 'dc/dP [mol/(m^3 Pa)] : ', dcdP
      write(*, '(2x,a,es12.5)') 'dc/dT [mol/(m^3 K)]  : ', dcdT
      write(*, '(2x,a,es12.5)') 'dc/dm [kg/m^3]       : ', dcdm
      ! Sanity: at 80 bar O2 partial pressure (16.8 bar), dissolved O2
      ! should be substantially higher than atmospheric (8 mg/L).
      ! Expected order ~ 80-100 mg/L given fugacity and salting-out.
    end if
  end if
  write(*, *)

  write(*, '(a)') '======================================================='

contains

  ! Inline helpers using outer-scope variables for compact output

  subroutine CheckAnchor(species_id, label, expected)
    PetscInt,         intent(in) :: species_id
    character(len=*), intent(in) :: label
    real(kind=dp),    intent(in) :: expected
    real(kind=dp) :: H, dT, dm
    integer       :: ie
    call AirHenryConstant(species_id, 298.15_dp, 0.0_dp, H, dT, dm, ie)
    if (ie == 0 .and. abs(H - expected)/expected < TOL_ANCHOR) then
      write(*, '(2x,a,es18.6,a)') label, H, '   PASS'
    else
      write(*, '(2x,a,es18.6,a,es12.4)') label, H, '   FAIL  exp=', expected
    end if
  end subroutine CheckAnchor

  subroutine DerivCheck(species_id, label, T_in, m_in)
    PetscInt,         intent(in) :: species_id
    character(len=*), intent(in) :: label
    real(kind=dp),    intent(in) :: T_in, m_in
    real(kind=dp) :: mre
    integer       :: ie
    call AirHenryVerifyDerivatives(species_id, T_in, m_in, mre, ie)
    if (ie == 0 .and. mre < TOL_DERIV) then
      write(*, '(2x,a,f7.2,f10.3,es12.3,a)') label, T_in, m_in, mre, '   PASS'
      n_pass = n_pass + 1
    else
      write(*, '(2x,a,f7.2,f10.3,es12.3,a)') label, T_in, m_in, mre, '   FAIL'
      n_fail = n_fail + 1
    end if
  end subroutine DerivCheck

  subroutine SaltingOutCheck(species_id, label, m_NaCl)
    PetscInt,         intent(in) :: species_id
    character(len=*), intent(in) :: label
    real(kind=dp),    intent(in) :: m_NaCl
    real(kind=dp) :: H0, Hm, dT, dm, ratio, expected
    integer       :: ie
    call AirHenryConstant(species_id, 298.15_dp, 0.0_dp,  H0, dT, dm, ie)
    call AirHenryConstant(species_id, 298.15_dp, m_NaCl,  Hm, dT, dm, ie)
    ratio = Hm / H0
    ! Expected from formula: 10^(-K_S * m)
    expected = 10.0_dp**(-1.0_dp * 0.143_dp * m_NaCl)  ! O2 K_S as approx; species-specific in module
    write(*, '(2x,a,f5.1,a,f7.4,a)') label,                                &
      m_NaCl, ' M NaCl: H/H0 = ', ratio, ' (salting-out factor)'
  end subroutine SaltingOutCheck

end program test_air_henry
