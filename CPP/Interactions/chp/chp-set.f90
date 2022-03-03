
! Preset for N2LO opt
subroutine chp_set_N2LOopt(parameter_array) bind(C,name="chp_set_N2LOopt")
  real (8), dimension (1:15), intent (in) :: parameter_array
  call chp_preset_NNLOopt(parameter_array)
end subroutine chp_set_N2LOopt

! Preset for Idaho N3LO
subroutine chp_set_Idaho_N3LO(parameter_array) bind(C,name="chp_set_Idaho_N3LO")
  real (8), dimension (1:35), intent (in) :: parameter_array
  call chp_preset_Idaho_N3LO(parameter_array)
end subroutine chp_set_Idaho_N3LO

! Preset for IS LO
subroutine chp_set_ispot_lo(parameter_array) bind(C,name="chp_set_ispot_lo")
  real (8), dimension (1:2), intent (in) :: parameter_array
  call chp_preset_ispot_lo(parameter_array)
end subroutine chp_set_ispot_lo

! Preset for IS NLO
subroutine chp_set_ispot_nlo(parameter_array) bind(C,name="chp_set_ispot_nlo")
  real (8), dimension (1:13), intent (in) :: parameter_array
  call chp_preset_ispot_nlo(parameter_array)
end subroutine chp_set_ispot_nlo

! Preset for IS N2LO
subroutine chp_set_ispot_n2lo(parameter_array) bind(C,name="chp_set_ispot_n2lo")
  real (8), dimension (1:16), intent (in) :: parameter_array
  call chp_preset_ispot_n2lo(parameter_array)
end subroutine chp_set_ispot_n2lo

! ############################################################################################################# !
! ############################################################################################################# !
! ############################################################################################################# !
! LO chiral potential
subroutine chp_preset_LO
  use idaho_chiral_potential
  
  implicit none
  
  ! real(8) :: p,pp, ME
  ! integer :: S,j,T,tz, n, l ,nn, ll
  ! real(8) :: pot(1:6)
   character(len=42) :: title
  ! logical :: coup
  
  ! p   = 0.2
  ! pp  = 0.3
  ! coup = .true.
  ! S = 1
  ! j = 3
  ! T = 1
  ! tz = 0

  call initialize_chiral_potential
  
  ! GENERAL PARAMETERS AND CONSTANTS
  ! (proton, nucleon, neutron)
  call chp_set_mass_nucleon((/9.38272000000000048D+02, 9.38918054847146095D+02, 9.39565000000000055D+02/))
  ! (pi-, pi, pi+)
  call chp_set_mass_pion((/1.38038999999999987D+02, 1.38038999999999987D+02, 1.38038999999999987D+02/))
  
  call chp_set_chiral_order(LO)
  call chp_set_chiral_mode(chiral_mode_EM2015)
  call chp_set_reg("SF", 7.00000000000000000D+02)
  call chp_set_gA(1.28899999999999992D+00)
  call chp_set_fpi(9.22000000000000028D+01)
  call chp_set_fine_structure(7.29735257000000033D-03)
  
  call chp_set_Lambda(4.50000000000000000D+02)
  call chp_set_contact_format("PW")
  
  ! SET INCLUDED DIAGRAMS / TERMS
  
  ! LO contacts
  call chp_set_chiral_Ct(1) ! Use non-CIB contacts
  call chp_set_LO_contact(1, -1.12926999999999986D-01) ! Ct_1S0
  call chp_set_LO_contact(2, -8.73400000000000010D-02) ! Ct_3S1
  
  ! NLO contacts
  call chp_set_chiral_C(0) ! Do not use
  
  ! N3LO contacts
  call chp_set_chiral_D(0) ! Do not use
  
  ! Set needed ci and di, if any
  
  ! Set regulator parameter n
  call chp_set_1PE_reg_par(3.0D0)
  
  ! Set pion exchange contributions
  ! Basic 1PE
  call chp_set_chiral_1PE(1) ! Use
  ! CIB effects in 1PE
  call chp_set_chiral_1PE_CIB(0) ! Do not use
  ! pion-gamma exchange
  call chp_set_chiral_1PE_gamma(0) ! Do not use
  ! Relativistic corrections to 1PE
  call chp_set_chiral_1PE_relcorr(0) ! Do not use
  
  ! Leading 2PE
  call chp_set_chiral_2PE_1loop_0(0) ! Do not use
  
  ! 1-loop 2PE proportional to ci
  call chp_set_chiral_2PE_1loop_d(0) ! Do not use
  
  ! 1-loop 2PE proportional to 1/M_N (relativistic corrections)
  call chp_set_chiral_2PE_1loop_r(0) ! Do not use
  
  ! 1-loop 2PE proportional to ci*cj)
  call chp_set_chiral_2PE_1loop_dd(0) ! Do not use
  
  ! 1-loop 2PE proportional to ci/M_N)
  call chp_set_chiral_2PE_1loop_dr(0) ! Do not use
  
  ! 1-loop 2PE proportional to 1/M_N^2)
  call chp_set_chiral_2PE_1loop_rr(0) ! Do not use
  
  ! 2-loop 2PE
  call chp_set_chiral_2PE_2loop(0) ! Do not use
  ! Contributions to 2-loop 2PE that do not have analytical expressions
  call chp_set_chiral_2PE_2loop_int(0) ! Do not use
  
  ! Use correct nucleon mass in 2PE relativistic corrections
  call chp_set_chiral_2PE_CSB_correct_mass(0) ! Do not use
  
  ! Use minimal relativity
  call chp_set_chiral_minimal_relativity(1) ! Use
  
  ! Use Kamada-Glockle transform
  call chp_set_chiral_kamada_glockle_transform(0) ! Do not use
  
  call chp_set_units_and_derive_constants
  
  !call chp_print_constants(title,6)
  !call chp(p,pp,coup,S,j,T,tz,pot)
  !write(6,*) p
  !write(6,*) pp
  !write(6,*) pot
end subroutine chp_preset_LO

! ############################################################################################################# !
! ############################################################################################################# !
! ############################################################################################################# !
! Idaho N2LO chiral potential
subroutine chp_preset_NNLOopt(parameter_array)
  use idaho_chiral_potential

  implicit none
  
  real (8), dimension (0:14), intent (in) :: parameter_array
  
  real (8) :: Ct_1S0pp
  real (8) :: Ct_3S1pp
  real (8) :: Ct_1S0np
  real (8) :: Ct_3S1np
  real (8) :: Ct_1S0nn
  real (8) :: Ct_3S1nn
  real (8) :: C_1S0
  real (8) :: C_3P0
  real (8) :: C_1P1
  real (8) :: C_3P1
  real (8) :: C_3S1
  real (8) :: C_3S1_3D1
  real (8) :: C_3P2
  real (8) :: gA
  real (8) :: c1
  real (8) :: c3
  real (8) :: c4
  
  ! LO contact terms
  Ct_1S0pp = parameter_array(0)
  Ct_1S0np = parameter_array(1)
  Ct_1S0nn = parameter_array(2)
  ! Ct_3S1np/pp/nn must be equal
  Ct_3S1pp = parameter_array(3)
  Ct_3S1np = parameter_array(3)
  Ct_3S1nn = parameter_array(3)
  
  ! NLO contact terms
  C_1S0 = parameter_array(4)
  C_3P0 = parameter_array(5)
  C_1P1 = parameter_array(6)
  C_3P1 = parameter_array(7)
  C_3S1 = parameter_array(8)
  C_3S1_3D1 = parameter_array(9)
  C_3P2 = parameter_array(10)
  
  ! Other constants
  gA = parameter_array(11)
  c1 = parameter_array(12)
  c3 = parameter_array(13)
  c4 = parameter_array(14)
  
  call initialize_chiral_potential

  ! GENERAL PARAMETERS AND CONSTANTS
  ! (proton, nucleon, neutron)
  call chp_set_mass_nucleon((/9.38272000000000048D+02, 9.38918400000000020D+02, 9.39565299999999979D+02/))
  ! (pi-, pi, pi+)
  call chp_set_mass_pion((/1.39570200000000000D+02, 1.34976599999999991D+02, 1.39570200000000000D+02/))

  call chp_set_chiral_order(NNLO)
  call chp_set_chiral_mode(chiral_mode_EM2011)
  call chp_set_reg("SF", 7.00000000000000000D+02)
  call chp_set_gA(gA)
  call chp_set_fpi(9.24000000000000057D+01)
  call chp_set_fine_structure(7.29735256979999972D-03)

  call chp_set_Lambda(5.00000000000000000D+02)
  call chp_set_contact_format("PW")

  ! SET INCLUDED DIAGRAMS / TERMS

  ! LO contacts
  call chp_set_chiral_Ct_CIB(1) ! Use CIB contacts
  call chp_set_CIB_LO_contact(1, -1, Ct_1S0pp) ! Ct_1S0pp
  call chp_set_CIB_LO_contact(2, -1, Ct_3S1pp) ! Ct_3S1pp
  call chp_set_CIB_LO_contact(1,  0, Ct_1S0np) ! Ct_1S0np
  call chp_set_CIB_LO_contact(2,  0, Ct_3S1np) ! Ct_3S1np
  call chp_set_CIB_LO_contact(1,  1, Ct_1S0nn) ! Ct_1S0nn
  call chp_set_CIB_LO_contact(2,  1, Ct_3S1nn) ! Ct_3S1nn

  ! NLO contacts
  call chp_set_chiral_C(1) ! Use
  call chp_set_NLO_contact(1, C_1S0) ! C_1S0
  call chp_set_NLO_contact(2, C_3P0) ! C_3P0
  call chp_set_NLO_contact(3, C_1P1) ! C_1P1
  call chp_set_NLO_contact(4, C_3P1) ! C_3P1
  call chp_set_NLO_contact(5, C_3S1) ! C_3S1
  call chp_set_NLO_contact(6, C_3S1_3D1) ! C_3S1-3D1
  call chp_set_NLO_contact(7, C_3P2) ! C_3P2

  ! N3LO contacts
  call chp_set_chiral_D(0) ! Do not use

  ! Set needed ci and di, if any
  call chp_set_c1(c1)
  call chp_set_c3(c3)
  call chp_set_c4(c4)

  ! Parameters for the NNN force
  !call chp_set_c_D(-2.00000000000000011D-01)
  !call chp_set_c_E(-3.59999999999999987D-01)

  ! Set regulator parameter n
  call chp_set_1PE_reg_par(3.0D0)
  call chp_set_2PE_reg_par(3.0D0)
  call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
  call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
  call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
  call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
  call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
  call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
  call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
  call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
  call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2
  call chp_set_N3LO_contact_reg_par(1, 3.0D0) ! Dh_1S0
  call chp_set_N3LO_contact_reg_par(2, 3.0D0) ! D_1S0
  call chp_set_N3LO_contact_reg_par(3, 3.0D0) ! D_3P0
  call chp_set_N3LO_contact_reg_par(4, 3.0D0) ! D_1P1
  call chp_set_N3LO_contact_reg_par(5, 3.0D0) ! D_3P1
  call chp_set_N3LO_contact_reg_par(6, 3.0D0) ! Dh_3S1
  call chp_set_N3LO_contact_reg_par(7, 3.0D0) ! D_3S1
  call chp_set_N3LO_contact_reg_par(8, 3.0D0) ! D_3D1
  call chp_set_N3LO_contact_reg_par(9, 3.0D0) ! Dh_3S1-3D1
  call chp_set_N3LO_contact_reg_par(10, 3.0D0) ! D_3S1-3D1
  call chp_set_N3LO_contact_reg_par(11, 3.0D0) ! D_1D2
  call chp_set_N3LO_contact_reg_par(12, 3.0D0) ! D_3D2
  call chp_set_N3LO_contact_reg_par(13, 3.0D0) ! D_3P2
  call chp_set_N3LO_contact_reg_par(14, 3.0D0) ! D_3P2-3F2
  call chp_set_N3LO_contact_reg_par(15, 3.0D0) ! D_3D3

  ! Set pion exchange contributions
  ! Basic 1PE
  call chp_set_chiral_1PE(1) ! Use
  ! CIB effects in 1PE
  call chp_set_chiral_1PE_CIB(1) ! Use
  ! pion-gamma exchange
  call chp_set_chiral_1PE_gamma(0) ! Do not use
  ! Relativistic corrections to 1PE
  call chp_set_chiral_1PE_relcorr(0) ! Do not use

  ! Leading 2PE
  call chp_set_chiral_2PE_1loop_0(1) ! Use

  ! 1-loop 2PE proportional to ci
  call chp_set_chiral_2PE_1loop_d(1) ! Use

  ! 1-loop 2PE proportional to 1/M_N (relativistic corrections)
  call chp_set_chiral_2PE_1loop_r(1) ! Use
  call chp_set_chiral_2PE_1loop_r_mode(201) ! mode 'EM 2011'

  ! 1-loop 2PE proportional to ci*cj)
  call chp_set_chiral_2PE_1loop_dd(0) ! Do not use

  ! 1-loop 2PE proportional to ci/M_N)
  call chp_set_chiral_2PE_1loop_dr(0) ! Do not use

  ! 1-loop 2PE proportional to 1/M_N^2)
  call chp_set_chiral_2PE_1loop_rr(0) ! Do not use

  ! 2-loop 2PE
  call chp_set_chiral_2PE_2loop(0) ! Do not use
  ! Contributions to 2-loop 2PE that do not have analytical expressions
  call chp_set_chiral_2PE_2loop_int(0) ! Do not use

  ! Use correct nucleon mass in 2PE relativistic corrections
  call chp_set_chiral_2PE_CSB_correct_mass(1) ! Use

  ! Use minimal relativity
  call chp_set_chiral_minimal_relativity(1) ! Use

  ! Use Kamada-Glockle transform
  call chp_set_chiral_kamada_glockle_transform(0) ! Do not use

  call chp_set_units_and_derive_constants
end subroutine

! ############################################################################################################# !
! ############################################################################################################# !
! ############################################################################################################# !
! Idaho N3LO chiral potential
subroutine chp_preset_Idaho_N3LO(parameter_array)
  use idaho_chiral_potential
  
  implicit none
  
  character(len=42) :: title
  
  real (8), dimension (0:34), intent (in) :: parameter_array
  ! logical :: use_parameter_array
  
  real (8) :: Ct_1S0pp
  real (8) :: Ct_3S1pp
  real (8) :: Ct_1S0np
  real (8) :: Ct_3S1np
  real (8) :: Ct_1S0nn
  real (8) :: Ct_3S1nn
  real (8) :: C_1S0
  real (8) :: C_3P0
  real (8) :: C_1P1
  real (8) :: C_3P1
  real (8) :: C_3S1
  real (8) :: C_3S1_3D1
  real (8) :: C_3P2
  real (8) :: Dh_1S0
  real (8) :: D_1S0
  real (8) :: D_3P0
  real (8) :: D_1P1
  real (8) :: D_3P1
  real (8) :: Dh_3S1
  real (8) :: D_3S1
  real (8) :: D_3D1
  real (8) :: Dh_3S1_3D1
  real (8) :: D_3S1_3D1
  real (8) :: D_1D2
  real (8) :: D_3D2
  real (8) :: D_3P2
  real (8) :: D_3P2_3F2
  real (8) :: D_3D3
  real (8) :: gA
  real (8) :: c1
  real (8) :: c2
  real (8) :: c3
  real (8) :: c4
  real (8) :: d1_plus_d2
  real (8) :: d3
  real (8) :: d5
  real (8) :: d14_minus_d15
  
  title = "Idaho-N3LO"
  
  ! if use_parameter_array==1 then
  ! LO contact terms
  Ct_1S0pp = parameter_array(0)
  Ct_1S0np = parameter_array(1)
  Ct_1S0nn = parameter_array(2)
  ! Ct_3S1np/pp/nn must be equal
  Ct_3S1pp = parameter_array(3)
  Ct_3S1np = parameter_array(3)
  Ct_3S1nn = parameter_array(3)
  
  ! NLO contact terms
  C_1S0 = parameter_array(4)
  C_3P0 = parameter_array(5)
  C_1P1 = parameter_array(6)
  C_3P1 = parameter_array(7)
  C_3S1 = parameter_array(8)
  C_3S1_3D1 = parameter_array(9)
  C_3P2 = parameter_array(10)
  
  ! N3LO contact terms
  Dh_1S0 = parameter_array(11)
  D_1S0 = parameter_array(12)
  D_3P0 = parameter_array(13)
  D_1P1 = parameter_array(14)
  D_3P1 = parameter_array(15)
  Dh_3S1 = parameter_array(16)
  D_3S1 = parameter_array(17)
  D_3D1 = parameter_array(18)
  Dh_3S1_3D1 = parameter_array(19)
  D_3S1_3D1 = parameter_array(20)
  D_1D2 = parameter_array(21)
  D_3D2 = parameter_array(22)
  D_3P2 = parameter_array(23)
  D_3P2_3F2 = parameter_array(24)
  D_3D3 = parameter_array(25)
  
  ! Other constants
  gA = parameter_array(26)
  c1 = parameter_array(27)
  c2 = parameter_array(28)
  c3 = parameter_array(29)
  c4 = parameter_array(30)
  d1_plus_d2 = parameter_array(31)
  d3 = parameter_array(32)
  d5 = parameter_array(33)
  d14_minus_d15 = parameter_array(34)
  
  call initialize_chiral_potential
  
  ! GENERAL PARAMETERS AND CONSTANTS
  ! (proton, nucleon, neutron)
  call chp_set_mass_nucleon((/9.38272000000000048D+02, 9.38918204640625618D+02, 9.39565299999999979D+02/))
  ! (pi-, pi, pi+)
  call chp_set_mass_pion((/1.39570200000000000D+02, 1.34976599999999991D+02, 1.39570200000000000D+02/))
  
  call chp_set_chiral_order(N3LO)
  call chp_set_chiral_mode(chiral_mode_EM2011)
  call chp_set_reg("DR", 0.00000000000000000D+00)
  call chp_set_gA(gA)
  call chp_set_fpi(9.24000000000000057D+01)
  call chp_set_fine_structure(7.29679570071830198D-03)
  
  call chp_set_Lambda(5.00000000000000000D+02)
  call chp_set_contact_format("PW")
  
  ! SET INCLUDED DIAGRAMS / TERMS
  
  ! LO contacts
  call chp_set_chiral_Ct_CIB(1) ! Use CIB contacts
  call chp_set_CIB_LO_contact(1, -1, Ct_1S0pp) ! Ct_1S0pp
  call chp_set_CIB_LO_contact(2, -1, Ct_3S1pp) ! Ct_3S1pp
  call chp_set_CIB_LO_contact(1,  0, Ct_1S0np) ! Ct_1S0np
  call chp_set_CIB_LO_contact(2,  0, Ct_3S1np) ! Ct_3S1np
  call chp_set_CIB_LO_contact(1,  1, Ct_1S0nn) ! Ct_1S0nn
  call chp_set_CIB_LO_contact(2,  1, Ct_3S1nn) ! Ct_3S1nn
  
  ! NLO contacts
  call chp_set_chiral_C(1) ! Use
  call chp_set_NLO_contact(1, C_1S0) ! C_1S0
  call chp_set_NLO_contact(2, C_3P0)! C_3P0
  call chp_set_NLO_contact(3, C_1P1)! C_1P1
  call chp_set_NLO_contact(4, C_3P1) ! C_3P1
  call chp_set_NLO_contact(5, C_3S1)! C_3S1
  call chp_set_NLO_contact(6, C_3S1_3D1)! C_3S1-3D1
  call chp_set_NLO_contact(7, C_3P2) ! C_3P2
  
  ! N3LO contacts
  call chp_set_chiral_D(1) ! Use
  call chp_set_N3LO_contact( 1, Dh_1S0)  ! Dh_1S0
  call chp_set_N3LO_contact( 2, D_1S0)  ! D_1S0
  call chp_set_N3LO_contact( 3, D_3P0) ! D_3P0
  call chp_set_N3LO_contact( 4, D_1P1) ! D_1P1
  call chp_set_N3LO_contact( 5, D_3P1) ! D_3P1
  call chp_set_N3LO_contact( 6, Dh_3S1) ! Dh_3S1
  call chp_set_N3LO_contact( 7, D_3S1) ! D_3S1
  call chp_set_N3LO_contact( 8, D_3D1)  ! D_3D1
  call chp_set_N3LO_contact( 9, Dh_3S1_3D1) ! Dh_3S1-3D1
  call chp_set_N3LO_contact(10, D_3S1_3D1) ! D_3S1-3D1
  call chp_set_N3LO_contact(11, D_1D2)  ! D_1D2
  call chp_set_N3LO_contact(12, D_3D2)  ! D_3D2
  call chp_set_N3LO_contact(13, D_3P2) ! D_3P2
  call chp_set_N3LO_contact(14, D_3P2_3F2)  ! D_3P2-3F2
  call chp_set_N3LO_contact(15, D_3D3) ! D_3D3
  
  ! Set needed ci and di, if any
  call chp_set_c1(c1)
  call chp_set_c2(c2)
  call chp_set_c3(c3)
  call chp_set_c4(c4)
  call chp_set_d1_plus_d2(d1_plus_d2)
  call chp_set_d3(d3)
  call chp_set_d5(d5)
  call chp_set_d14_minus_d15(d14_minus_d15)
  
  ! Set regulator parameter n
  call chp_set_1PE_reg_par(4.0D0)
  call chp_set_2PE_reg_par(2.0D0)
  call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
  call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
  call chp_set_NLO_contact_reg_par(1, 2.0D0) ! C_1S0
  call chp_set_NLO_contact_reg_par(2, 2.0D0) ! C_3P0
  call chp_set_NLO_contact_reg_par(3, 2.0D0) ! C_1P1
  call chp_set_NLO_contact_reg_par(4, 2.0D0) ! C_3P1
  call chp_set_NLO_contact_reg_par(5, 2.0D0) ! C_3S1
  call chp_set_NLO_contact_reg_par(6, 2.0D0) ! C_3S1-3D1
  call chp_set_NLO_contact_reg_par(7, 2.0D0) ! C_3P2
  call chp_set_N3LO_contact_reg_par(1, 2.0D0) ! Dh_1S0
  call chp_set_N3LO_contact_reg_par(2, 2.0D0) ! D_1S0
  call chp_set_N3LO_contact_reg_par(3, 3.0D0) ! D_3P0
  call chp_set_N3LO_contact_reg_par(4, 2.0D0) ! D_1P1
  call chp_set_N3LO_contact_reg_par(5, 4.0D0) ! D_3P1
  call chp_set_N3LO_contact_reg_par(6, 2.0D0) ! Dh_3S1
  call chp_set_N3LO_contact_reg_par(7, 2.0D0) ! D_3S1
  call chp_set_N3LO_contact_reg_par(8, 2.0D0) ! D_3D1
  call chp_set_N3LO_contact_reg_par(9, 2.0D0) ! Dh_3S1-3D1
  call chp_set_N3LO_contact_reg_par(10, 2.0D0) ! D_3S1-3D1
  call chp_set_N3LO_contact_reg_par(11, 4.0D0) ! D_1D2
  call chp_set_N3LO_contact_reg_par(12, 2.0D0) ! D_3D2
  call chp_set_N3LO_contact_reg_par(13, 2.0D0) ! D_3P2
  call chp_set_N3LO_contact_reg_par(14, 4.0D0) ! D_3P2-3F2
  call chp_set_N3LO_contact_reg_par(15, -1.0D0) ! D_3D3
  
  ! Set pion exchange contributions
  ! Basic 1PE
  call chp_set_chiral_1PE(1) ! Use
  ! CIB effects in 1PE
  call chp_set_chiral_1PE_CIB(1) ! Use
  ! pion-gamma exchange
  call chp_set_chiral_1PE_gamma(1) ! Use
  ! Relativistic corrections to 1PE
  call chp_set_chiral_1PE_relcorr(0) ! Do not use
  
  ! Leading 2PE
  call chp_set_chiral_2PE_1loop_0(1) ! Use
  
  ! 1-loop 2PE proportional to ci
  call chp_set_chiral_2PE_1loop_d(1) ! Use
  
  ! 1-loop 2PE proportional to 1/M_N (relativistic corrections)
  call chp_set_chiral_2PE_1loop_r(1) ! Use
  call chp_set_chiral_2PE_1loop_r_mode(201) ! mode 'EM 2011'
  
  ! 1-loop 2PE proportional to ci*cj)
  call chp_set_chiral_2PE_1loop_dd(1) ! Use
  
  ! 1-loop 2PE proportional to ci/M_N)
  call chp_set_chiral_2PE_1loop_dr(1) ! Use
  
  ! 1-loop 2PE proportional to 1/M_N^2)
  call chp_set_chiral_2PE_1loop_rr(1) ! Use
  
  ! 2-loop 2PE
  call chp_set_chiral_2PE_2loop(1) ! Use
  ! Contributions to 2-loop 2PE that do not have analytical expressions
  call chp_set_chiral_2PE_2loop_int(0) ! Do not use
  
  ! Use correct nucleon mass in 2PE relativistic corrections
  call chp_set_chiral_2PE_CSB_correct_mass(1) ! Use
  
  ! Use minimal relativity
  call chp_set_chiral_minimal_relativity(1) ! Use
  
  ! Use Kamada-Glockle transform
  call chp_set_chiral_kamada_glockle_transform(0) ! Do not use
  
  call chp_set_units_and_derive_constants
  
  !call chp_print_constants(title,6)
end subroutine chp_preset_Idaho_N3LO

! ############################################################################################################# !
! ############################################################################################################# !
! ############################################################################################################# !
! I. Svensson's leading-order potential (Phys. Rev. C 105, 014004 (2022))
subroutine chp_preset_ispot_lo(parameter_array)
  use idaho_chiral_potential

  implicit none

  real (8), dimension (0:1), intent (in) :: parameter_array
  
  ! LO contacts
  real (8) :: Ct_1S0
  real (8) :: Ct_3S1

  !if (use_parameter_array) then
  !  ! LO contact terms
  !  Ct_1S0 = parameter_array(0)
  !  Ct_3S1 = parameter_array(1)
  !else
  Ct_1S0 = -1.12073397972979993D-01
  Ct_3S1 = -3.60868520937899992D-02
  !end if

  call initialize_chiral_potential

  ! GENERAL PARAMETERS AND CONSTANTS
  ! (proton, nucleon, neutron)
  call chp_set_mass_nucleon((/9.38272046000000046D+02, 9.38918267117927371D+02, 9.39565379000000007D+02/))
  ! (pi-, pi, pi+)
  call chp_set_mass_pion((/1.39570179999999993D+02, 1.34976599999999991D+02, 1.39570179999999993D+02/))

  call chp_set_chiral_order(LO)
  call chp_set_chiral_mode(chiral_mode_EM2015)
  call chp_set_reg("SF", 7.00000000000000000D+02)
  call chp_set_gA(1.28899999999999992D+00)
  call chp_set_fpi(9.22000000000000028D+01)
  call chp_set_fine_structure(7.29735256899999990D-03)

  call chp_set_Lambda(4.50000000000000000D+02)
  call chp_set_contact_format("PW")

  ! SET INCLUDED DIAGRAMS / TERMS

  ! LO contacts
  call chp_set_chiral_Ct(1) ! Use non-CIB contacts
  call chp_set_LO_contact(1, Ct_1S0) ! Ct_1S0
  call chp_set_LO_contact(2, Ct_3S1) ! Ct_3S1

  ! NLO contacts
  call chp_set_chiral_C(0) ! Do not use

  ! N3LO contacts
  call chp_set_chiral_D(0) ! Do not use

  ! Set needed ci and di, if any

  ! Set regulator parameter n
  call chp_set_1PE_reg_par(3.0D0)
  call chp_set_2PE_reg_par(3.0D0)
  call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
  call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
  call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
  call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
  call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
  call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
  call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
  call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
  call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2
  call chp_set_N3LO_contact_reg_par(1, 3.0D0) ! Dh_1S0
  call chp_set_N3LO_contact_reg_par(2, 3.0D0) ! D_1S0
  call chp_set_N3LO_contact_reg_par(3, 3.0D0) ! D_3P0
  call chp_set_N3LO_contact_reg_par(4, 3.0D0) ! D_1P1
  call chp_set_N3LO_contact_reg_par(5, 3.0D0) ! D_3P1
  call chp_set_N3LO_contact_reg_par(6, 3.0D0) ! Dh_3S1
  call chp_set_N3LO_contact_reg_par(7, 3.0D0) ! D_3S1
  call chp_set_N3LO_contact_reg_par(8, 3.0D0) ! D_3D1
  call chp_set_N3LO_contact_reg_par(9, 3.0D0) ! Dh_3S1-3D1
  call chp_set_N3LO_contact_reg_par(10, 3.0D0) ! D_3S1-3D1
  call chp_set_N3LO_contact_reg_par(11, 3.0D0) ! D_1D2
  call chp_set_N3LO_contact_reg_par(12, 3.0D0) ! D_3D2
  call chp_set_N3LO_contact_reg_par(13, 3.0D0) ! D_3P2
  call chp_set_N3LO_contact_reg_par(14, 3.0D0) ! D_3P2-3F2
  call chp_set_N3LO_contact_reg_par(15, 3.0D0) ! D_3D3

  ! Set pion exchange contributions
  ! Basic 1PE
  call chp_set_chiral_1PE(1) ! Use
  ! CIB effects in 1PE
  call chp_set_chiral_1PE_CIB(0) ! Do not use
  ! pion-gamma exchange
  call chp_set_chiral_1PE_gamma(0) ! Do not use
  ! Relativistic corrections to 1PE
  call chp_set_chiral_1PE_relcorr(0) ! Do not use

  ! Leading 2PE
  call chp_set_chiral_2PE_1loop_0(0) ! Do not use

  ! 1-loop 2PE proportional to ci
  call chp_set_chiral_2PE_1loop_d(0) ! Do not use

  ! 1-loop 2PE proportional to 1/M_N (relativistic corrections)
  call chp_set_chiral_2PE_1loop_r(0) ! Do not use

  ! 1-loop 2PE proportional to ci*cj)
  call chp_set_chiral_2PE_1loop_dd(0) ! Do not use

  ! 1-loop 2PE proportional to ci/M_N)
  call chp_set_chiral_2PE_1loop_dr(0) ! Do not use

  ! 1-loop 2PE proportional to 1/M_N^2)
  call chp_set_chiral_2PE_1loop_rr(0) ! Do not use

  ! 2-loop 2PE
  call chp_set_chiral_2PE_2loop(0) ! Do not use
  ! Contributions to 2-loop 2PE that do not have analytical expressions
  call chp_set_chiral_2PE_2loop_int(0) ! Do not use

  ! Use correct nucleon mass in 2PE relativistic corrections
  call chp_set_chiral_2PE_CSB_correct_mass(0) ! Do not use

  ! Use minimal relativity
  call chp_set_chiral_minimal_relativity(1) ! Use

  ! Use Kamada-Glockle transform
  call chp_set_chiral_kamada_glockle_transform(0) ! Do not use

  call chp_set_units_and_derive_constants

end subroutine chp_preset_ispot_lo

! ############################################################################################################# !
! ############################################################################################################# !
! ############################################################################################################# !
! I. Svensson's nex-to-leading-order potential (Phys. Rev. C 105, 014004 (2022))
subroutine chp_preset_ispot_nlo(parameter_array)
  use idaho_chiral_potential

  implicit none

  real (8), dimension (0:12), intent (in) :: parameter_array
  
  ! LO contacts
  real (8) :: Ct_1S0pp
  real (8) :: Ct_3S1pp
  real (8) :: Ct_1S0np
  real (8) :: Ct_3S1np
  real (8) :: Ct_1S0nn
  real (8) :: Ct_3S1nn
  ! NLO contacts
  real (8) :: C_1S0
  real (8) :: C_3P0
  real (8) :: C_1P1
  real (8) :: C_3P1
  real (8) :: C_3S1
  real (8) :: C_3S1_3D1
  real (8) :: C_3P2

  !if (use_parameter_array) then
  Ct_1S0pp  = parameter_array(0)
  Ct_3S1pp  = parameter_array(1)
  Ct_1S0np  = parameter_array(2)
  Ct_3S1np  = parameter_array(3)
  Ct_1S0nn  = parameter_array(4)
  Ct_3S1nn  = parameter_array(5)
  C_1S0     = parameter_array(6)
  C_3P0     = parameter_array(7)
  C_1P1     = parameter_array(8)
  C_3P1     = parameter_array(9)
  C_3S1     = parameter_array(10)
  C_3S1_3D1 = parameter_array(11)
  C_3P2     = parameter_array(12)
  !else
  !  Ct_1S0pp  = -1.54250208697730012D-01
  !  Ct_3S1pp  = -1.48057121689049997D-01
  !  Ct_1S0np  = -1.55222568956019991D-01
  !  Ct_3S1np  = -1.48057121689049997D-01
  !  Ct_1S0nn  = -1.54958402671070000D-01
  !  Ct_3S1nn  = -1.48057121689049997D-01
  !  C_1S0     =  1.60603274119125006D+00
  !  C_3P0     =  1.19584298162829006D+00
  !  C_1P1     =  6.06776291764649978D-01
  !  C_3P1     = -2.49343337260960030D-01
  !  C_3S1     = -7.46664535867350043D-01
  !  C_3S1_3D1 =  1.56125611619209986D-01
  !  C_3P2     = -1.86669728671579993D-01
  !end if

  call initialize_chiral_potential

  ! GENERAL PARAMETERS AND CONSTANTS
  ! (proton, nucleon, neutron)
  call chp_set_mass_nucleon((/9.38272046000000046D+02, 9.38918267117927371D+02, 9.39565379000000007D+02/))
  ! (pi-, pi, pi+)
  call chp_set_mass_pion((/1.39570179999999993D+02, 1.34976599999999991D+02, 1.39570179999999993D+02/))

  call chp_set_chiral_order(NLO)
  call chp_set_chiral_mode(chiral_mode_EM2015)
  call chp_set_reg("SF", 7.00000000000000000D+02)
  call chp_set_gA(1.28899999999999992D+00)
  call chp_set_fpi(9.22000000000000028D+01)
  call chp_set_fine_structure(7.29735256899999990D-03)

  call chp_set_Lambda(4.50000000000000000D+02)
  call chp_set_contact_format("PW")

  ! SET INCLUDED DIAGRAMS / TERMS

  ! LO contacts
  call chp_set_chiral_Ct_CIB(1) ! Use CIB contacts
  call chp_set_CIB_LO_contact(1, -1, Ct_1S0pp) ! Ct_1S0pp
  call chp_set_CIB_LO_contact(2, -1, Ct_3S1pp) ! Ct_3S1pp
  call chp_set_CIB_LO_contact(1,  0, Ct_1S0np) ! Ct_1S0np
  call chp_set_CIB_LO_contact(2,  0, Ct_3S1np) ! Ct_3S1np
  call chp_set_CIB_LO_contact(1,  1, Ct_1S0nn) ! Ct_1S0nn
  call chp_set_CIB_LO_contact(2,  1, Ct_3S1nn) ! Ct_3S1nn

  ! NLO contacts
  call chp_set_chiral_C(1) ! Use
  call chp_set_NLO_contact(1, C_1S0) ! C_1S0
  call chp_set_NLO_contact(2, C_3P0) ! C_3P0
  call chp_set_NLO_contact(3, C_1P1) ! C_1P1
  call chp_set_NLO_contact(4, C_3P1) ! C_3P1
  call chp_set_NLO_contact(5, C_3S1) ! C_3S1
  call chp_set_NLO_contact(6, C_3S1_3D1) ! C_3S1-3D1
  call chp_set_NLO_contact(7, C_3P2) ! C_3P2

  ! N3LO contacts
  call chp_set_chiral_D(0) ! Do not use

  ! Set needed ci and di, if any

  ! Set regulator parameter n
  call chp_set_1PE_reg_par(3.0D0)
  call chp_set_2PE_reg_par(3.0D0)
  call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
  call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
  call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
  call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
  call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
  call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
  call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
  call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
  call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2
  call chp_set_N3LO_contact_reg_par(1, 3.0D0) ! Dh_1S0
  call chp_set_N3LO_contact_reg_par(2, 3.0D0) ! D_1S0
  call chp_set_N3LO_contact_reg_par(3, 3.0D0) ! D_3P0
  call chp_set_N3LO_contact_reg_par(4, 3.0D0) ! D_1P1
  call chp_set_N3LO_contact_reg_par(5, 3.0D0) ! D_3P1
  call chp_set_N3LO_contact_reg_par(6, 3.0D0) ! Dh_3S1
  call chp_set_N3LO_contact_reg_par(7, 3.0D0) ! D_3S1
  call chp_set_N3LO_contact_reg_par(8, 3.0D0) ! D_3D1
  call chp_set_N3LO_contact_reg_par(9, 3.0D0) ! Dh_3S1-3D1
  call chp_set_N3LO_contact_reg_par(10, 3.0D0) ! D_3S1-3D1
  call chp_set_N3LO_contact_reg_par(11, 3.0D0) ! D_1D2
  call chp_set_N3LO_contact_reg_par(12, 3.0D0) ! D_3D2
  call chp_set_N3LO_contact_reg_par(13, 3.0D0) ! D_3P2
  call chp_set_N3LO_contact_reg_par(14, 3.0D0) ! D_3P2-3F2
  call chp_set_N3LO_contact_reg_par(15, 3.0D0) ! D_3D3

  ! Set pion exchange contributions
  ! Basic 1PE
  call chp_set_chiral_1PE(1) ! Use
  ! CIB effects in 1PE
  call chp_set_chiral_1PE_CIB(1) ! Use
  ! pion-gamma exchange
  call chp_set_chiral_1PE_gamma(0) ! Do not use
  ! Relativistic corrections to 1PE
  call chp_set_chiral_1PE_relcorr(0) ! Do not use

  ! Leading 2PE
  call chp_set_chiral_2PE_1loop_0(1) ! Use

  ! 1-loop 2PE proportional to ci
  call chp_set_chiral_2PE_1loop_d(0) ! Do not use

  ! 1-loop 2PE proportional to 1/M_N (relativistic corrections)
  call chp_set_chiral_2PE_1loop_r(0) ! Do not use

  ! 1-loop 2PE proportional to ci*cj)
  call chp_set_chiral_2PE_1loop_dd(0) ! Do not use

  ! 1-loop 2PE proportional to ci/M_N)
  call chp_set_chiral_2PE_1loop_dr(0) ! Do not use

  ! 1-loop 2PE proportional to 1/M_N^2)
  call chp_set_chiral_2PE_1loop_rr(0) ! Do not use

  ! 2-loop 2PE
  call chp_set_chiral_2PE_2loop(0) ! Do not use
  ! Contributions to 2-loop 2PE that do not have analytical expressions
  call chp_set_chiral_2PE_2loop_int(0) ! Do not use

  ! Use correct nucleon mass in 2PE relativistic corrections
  call chp_set_chiral_2PE_CSB_correct_mass(0) ! Do not use

  ! Use minimal relativity
  call chp_set_chiral_minimal_relativity(1) ! Use

  ! Use Kamada-Glockle transform
  call chp_set_chiral_kamada_glockle_transform(0) ! Do not use

  call chp_set_units_and_derive_constants

end subroutine chp_preset_ispot_nlo

! ############################################################################################################# !
! ############################################################################################################# !
! ############################################################################################################# !
! I. Svensson's next-to-next-to-leading-order potential (Phys. Rev. C 105, 014004 (2022))
subroutine chp_preset_ispot_n2lo(parameter_array)
  use idaho_chiral_potential

  implicit none

  real (8), dimension (0:15), intent (in) :: parameter_array
  
  ! LO contacts
  real (8) :: Ct_1S0pp
  real (8) :: Ct_3S1pp
  real (8) :: Ct_1S0np
  real (8) :: Ct_3S1np
  real (8) :: Ct_1S0nn
  real (8) :: Ct_3S1nn
  ! NLO contacts
  real (8) :: C_1S0
  real (8) :: C_3P0
  real (8) :: C_1P1
  real (8) :: C_3P1
  real (8) :: C_3S1
  real (8) :: C_3S1_3D1
  real (8) :: C_3P2
  ! N2LO ci
  real (8) :: c1
  real (8) :: c3
  real (8) :: c4

  !if (use_parameter_array) then
  Ct_1S0pp  = parameter_array(0)
  Ct_3S1pp  = parameter_array(1)
  Ct_1S0np  = parameter_array(2)
  Ct_3S1np  = parameter_array(3)
  Ct_1S0nn  = parameter_array(4)
  Ct_3S1nn  = parameter_array(5)
  C_1S0     = parameter_array(6)
  C_3P0     = parameter_array(7)
  C_1P1     = parameter_array(8)
  C_3P1     = parameter_array(9)
  C_3S1     = parameter_array(10)
  C_3S1_3D1 = parameter_array(11)
  C_3P2     = parameter_array(12)
  c1        = parameter_array(13)
  c3        = parameter_array(14)
  c4        = parameter_array(15)
  !else
  !  Ct_1S0pp  = -1.54250208697730012D-01
  !  Ct_3S1pp  = -1.48057121689049997D-01
  !  Ct_1S0np  = -1.55222568956019991D-01
  !  Ct_3S1np  = -1.48057121689049997D-01
  !  Ct_1S0nn  = -1.54958402671070000D-01
  !  Ct_3S1nn  = -1.48057121689049997D-01
  !  C_1S0     =  1.60603274119125006D+00
  !  C_3P0     =  1.19584298162829006D+00
  !  C_1P1     =  6.06776291764649978D-01
  !  C_3P1     = -2.49343337260960030D-01
  !  C_3S1     = -7.46664535867350043D-01
  !  C_3S1_3D1 =  1.56125611619209986D-01
  !  C_3P2     = -1.86669728671579993D-01
  !  c1        = -7.39616066183449994D-01
  !  c3        = -3.61062504021943020D+00
  !  c4        =  2.43675828100628999D+00
  !end if

  call initialize_chiral_potential

  ! GENERAL PARAMETERS AND CONSTANTS
  ! (proton, nucleon, neutron)
  call chp_set_mass_nucleon((/9.38272046000000046D+02, 9.38918267117927371D+02, 9.39565379000000007D+02/))
  ! (pi-, pi, pi+)
  call chp_set_mass_pion((/1.39570179999999993D+02, 1.34976599999999991D+02, 1.39570179999999993D+02/))

  call chp_set_chiral_order(NNLO)
  call chp_set_chiral_mode(chiral_mode_EM2015)
  call chp_set_reg("SF", 7.00000000000000000D+02)
  call chp_set_gA(1.28899999999999992D+00)
  call chp_set_fpi(9.22000000000000028D+01)
  call chp_set_fine_structure(7.29735256899999990D-03)

  call chp_set_Lambda(4.50000000000000000D+02)
  call chp_set_contact_format("PW")

  ! SET INCLUDED DIAGRAMS / TERMS

  ! LO contacts
  call chp_set_chiral_Ct_CIB(1) ! Use CIB contacts
  call chp_set_CIB_LO_contact(1, -1, Ct_1S0pp) ! Ct_1S0pp
  call chp_set_CIB_LO_contact(2, -1, Ct_3S1pp) ! Ct_3S1pp
  call chp_set_CIB_LO_contact(1,  0, Ct_1S0np) ! Ct_1S0np
  call chp_set_CIB_LO_contact(2,  0, Ct_3S1np) ! Ct_3S1np
  call chp_set_CIB_LO_contact(1,  1, Ct_1S0nn) ! Ct_1S0nn
  call chp_set_CIB_LO_contact(2,  1, Ct_3S1nn) ! Ct_3S1nn

  ! NLO contacts
  call chp_set_chiral_C(1) ! Use
  call chp_set_NLO_contact(1, C_1S0) ! C_1S0
  call chp_set_NLO_contact(2, C_3P0) ! C_3P0
  call chp_set_NLO_contact(3, C_1P1) ! C_1P1
  call chp_set_NLO_contact(4, C_3P1) ! C_3P1
  call chp_set_NLO_contact(5, C_3S1) ! C_3S1
  call chp_set_NLO_contact(6, C_3S1_3D1) ! C_3S1-3D1
  call chp_set_NLO_contact(7, C_3P2) ! C_3P2

  ! N3LO contacts
  call chp_set_chiral_D(0) ! Do not use

  ! Set needed ci and di, if any
  call chp_set_c1(c1)
  call chp_set_c3(c3)
  call chp_set_c4(c4)

  ! Set regulator parameter n
  call chp_set_1PE_reg_par(3.0D0)
  call chp_set_2PE_reg_par(3.0D0)
  call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
  call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
  call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
  call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
  call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
  call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
  call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
  call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
  call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2
  call chp_set_N3LO_contact_reg_par(1, 3.0D0) ! Dh_1S0
  call chp_set_N3LO_contact_reg_par(2, 3.0D0) ! D_1S0
  call chp_set_N3LO_contact_reg_par(3, 3.0D0) ! D_3P0
  call chp_set_N3LO_contact_reg_par(4, 3.0D0) ! D_1P1
  call chp_set_N3LO_contact_reg_par(5, 3.0D0) ! D_3P1
  call chp_set_N3LO_contact_reg_par(6, 3.0D0) ! Dh_3S1
  call chp_set_N3LO_contact_reg_par(7, 3.0D0) ! D_3S1
  call chp_set_N3LO_contact_reg_par(8, 3.0D0) ! D_3D1
  call chp_set_N3LO_contact_reg_par(9, 3.0D0) ! Dh_3S1-3D1
  call chp_set_N3LO_contact_reg_par(10, 3.0D0) ! D_3S1-3D1
  call chp_set_N3LO_contact_reg_par(11, 3.0D0) ! D_1D2
  call chp_set_N3LO_contact_reg_par(12, 3.0D0) ! D_3D2
  call chp_set_N3LO_contact_reg_par(13, 3.0D0) ! D_3P2
  call chp_set_N3LO_contact_reg_par(14, 3.0D0) ! D_3P2-3F2
  call chp_set_N3LO_contact_reg_par(15, 3.0D0) ! D_3D3

  ! Set pion exchange contributions
  ! Basic 1PE
  call chp_set_chiral_1PE(1) ! Use
  ! CIB effects in 1PE
  call chp_set_chiral_1PE_CIB(1) ! Use
  ! pion-gamma exchange
  call chp_set_chiral_1PE_gamma(0) ! Do not use
  ! Relativistic corrections to 1PE
  call chp_set_chiral_1PE_relcorr(0) ! Do not use

  ! Leading 2PE
  call chp_set_chiral_2PE_1loop_0(1) ! Use

  ! 1-loop 2PE proportional to ci
  call chp_set_chiral_2PE_1loop_d(1) ! Use

  ! 1-loop 2PE proportional to 1/M_N (relativistic corrections)
  call chp_set_chiral_2PE_1loop_r(0) ! Do not use

  ! 1-loop 2PE proportional to ci*cj)
  call chp_set_chiral_2PE_1loop_dd(0) ! Do not use

  ! 1-loop 2PE proportional to ci/M_N)
  call chp_set_chiral_2PE_1loop_dr(0) ! Do not use

  ! 1-loop 2PE proportional to 1/M_N^2)
  call chp_set_chiral_2PE_1loop_rr(0) ! Do not use

  ! 2-loop 2PE
  call chp_set_chiral_2PE_2loop(0) ! Do not use
  ! Contributions to 2-loop 2PE that do not have analytical expressions
  call chp_set_chiral_2PE_2loop_int(0) ! Do not use

  ! Use correct nucleon mass in 2PE relativistic corrections
  call chp_set_chiral_2PE_CSB_correct_mass(0) ! Do not use

  ! Use minimal relativity
  call chp_set_chiral_minimal_relativity(1) ! Use

  ! Use Kamada-Glockle transform
  call chp_set_chiral_kamada_glockle_transform(0) ! Do not use

  call chp_set_units_and_derive_constants

end subroutine chp_preset_ispot_n2lo
