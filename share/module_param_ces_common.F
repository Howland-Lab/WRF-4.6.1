MODULE module_param_ces_common

  USE ccpp_kind_types, ONLY : kind_phys

  IMPLICIT NONE

  SAVE

  ! MYNN only
  REAL(kind=kind_phys), PARAMETER :: mynn_pr = 0.74
  REAL(kind=kind_phys), PARAMETER :: mynn_b2 = 15.0
  REAL(kind=kind_phys), PARAMETER :: mynn_c2 = 0.729
  REAL(kind=kind_phys), PARAMETER :: mynn_c4 = 0.0

  ! Declare variables without PARAMETER attribute
  ! YSU 
  REAL(kind=kind_phys) :: ysu_afac, ysu_bfac, ysu_pfac, ysu_d1, ysu_d2
  ! MYNN
  REAL(kind=kind_phys) :: mynn_g1, mynn_b1, mynn_cns
  REAL(kind=kind_phys) :: mynn_c3, mynn_c5, mynn_beta
  REAL(kind=kind_phys) :: mynn_a1, mynn_c1, mynn_a2, mynn_g2
  REAL(kind=kind_phys) :: mynn_cc2, mynn_cc3, mynn_e1c, mynn_e2c
  REAL(kind=kind_phys) :: mynn_e3c, mynn_e4c, mynn_e5c
  ! SH
  REAL(kind=kind_phys) :: sh_afac, sh_bfac, sh_pfac, sh_d1, sh_d2
  ! REV MM5 - WSDF Charnock relation
  REAL(kind=kind_phys) :: wsdf_m, wsdf_b

CONTAINS

  SUBROUTINE init_module_param_ces_common

  IMPLICIT NONE

    INTEGER :: io_status
    CHARACTER*256 :: filename

    ! Open the file and handle errors if it does not exist
    filename = "wrf_ces_params"
    OPEN(FILE=TRIM(filename)//".txt", UNIT=15, FORM='FORMATTED', STATUS='OLD', IOSTAT=io_status)
    IF(io_status /= 0) THEN
      CALL wrf_error_fatal("Error opening file ./wrf_ces_params.txt")
    ENDIF
    ! Read parameters
! YSU
    READ(15, *) ysu_afac
    READ(15, *) ysu_bfac
    READ(15, *) ysu_pfac
    READ(15, *) ysu_d1
    READ(15, *) ysu_d2
! MYNN
    READ(15, *) mynn_g1
    READ(15, *) mynn_b1
    READ(15, *) mynn_c3
    READ(15, *) mynn_c5
    READ(15, *) mynn_cns
    READ(15, *) mynn_beta
! SH
    READ(15, *) sh_afac
    READ(15, *) sh_bfac
    READ(15, *) sh_pfac
    READ(15, *) sh_d1
    READ(15, *) sh_d2
! REV MM5 WSDF Charnock params
    READ(15, *) wsdf_m
    READ(15, *) wsdf_b
    ! Close the file
16  CLOSE(15)

    ! set global parameters for mynn
    CALL set_mynn_params()

  END SUBROUTINE init_module_param_ces_common

  SUBROUTINE set_mynn_params()

  IMPLICIT NONE

    mynn_a1  = mynn_b1*( 1.0-3.0*mynn_g1 )/6.0
!    c1  = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) )
    mynn_c1  = mynn_g1 -1.0/( 3.0*mynn_a1*2.88449914061481660)
    mynn_a2  = mynn_a1*( mynn_g1-mynn_c1 )/( mynn_g1*mynn_pr )
    mynn_g2  = mynn_b2/mynn_b1*( 1.0-mynn_c3 ) +2.0*mynn_a1/mynn_b1*( 3.0-2.0*mynn_c2 )
    mynn_cc2 =  1.0-mynn_c2
    mynn_cc3 =  1.0-mynn_c3
    mynn_e1c =  3.0*mynn_a2*mynn_b2*mynn_cc3
    mynn_e2c =  9.0*mynn_a1*mynn_a2*mynn_cc2
    mynn_e3c =  9.0*mynn_a2*mynn_a2*mynn_cc2*( 1.0-mynn_c5 )
    mynn_e4c = 12.0*mynn_a1*mynn_a2*mynn_cc2
    mynn_e5c =  6.0*mynn_a1*mynn_a1

  END SUBROUTINE set_mynn_params

END MODULE module_param_ces_common

