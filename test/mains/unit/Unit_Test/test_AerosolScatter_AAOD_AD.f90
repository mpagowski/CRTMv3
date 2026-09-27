!
! test_AerosolScatter_AAOD_AD
!
! Unit test for the absorption aerosol optical depth (AAOD) adjoint
! computed by CRTM_AOD_AD with AAOD = .TRUE. (GOCART-GEOS5 table).
! The profile contains one aerosol of every type in the table.
!
! Verifies adjoint consistency with the tangent linear model,
!   <TL(dx), y> = <dx, AD(y)>,  y = TL(dx)
! for a perturbation dx of the aerosol concentration and relative humidity.
!
! CREATION HISTORY:
!       Written by:     Mariusz Pagowski, 2026/9/26
!                       mariusz.pagowski@colorado.edu
!

PROGRAM test_AerosolScatter_AAOD_AD

  USE CRTM_Module
  USE CRTM_AerosolCoeff, ONLY: AeroC

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_AAOD_AD'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  INTEGER,  PARAMETER :: N_LAYERS = 5
  REAL(fp), PARAMETER :: RH(N_LAYERS) = (/ 0.43_fp, 0.57_fp, 0.66_fp, 0.72_fp, 0.34_fp /)
  REAL(fp), PARAMETER :: TOLERANCE = 1.0e-12_fp

  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  INTEGER :: Error_Status, n_Channels, n_Aerosols, l, k, n
  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Options_type)     :: Options(1)
  TYPE(CRTM_Atmosphere_type)  :: Atm(1), Atm_TL(1), Atm_AD(1)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:), RTSolution_AD(:,:)
  REAL(fp) :: sum_tl, sum_ad, rel_err

  PRINT *, 'Starting test_AerosolScatter_AAOD_AD...'

  ! Initialize CRTM with the GOCART-GEOS5 aerosol table
  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, &
                            Aerosol_Model     = 'GOCART-GEOS5', &
                            AerosolCoeff_File = 'AerosolCoeff.GOCART-GEOS5.bin', &
                            File_Path         = COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) STOP 1
  n_Channels = ChannelInfo(1)%n_Channels
  n_Aerosols = AeroC%n_Types

  ALLOCATE( RTSolution(n_Channels,1), RTSolution_TL(n_Channels,1), RTSolution_AD(n_Channels,1) )
  CALL CRTM_RTSolution_Create( RTSolution   , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_AD, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm   , N_LAYERS, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, 1, 0, n_Aerosols )
  ! The AOD calculation only uses the aerosol and RH profiles
  Options(1)%Check_Input = .FALSE.

  ! Populate Input: every aerosol type in every layer
  Atm(1)%Absorber_ID(1) = H2O_ID
  Atm(1)%Relative_Humidity = RH
  DO n = 1, n_Aerosols
    Atm(1)%Aerosol(n)%Type = AeroC%Type(n)
    Atm(1)%Aerosol(n)%Effective_Radius = AeroC%Reff(1,n)
    Atm(1)%Aerosol(n)%Concentration = 1.0e-4_fp*(ONE + 0.1_fp*REAL(n,fp)) * &
                                      (/ 1.0_fp, 2.0_fp, 0.5_fp, 1.5_fp, 0.8_fp /)
  END DO

  ! Populate TL Input: perturb concentration and relative humidity
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  DO n = 1, n_Aerosols
    DO k = 1, N_LAYERS
      Atm_TL(1)%Aerosol(n)%Concentration(k) = 0.1_fp * Atm(1)%Aerosol(n)%Concentration(k) * &
                                              SIN(0.731_fp*REAL(k,fp) + 1.37_fp*REAL(n,fp))
    END DO
  END DO
  Atm_TL(1)%Relative_Humidity = (/ (SIN(0.731_fp*REAL(k,fp) + 0.5_fp), k=1,N_LAYERS) /)

  ! TL Call
  Error_Status = CRTM_AOD_TL( Atm, Atm_TL, ChannelInfo, RTSolution, RTSolution_TL, &
                              Options=Options, AAOD=.TRUE. )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Adjoint Input: y = AAOD_TL
  sum_tl = ZERO
  DO l = 1, n_Channels
    RTSolution_AD(l,1)%Layer_Optical_Depth = RTSolution_TL(l,1)%Layer_Optical_Depth
    sum_tl = sum_tl + SUM(RTSolution_TL(l,1)%Layer_Optical_Depth**2)
  END DO

  ! AD Call
  CALL CRTM_Atmosphere_Zero( Atm_AD )
  Error_Status = CRTM_AOD_AD( Atm, RTSolution_AD, ChannelInfo, RTSolution, Atm_AD, &
                              Options=Options, AAOD=.TRUE. )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Adjoint Test: <TL(dx), y> == <dx, AD(y)>
  sum_ad = SUM(Atm_TL(1)%Relative_Humidity * Atm_AD(1)%Relative_Humidity)
  DO n = 1, n_Aerosols
    sum_ad = sum_ad + &
             SUM(Atm_TL(1)%Aerosol(n)%Concentration      * Atm_AD(1)%Aerosol(n)%Concentration) + &
             SUM(Atm_TL(1)%Aerosol(n)%Effective_Radius   * Atm_AD(1)%Aerosol(n)%Effective_Radius) + &
             SUM(Atm_TL(1)%Aerosol(n)%Effective_Variance * Atm_AD(1)%Aerosol(n)%Effective_Variance)
  END DO

  PRINT *, 'LHS (TL sum): ', sum_tl
  PRINT *, 'RHS (AD sum): ', sum_ad
  PRINT *, 'Diff:         ', ABS(sum_tl - sum_ad)

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_AD )
  CALL CRTM_RTSolution_Destroy( RTSolution )
  CALL CRTM_RTSolution_Destroy( RTSolution_TL )
  CALL CRTM_RTSolution_Destroy( RTSolution_AD )
  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_AD )
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( sum_tl > ZERO ) THEN
    rel_err = ABS(sum_tl - sum_ad) / sum_tl
    PRINT *, 'Rel. Diff:    ', rel_err
    IF ( rel_err < TOLERANCE ) THEN
       PRINT *, 'SUCCESS: AAOD adjoint test passed.'
    ELSE
       PRINT *, 'FAIL: AAOD adjoint test failed.'
       STOP 1
    END IF
  ELSE
    PRINT *, 'FAIL: LHS is ZERO, test is invalid.'
    STOP 1
  END IF

END PROGRAM test_AerosolScatter_AAOD_AD
