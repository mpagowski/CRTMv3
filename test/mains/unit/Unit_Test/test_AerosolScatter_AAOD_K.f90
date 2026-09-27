!
! test_AerosolScatter_AAOD_K
!
! Unit test for the absorption aerosol optical depth (AAOD) Jacobians
! computed by CRTM_AOD_K with AAOD = .TRUE. (GOCART-GEOS5 table).
! The profile contains one aerosol of every type in the table.
!
! The K-matrix input is the column AAOD, RTSolution_K%Layer_Optical_Depth = 1,
! so Atm_K holds d(column AAOD)/dx for every channel. Every concentration
! and relative humidity Jacobian element is checked against the column
! AAOD TL (CRTM_AOD_TL with AAOD = .TRUE.) for a unit perturbation.
!
! CREATION HISTORY:
!       Written by:     Mariusz Pagowski, 2026/9/26
!                       mariusz.pagowski@colorado.edu
!

PROGRAM test_AerosolScatter_AAOD_K

  USE CRTM_Module
  USE CRTM_AerosolCoeff, ONLY: AeroC

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_AAOD_K'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  INTEGER,  PARAMETER :: N_LAYERS = 5
  REAL(fp), PARAMETER :: RH(N_LAYERS) = (/ 0.43_fp, 0.57_fp, 0.66_fp, 0.72_fp, 0.34_fp /)
  REAL(fp), PARAMETER :: TOLERANCE = 1.0e-12_fp

  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  INTEGER :: Error_Status, n_Channels, n_Aerosols, l, k, n
  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Options_type)     :: Options(1)
  TYPE(CRTM_Atmosphere_type)  :: Atm(1), Atm_TL(1)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:), RTSolution_K(:,:)
  REAL(fp) :: col_tl, jac, max_err, max_jac

  PRINT *, 'Starting test_AerosolScatter_AAOD_K...'

  ! Initialize CRTM with the GOCART-GEOS5 aerosol table
  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, &
                            Aerosol_Model     = 'GOCART-GEOS5', &
                            AerosolCoeff_File = 'AerosolCoeff.GOCART-GEOS5.bin', &
                            File_Path         = COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) STOP 1
  n_Channels = ChannelInfo(1)%n_Channels
  n_Aerosols = AeroC%n_Types

  ALLOCATE( RTSolution(n_Channels,1), RTSolution_TL(n_Channels,1), RTSolution_K(n_Channels,1), &
            Atm_K(n_Channels,1) )
  CALL CRTM_RTSolution_Create( RTSolution   , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K , N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm   , N_LAYERS, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_K , N_LAYERS, 1, 0, n_Aerosols )
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

  ! K-matrix input: column AAOD for every channel
  DO l = 1, n_Channels
    RTSolution_K(l,1)%Layer_Optical_Depth = ONE
  END DO

  ! Call K-Matrix
  CALL CRTM_Atmosphere_Zero( Atm_K )
  Error_Status = CRTM_AOD_K( Atm, RTSolution_K, ChannelInfo, RTSolution, Atm_K, &
                             Options=Options, AAOD=.TRUE. )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Inspect Aerosol Jacobians for Channel 1
  PRINT *, 'Column AAOD Jacobians for Channel 1, aerosol type 1:'
  PRINT *, 'Layer, RH_Jacobian, Concentration_Jacobian'
  DO k = 1, N_LAYERS
    WRITE(*, '(i3, 2(1x,es12.5))') k, Atm_K(1,1)%Relative_Humidity(k), Atm_K(1,1)%Aerosol(1)%Concentration(k)
  END DO

  ! Check every Jacobian element against the TL for a unit perturbation
  max_err = ZERO
  max_jac = ZERO
  DO k = 1, N_LAYERS
    DO n = 0, n_Aerosols        ! n = 0: relative humidity
      CALL CRTM_Atmosphere_Zero( Atm_TL )
      IF ( n == 0 ) THEN
        Atm_TL(1)%Relative_Humidity(k) = ONE
      ELSE
        Atm_TL(1)%Aerosol(n)%Concentration(k) = ONE
      END IF
      Error_Status = CRTM_AOD_TL( Atm, Atm_TL, ChannelInfo, RTSolution, RTSolution_TL, &
                                  Options=Options, AAOD=.TRUE. )
      IF ( Error_Status /= SUCCESS ) STOP 1
      DO l = 1, n_Channels
        col_tl = SUM(RTSolution_TL(l,1)%Layer_Optical_Depth)
        IF ( n == 0 ) THEN
          jac = Atm_K(l,1)%Relative_Humidity(k)
        ELSE
          jac = Atm_K(l,1)%Aerosol(n)%Concentration(k)
        END IF
        max_err = MAX( max_err, ABS(jac - col_tl) )
        max_jac = MAX( max_jac, ABS(col_tl) )
      END DO
    END DO
  END DO
  PRINT *, 'Max |K - TL|: ', max_err
  PRINT *, 'Max |TL|:     ', max_jac

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_K )
  CALL CRTM_RTSolution_Destroy( RTSolution )
  CALL CRTM_RTSolution_Destroy( RTSolution_TL )
  CALL CRTM_RTSolution_Destroy( RTSolution_K )
  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_K, Atm_K )
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( max_jac <= ZERO ) THEN
    PRINT *, 'FAIL: AAOD Jacobians are ZERO, test is invalid.'
    STOP 1
  ELSE IF ( max_err > TOLERANCE*max_jac ) THEN
    PRINT *, 'FAIL: AAOD Jacobians differ from the TL.'
    STOP 1
  END IF
  PRINT *, 'SUCCESS: AAOD Jacobians agree with the TL.'

END PROGRAM test_AerosolScatter_AAOD_K
