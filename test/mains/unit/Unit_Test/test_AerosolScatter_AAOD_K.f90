
!
! test_AerosolScatter_AAOD_K
!
! AAOD version of test_AerosolScatter_K: unit test for the aerosol absorption
! optical depth K-matrix, CRTM_AAOD_K (CRTM_AAOD_Module), driven by
! RTSolution_K%Layer_Absorption_Optical_Depth = 1 in every layer so that
! Atmosphere_K holds the Jacobian of the column AAOD.
! (The original test drives the radiance K-matrix, CRTM_K_Matrix, with
! Brightness_Temperature_K = 1; there is no radiance path for AAOD.)
! Same sensor (modis_aqua, all channels; channel 1 inspected), aerosol
! (GOCART-GEOS5 sea salt 1, hygroscopic) and 10-layer profile as the original,
! except that RH and concentration differ from layer to layer (RH 0.53 .. 0.755
! in steps of 0.025, all off the LUT RH nodes at 0.05 spacing; concentration
! 1e-5 x (1 + 0.1 (k-1))) so that the per-layer finite-difference check is
! not ten copies of one check and a layer-index error would show.
! Verifies that the aerosol RH Jacobians are non-zero, that the concentration
! and RH Jacobians of channel 1 match central finite differences of CRTM_AAOD
! layer by layer, and that CRTM_AAOD_K driven by Layer_Optical_Depth = 1
! reproduces the untouched CRTM_AOD_K bitwise.
! Not covered here (single aerosol, scattering on, no clamp active): the
! multi-aerosol accumulation, the scattering-off branch and the QC clamps —
! see the k_matrix regression test test_AAOD (6 sensors, 2 profiles, 3 aerosols).
!

PROGRAM test_AerosolScatter_AAOD_K

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,        ONLY: ZERO, ONE, TWO
  USE CRTM_SpcCoeff,          ONLY: SC, CRTM_SpcCoeff_Load, CRTM_SpcCoeff_Destroy
  USE CRTM_AerosolCoeff,      ONLY: AeroC, CRTM_AerosolCoeff_Load, CRTM_AerosolCoeff_Destroy
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type, &
                                    CRTM_Atmosphere_Create, &
                                    CRTM_Atmosphere_Destroy, &
                                    CRTM_Atmosphere_Zero, &
                                    H2O_ID, O3_ID, CO2_ID, N2O_ID, CH4_ID, CO_ID, &
                                    MASS_MIXING_RATIO_UNITS
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type
  USE CRTM_RTSolution_Define,  ONLY: CRTM_RTSolution_type, &
                                     CRTM_RTSolution_Create, &
                                     CRTM_RTSolution_Destroy, &
                                     CRTM_RTSolution_Zero
  USE CRTM_LifeCycle,          ONLY: CRTM_Init, CRTM_Destroy
  USE CRTM_AOD_Module,         ONLY: CRTM_AOD_K
  USE CRTM_AAOD_Module,        ONLY: CRTM_AAOD, CRTM_AAOD_K

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_AAOD_K'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  INTEGER :: Error_Status
  CHARACTER(256) :: Sensor_Id = 'modis_aqua'

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Atmosphere_type)  :: Atm(1), Atm_p(1), Atm_m(1)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:), Atm_K2(:,:), Atm_KOD(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_K(:,:), RTSolution_K2(:,:), RTSolution_KOD(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_p(:,:), RTSolution_m(:,:)

  INTEGER :: n_Channels
  INTEGER :: n_Layers = 10
  INTEGER :: n_Absorbers = 6
  INTEGER :: n_Aerosols = 1
  INTEGER :: n_Clouds = 0
  REAL(fp), PARAMETER :: RH_0      = 0.53_fp    ! layer-1 RH; layer k: RH_0 + 0.025 (k-1), all 0.005 off the LUT nodes
  REAL(fp), PARAMETER :: DRH_LAYER = 0.025_fp
  REAL(fp), PARAMETER :: EPS_C     = 1.0e-3_fp  ! relative concentration step
  REAL(fp), PARAMETER :: H_RH      = 2.0e-6_fp  ! RH step; per-layer differences keep the round-off at ~1e-11
  REAL(fp), PARAMETER :: FD_TOL_C  = 1.0e-10_fp ! AAOD is linear in the concentration (observed ~1e-13)
  REAL(fp), PARAMETER :: FD_TOL_RH = 1.0e-8_fp  ! h^2 truncation ~2e-11 + round-off ~3e-11 at h = 2e-6 (observed ~1e-10)
  REAL(fp), PARAMETER :: FLOOR_REL = 1.0e-10_fp ! division guard of the relative comparison only (never active here)

  INTEGER :: l, k
  REAL(fp) :: jac, fd, denom, floor_val
  LOGICAL :: failed

  PRINT *, 'Starting test_AerosolScatter_AAOD_K...'

  ! Initialize CRTM
  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Load GOCART-GEOS5
  Error_Status = CRTM_AerosolCoeff_Load('GOCART-GEOS5', 'AerosolCoeff.GOCART-GEOS5.bin', File_Path=COEFFICIENTS_PATH)
  IF ( Error_Status /= SUCCESS ) STOP 1

  n_Channels = ChannelInfo(1)%n_Channels
  ALLOCATE( RTSolution( n_Channels, 1 ), RTSolution_K( n_Channels, 1 ), &
            RTSolution_K2( n_Channels, 1 ), RTSolution_KOD( n_Channels, 1 ), &
            RTSolution_p( n_Channels, 1 ), RTSolution_m( n_Channels, 1 ), &
            Atm_K( n_Channels, 1 ), Atm_K2( n_Channels, 1 ), Atm_KOD( n_Channels, 1 ) )

  ! Create structures
  CALL CRTM_Atmosphere_Create( Atm, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_K, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_K2, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_KOD, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_RTSolution_Create( RTSolution, n_Layers )
  CALL CRTM_RTSolution_Create( RTSolution_K, n_Layers )
  CALL CRTM_RTSolution_Create( RTSolution_K2, n_Layers )
  CALL CRTM_RTSolution_Create( RTSolution_KOD, n_Layers )
  CALL CRTM_RTSolution_Create( RTSolution_p, n_Layers )
  CALL CRTM_RTSolution_Create( RTSolution_m, n_Layers )

  ! Populate Input (as in test_AerosolScatter_K)
  Atm(1)%Absorber_ID(1) = H2O_ID
  Atm(1)%Absorber_ID(2) = O3_ID
  Atm(1)%Absorber_ID(3) = CO2_ID
  Atm(1)%Absorber_ID(4) = N2O_ID
  Atm(1)%Absorber_ID(5) = CH4_ID
  Atm(1)%Absorber_ID(6) = CO_ID
  Atm(1)%Aerosol(1)%Type = 6 ! Sea Salt (hygroscopic)

  ! Standard Atmosphere (simplified)
  Atm(1)%Pressure(1:n_Layers) = (/ (1000.0_fp - REAL(k-1,fp)*100.0_fp, k=1,n_Layers) /)
  Atm(1)%Temperature(1:n_Layers) = 280.0_fp
  Atm(1)%Absorber_Units = MASS_MIXING_RATIO_UNITS
  Atm(1)%Absorber(1:n_Layers, 1) = 1.0e-3_fp ! H2O
  Atm(1)%Absorber(1:n_Layers, 2) = 1.0e-6_fp ! O3
  Atm(1)%Absorber(1:n_Layers, 3:6) = 1.0e-7_fp
  Atm(1)%Aerosol(1)%Concentration(1:n_Layers) = (/ (1.0e-5_fp*(ONE + 0.1_fp*REAL(k-1,fp)), k=1,n_Layers) /)
  Atm(1)%Aerosol(1)%Effective_Radius(1:n_Layers) = 1.0_fp
  Atm(1)%Relative_Humidity(1:n_Layers) = (/ (RH_0 + DRH_LAYER*REAL(k-1,fp), k=1,n_Layers) /)

  ! Initialize K-Matrix inputs: sensitivity to the column AAOD
  CALL CRTM_Atmosphere_Zero( Atm_K )
  CALL CRTM_RTSolution_Zero( RTSolution_K )
  DO l = 1, n_Channels
    RTSolution_K(l,1)%Layer_Absorption_Optical_Depth = ONE
  END DO

  ! Call K-Matrix
  Error_Status = CRTM_AAOD_K( Atm, RTSolution_K, ChannelInfo, RTSolution, Atm_K )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Inspect Aerosol Jacobians for Channel 1
  PRINT *, 'Column AAOD Jacobians for Channel 1 (', TRIM(Sensor_Id), ', wavenumber ', &
           SC(ChannelInfo(1)%Sensor_Index)%Wavenumber(ChannelInfo(1)%Channel_Index(1)), ' cm-1):'
  PRINT *, 'Layer, RH, RH_Jacobian, Concentration_Jacobian, AAOD(k), OD(k)'
  DO k = 1, n_Layers
    WRITE(*, '(i3, 1x,f6.3, 4(1x,es12.5))') k, Atm(1)%Relative_Humidity(k), Atm_K(1,1)%Relative_Humidity(k), &
                                   Atm_K(1,1)%Aerosol(1)%Concentration(k), &
                                   RTSolution(1,1)%Layer_Absorption_Optical_Depth(k), RTSolution(1,1)%Layer_Optical_Depth(k)
  END DO

  ! Check if RH Jacobian is non-zero
  IF ( ANY(ABS(Atm_K(1,1)%Relative_Humidity) > 1.0e-12_fp) ) THEN
    PRINT *, 'SUCCESS: AAOD RH Jacobians are NON-ZERO.'
  ELSE
    PRINT *, 'FAIL: AAOD RH Jacobians are ZERO.'
    STOP 1
  END IF

  ! --------------------------------------------------------------------------
  ! Jacobians of channel 1 vs central finite differences of CRTM_AAOD, per layer
  ! (perturbing layer k changes only layer k of the AAOD profile)
  ! --------------------------------------------------------------------------
  failed = .FALSE.
  floor_val = FLOOR_REL * SUM( RTSolution(1,1)%Layer_Optical_Depth(1:n_Layers) )
  DO k = 1, n_Layers
    ! ...concentration
    Atm_p(1) = Atm(1); Atm_m(1) = Atm(1)
    Atm_p(1)%Aerosol(1)%Concentration(k) = Atm(1)%Aerosol(1)%Concentration(k) * ( ONE + EPS_C )
    Atm_m(1)%Aerosol(1)%Concentration(k) = Atm(1)%Aerosol(1)%Concentration(k) * ( ONE - EPS_C )
    Error_Status = CRTM_AAOD( Atm_p, ChannelInfo, RTSolution_p ); IF ( Error_Status /= SUCCESS ) STOP 1
    Error_Status = CRTM_AAOD( Atm_m, ChannelInfo, RTSolution_m ); IF ( Error_Status /= SUCCESS ) STOP 1
    fd  = ( RTSolution_p(1,1)%Layer_Absorption_Optical_Depth(k) - RTSolution_m(1,1)%Layer_Absorption_Optical_Depth(k) ) / &
          ( TWO * EPS_C * Atm(1)%Aerosol(1)%Concentration(k) )
    jac = Atm_K(1,1)%Aerosol(1)%Concentration(k)
    denom = MAX( ABS(jac), ABS(fd), floor_val )
    WRITE(*,'(" dAAOD/dC  layer ",i2,": K=",es22.15," FD=",es22.15," rel=",es9.2)') k, jac, fd, ABS(jac-fd)/denom
    IF ( .NOT. ( ABS(jac-fd) <= FD_TOL_C*denom ) ) failed = .TRUE.
    ! ...relative humidity
    Atm_p(1) = Atm(1); Atm_m(1) = Atm(1)
    Atm_p(1)%Relative_Humidity(k) = Atm(1)%Relative_Humidity(k) + H_RH
    Atm_m(1)%Relative_Humidity(k) = Atm(1)%Relative_Humidity(k) - H_RH
    Error_Status = CRTM_AAOD( Atm_p, ChannelInfo, RTSolution_p ); IF ( Error_Status /= SUCCESS ) STOP 1
    Error_Status = CRTM_AAOD( Atm_m, ChannelInfo, RTSolution_m ); IF ( Error_Status /= SUCCESS ) STOP 1
    fd  = ( RTSolution_p(1,1)%Layer_Absorption_Optical_Depth(k) - RTSolution_m(1,1)%Layer_Absorption_Optical_Depth(k) ) / &
          ( TWO * H_RH )
    jac = Atm_K(1,1)%Relative_Humidity(k)
    denom = MAX( ABS(jac), ABS(fd), floor_val )
    WRITE(*,'(" dAAOD/dRH layer ",i2,": K=",es22.15," FD=",es22.15," rel=",es9.2)') k, jac, fd, ABS(jac-fd)/denom
    IF ( .NOT. ( ABS(jac-fd) <= FD_TOL_RH*denom ) ) failed = .TRUE.
  END DO
  IF ( failed ) THEN
    PRINT *, 'FAIL: AAOD K-matrix vs finite differences.'
    STOP 1
  END IF
  PRINT *, 'SUCCESS: AAOD K-matrix (concentration, RH) matches finite differences.'

  ! --------------------------------------------------------------------------
  ! CRTM_AAOD_K driven by Layer_Optical_Depth = 1 (AAOD input 0) must
  ! reproduce the untouched CRTM_AOD_K bitwise
  ! --------------------------------------------------------------------------
  CALL CRTM_RTSolution_Zero( RTSolution_K2 )
  CALL CRTM_RTSolution_Zero( RTSolution_KOD )
  DO l = 1, n_Channels
    RTSolution_K2(l,1)%Layer_Optical_Depth  = ONE
    RTSolution_KOD(l,1)%Layer_Optical_Depth = ONE
  END DO
  CALL CRTM_Atmosphere_Zero( Atm_K2 )
  CALL CRTM_Atmosphere_Zero( Atm_KOD )
  Error_Status = CRTM_AAOD_K( Atm, RTSolution_K2 , ChannelInfo, RTSolution_p, Atm_K2  ); IF ( Error_Status /= SUCCESS ) STOP 1
  Error_Status = CRTM_AOD_K ( Atm, RTSolution_KOD, ChannelInfo, RTSolution_m, Atm_KOD ); IF ( Error_Status /= SUCCESS ) STOP 1
  failed = .FALSE.
  DO l = 1, n_Channels
    IF ( ANY(Atm_K2(l,1)%Relative_Humidity        /= Atm_KOD(l,1)%Relative_Humidity       ) .OR. &
         ANY(Atm_K2(l,1)%Aerosol(1)%Concentration /= Atm_KOD(l,1)%Aerosol(1)%Concentration) .OR. &
         ANY(RTSolution_p(l,1)%Layer_Optical_Depth /= RTSolution_m(l,1)%Layer_Optical_Depth) ) failed = .TRUE.
  END DO
  IF ( failed ) THEN
    PRINT *, 'FAIL: CRTM_AAOD_K with Layer_Optical_Depth input is not bitwise identical to CRTM_AOD_K.'
    STOP 1
  END IF
  PRINT *, 'SUCCESS: CRTM_AAOD_K (Layer_Optical_Depth input) == CRTM_AOD_K bitwise.'

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_p )
  CALL CRTM_Atmosphere_Destroy( Atm_m )
  CALL CRTM_Atmosphere_Destroy( Atm_K )
  CALL CRTM_Atmosphere_Destroy( Atm_K2 )
  CALL CRTM_Atmosphere_Destroy( Atm_KOD )
  CALL CRTM_RTSolution_Destroy( RTSolution )
  CALL CRTM_RTSolution_Destroy( RTSolution_K )
  CALL CRTM_RTSolution_Destroy( RTSolution_K2 )
  CALL CRTM_RTSolution_Destroy( RTSolution_KOD )
  CALL CRTM_RTSolution_Destroy( RTSolution_p )
  CALL CRTM_RTSolution_Destroy( RTSolution_m )
  DEALLOCATE( RTSolution, RTSolution_K, RTSolution_K2, RTSolution_KOD, RTSolution_p, RTSolution_m, Atm_K, Atm_K2, Atm_KOD )
  Error_Status = CRTM_Destroy( ChannelInfo )
  Error_Status = CRTM_AerosolCoeff_Destroy()

  PRINT *, 'SUCCESS: test_AerosolScatter_AAOD_K passed.'

END PROGRAM test_AerosolScatter_AAOD_K
