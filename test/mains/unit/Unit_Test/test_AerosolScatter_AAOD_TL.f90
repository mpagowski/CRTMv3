!
! test_AerosolScatter_AAOD_TL
!
! Unit test for the absorption aerosol optical depth (AAOD) computed by
! CRTM_AOD and CRTM_AOD_TL with AAOD = .TRUE. (GOCART-GEOS5 table).
! The profile contains one aerosol of every type in the table.
!
! Verifies that
!   - the forward AAOD equals SUM_n (1-w).ke.rho, where the extinction
!     coefficient ke and single scatter albedo w are taken from
!     CRTM_Compute_AerosolScatter WITH scattering
!   - the TL AAOD agrees with a central finite difference of the forward
!     AAOD for perturbations of concentration, relative humidity, and both.
!
! CREATION HISTORY:
!       Written by:     Mariusz Pagowski, 2026/9/26
!                       mariusz.pagowski@colorado.edu
!

PROGRAM test_AerosolScatter_AAOD_TL

  USE CRTM_Module
  USE CRTM_Parameters,       ONLY: MAX_N_LEGENDRE_TERMS, MAX_N_PHASE_ELEMENTS
  USE CRTM_AerosolCoeff,     ONLY: AeroC
  USE CRTM_AtmOptics_Define, ONLY: CRTM_AtmOptics_type, CRTM_AtmOptics_Create, &
                                   CRTM_AtmOptics_Destroy, CRTM_AtmOptics_Zero
  USE CRTM_AerosolScatter,   ONLY: CRTM_Compute_AerosolScatter
  USE ASvar_Define,          ONLY: ASvar_type, ASvar_Create, ASvar_Destroy

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_AAOD_TL'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  INTEGER,  PARAMETER :: N_LAYERS = 5
  ! Relative humidities away from the 0.05-spaced LUT RH nodes
  REAL(fp), PARAMETER :: RH(N_LAYERS) = (/ 0.43_fp, 0.57_fp, 0.66_fp, 0.72_fp, 0.34_fp /)
  ! AAOD vs SUM (1-w).ke.rho: the LUT interpolation can give a slightly
  ! negative w (e.g. ~ -8e-5 for organic carbon near 1210 cm-1). The reference
  ! clips w to 0, but the AAOD uses ke*(1-w) before clipping.
  REAL(fp), PARAMETER :: FWD_RTOL = 1.0e-4_fp
  REAL(fp), PARAMETER :: FD_STEP  = 1.0e-4_fp
  REAL(fp), PARAMETER :: FD_RTOL  = 1.0e-5_fp

  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  INTEGER :: Error_Status, n_Channels, n_Aerosols, l, k, n, id, n_Failures
  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Options_type)     :: Options(1)
  TYPE(CRTM_Atmosphere_type)  :: Atm(1), Atm_TL(1), Atm_p(1), Atm_m(1)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:), RTS_p(:,:), RTS_m(:,:)
  REAL(fp), ALLOCATABLE :: aaod(:,:), ref(:,:), aaod_tl(:,:), aaod_fd(:,:)
  REAL(fp) :: rel_err
  CHARACTER(20) :: dname(3) = (/ 'concentration       ', 'relative humidity   ', 'concentration + RH  ' /)

  PRINT *, 'Starting test_AerosolScatter_AAOD_TL...'
  n_Failures = 0

  ! Initialize CRTM with the GOCART-GEOS5 aerosol table
  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, &
                            Aerosol_Model     = 'GOCART-GEOS5', &
                            AerosolCoeff_File = 'AerosolCoeff.GOCART-GEOS5.bin', &
                            File_Path         = COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) STOP 1
  n_Channels = ChannelInfo(1)%n_Channels
  n_Aerosols = AeroC%n_Types

  ALLOCATE( RTSolution(n_Channels,1), RTSolution_TL(n_Channels,1), RTS_p(n_Channels,1), RTS_m(n_Channels,1), &
            aaod(N_LAYERS,n_Channels), ref(N_LAYERS,n_Channels), &
            aaod_tl(N_LAYERS,n_Channels), aaod_fd(N_LAYERS,n_Channels) )
  CALL CRTM_RTSolution_Create( RTSolution   , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_p        , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_m        , N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, 1, 0, n_Aerosols )
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


  ! ------------------------------------------
  ! 1. Forward: AAOD = SUM_n (1-w).ke.rho
  ! ------------------------------------------
  Error_Status = CRTM_AOD( Atm, ChannelInfo, RTSolution, Options=Options, AAOD=.TRUE. )
  IF ( Error_Status /= SUCCESS ) STOP 1
  DO l = 1, n_Channels
    aaod(:,l) = RTSolution(l,1)%Layer_Optical_Depth
  END DO
  CALL Reference_AAOD( ref )
  PRINT '(1x,"Channel 1 column AAOD: ",es16.8,", reference: ",es16.8)', SUM(aaod(:,1)), SUM(ref(:,1))
  CALL Check( 'AAOD is non-zero', MAXVAL(ref) > ZERO )
  rel_err = MAXVAL(ABS(aaod - ref)) / MAXVAL(ABS(ref))
  PRINT '(1x,"Rel. Diff: ",es12.4)', rel_err
  CALL Check( 'Forward AAOD = SUM (1-w).ke.rho', rel_err <= FWD_RTOL )


  ! ------------------------------------------
  ! 2. TL vs central finite difference
  ! ------------------------------------------
  DO id = 1, 3
    CALL Set_Perturbation( id )
    Error_Status = CRTM_AOD_TL( Atm, Atm_TL, ChannelInfo, RTSolution, RTSolution_TL, &
                                Options=Options, AAOD=.TRUE. )
    IF ( Error_Status /= SUCCESS ) STOP 1
    DO l = 1, n_Channels
      aaod_tl(:,l) = RTSolution_TL(l,1)%Layer_Optical_Depth
    END DO

    Atm_p = Atm
    Atm_m = Atm
    Atm_p(1)%Relative_Humidity = Atm(1)%Relative_Humidity + FD_STEP*Atm_TL(1)%Relative_Humidity
    Atm_m(1)%Relative_Humidity = Atm(1)%Relative_Humidity - FD_STEP*Atm_TL(1)%Relative_Humidity
    DO n = 1, n_Aerosols
      Atm_p(1)%Aerosol(n)%Concentration = Atm(1)%Aerosol(n)%Concentration + &
                                          FD_STEP*Atm_TL(1)%Aerosol(n)%Concentration
      Atm_m(1)%Aerosol(n)%Concentration = Atm(1)%Aerosol(n)%Concentration - &
                                          FD_STEP*Atm_TL(1)%Aerosol(n)%Concentration
    END DO
    Error_Status = CRTM_AOD( Atm_p, ChannelInfo, RTS_p, Options=Options, AAOD=.TRUE. )
    IF ( Error_Status /= SUCCESS ) STOP 1
    Error_Status = CRTM_AOD( Atm_m, ChannelInfo, RTS_m, Options=Options, AAOD=.TRUE. )
    IF ( Error_Status /= SUCCESS ) STOP 1
    DO l = 1, n_Channels
      aaod_fd(:,l) = (RTS_p(l,1)%Layer_Optical_Depth - RTS_m(l,1)%Layer_Optical_Depth) / (2.0_fp*FD_STEP)
    END DO

    PRINT '(1x,"TL vs FD, ",a)', TRIM(dname(id))
    PRINT '(3x,"max |AAOD TL| = ",es12.4)', MAXVAL(ABS(aaod_tl))
    IF ( MAXVAL(ABS(aaod_tl)) > ZERO ) THEN
      rel_err = MAXVAL(ABS(aaod_tl - aaod_fd)) / MAXVAL(ABS(aaod_tl))
      PRINT '(3x,"Rel. Diff: ",es12.4)', rel_err
      CALL Check( 'TL vs FD, '//TRIM(dname(id)), rel_err <= FD_RTOL )
    ELSE
      CALL Check( 'TL vs FD (AAOD TL is ZERO), '//TRIM(dname(id)), .FALSE. )
    END IF
  END DO


  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_p )
  CALL CRTM_Atmosphere_Destroy( Atm_m )
  CALL CRTM_RTSolution_Destroy( RTSolution )
  CALL CRTM_RTSolution_Destroy( RTSolution_TL )
  CALL CRTM_RTSolution_Destroy( RTS_p )
  CALL CRTM_RTSolution_Destroy( RTS_m )
  DEALLOCATE( RTSolution, RTSolution_TL, RTS_p, RTS_m, aaod, ref, aaod_tl, aaod_fd )
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( n_Failures > 0 ) THEN
    PRINT '(1x,"FAIL: ",i0," AAOD check(s) failed.")', n_Failures
    STOP 1
  END IF
  PRINT *, 'SUCCESS: All AAOD forward and TL checks passed.'

CONTAINS

  ! Reference AAOD = SUM_n (1-w).ke.rho, with ke (extinction) and w from
  ! the aerosol scattering calculation with scattering included
  SUBROUTINE Reference_AAOD( ref )
    REAL(fp), INTENT(OUT) :: ref(:,:)
    TYPE(CRTM_AtmOptics_type) :: AScat
    TYPE(ASvar_type)          :: ASvar
    INTEGER :: l, k, n
    CALL CRTM_AtmOptics_Create( AScat, N_LAYERS, MAX_N_LEGENDRE_TERMS, AeroC%n_Phase_Elements )
    AScat%n_Legendre_Terms   = 4
    AScat%Include_Scattering = .TRUE.
    CALL ASvar_Create( ASvar, MAX_N_LEGENDRE_TERMS, MAX_N_PHASE_ELEMENTS, N_LAYERS, n_Aerosols )
    DO l = 1, n_Channels
      CALL CRTM_AtmOptics_Zero( AScat )
      Error_Status = CRTM_Compute_AerosolScatter( Atm(1), ChannelInfo(1)%Sensor_Index, &
                                                  ChannelInfo(1)%Channel_Index(l), AScat, ASvar )
      IF ( Error_Status /= SUCCESS ) STOP 1
      DO k = 1, N_LAYERS
        ref(k,l) = ZERO
        DO n = 1, n_Aerosols
          ref(k,l) = ref(k,l) + (ONE - ASvar%w(k,n)) * ASvar%ke(k,n) * Atm(1)%Aerosol(n)%Concentration(k)
        END DO
      END DO
    END DO
    CALL CRTM_AtmOptics_Destroy( AScat )
    CALL ASvar_Destroy( ASvar )
  END SUBROUTINE Reference_AAOD


  ! Perturbation directions: 1 concentration, 2 RH, 3 both
  SUBROUTINE Set_Perturbation( id )
    INTEGER, INTENT(IN) :: id
    INTEGER :: n, k
    Atm_TL = Atm
    CALL CRTM_Atmosphere_Zero( Atm_TL )
    IF ( id /= 2 ) THEN
      DO n = 1, n_Aerosols
        DO k = 1, N_LAYERS
          Atm_TL(1)%Aerosol(n)%Concentration(k) = 0.1_fp * Atm(1)%Aerosol(n)%Concentration(k) * &
                                                  SIN(0.731_fp*REAL(k,fp) + 1.37_fp*REAL(n,fp))
        END DO
      END DO
    END IF
    IF ( id /= 1 ) Atm_TL(1)%Relative_Humidity = (/ (SIN(0.731_fp*REAL(k,fp) + 0.5_fp), k=1,N_LAYERS) /)
  END SUBROUTINE Set_Perturbation


  SUBROUTINE Check( label, passed )
    CHARACTER(*), INTENT(IN) :: label
    LOGICAL     , INTENT(IN) :: passed
    IF ( passed ) THEN
      PRINT '(1x,"PASS: ",a)', label
    ELSE
      PRINT '(1x,"FAIL: ",a)', label
      n_Failures = n_Failures + 1
    END IF
  END SUBROUTINE Check

END PROGRAM test_AerosolScatter_AAOD_TL
