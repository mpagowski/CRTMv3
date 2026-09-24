!
! test_AAOD
!
! Test program for the CRTM aerosol Absorption Optical Depth (AAOD) K-matrix
! function CRTM_AAOD_K (CRTM_AAOD_Module, parallel to CRTM_AOD_Module).
! The K-matrix input is RTSolution_K%Layer_Absorption_Optical_Depth = 1 in every
! layer, so Atmosphere_K holds the Jacobian of the column AAOD.  The test_AOD
! profiles are used with an interior, off-node relative-humidity profile added
! (Load_Atm_Data.inc leaves RH = 0, the lower LUT boundary node).  Checks:
!   - CRTM_AAOD_K driven by Layer_Optical_Depth = 1 (AAOD input 0) reproduces
!     the untouched CRTM_AOD_K bitwise (Atmosphere_K and RTSolution)
!   - 0 <= dAAOD/dC <= dAOD/dC in every layer/aerosol, and 0 <= AAOD <= AOD
!   - Euler identity sum_{k,n} dAAOD/dC(k,n) * C(k,n) == column AAOD (linear in C)
!   - K-matrix vs central finite differences of CRTM_AAOD per profile and channel
!     along a NON-uniform concentration direction (a uniform relative direction
!     would only repeat the Euler identity) and along an RH direction
!   - Atmosphere_K against the saved reference results/k_matrix/test_AAOD_<sensor>.Atmosphere.bin
!     (created by the first run: a drift guard, not a validation)
!

PROGRAM test_AAOD

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  USE File_Utility, ONLY: File_Exists
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_AAOD'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/k_matrix/'
  CHARACTER(*), PARAMETER :: AEROSOLCOEFF_FILE = 'AerosolCoeff.GOCART-GEOS5.nc4'
  REAL(fp),     PARAMETER :: EPS_FD    = 1.0e-3_fp   ! relative concentration step (times f(k,n))
  REAL(fp),     PARAMETER :: H_RH      = 5.0e-6_fp   ! RH step (times d(k))
  REAL(fp),     PARAMETER :: FD_TOL    = 1.0e-6_fp   ! K vs finite difference (relative)
  REAL(fp),     PARAMETER :: RH_SIGNAL_MIN = 1.0e-6_fp   ! the RH check must be non-trivial somewhere
  REAL(fp),     PARAMETER :: REL_TOL   = 1.0e-12_fp  ! Euler identity (relative)
  REAL(fp),     PARAMETER :: FLOOR_REL = 1.0e-10_fp  ! zero-safe floor, relative to column OD
  REAL(fp),     PARAMETER :: FLOOR_RH  = 1.0e-8_fp   ! zero-safe floor of the per-layer RH check, relative to column AAOD:
                                                     ! the central FD of a layer carries a round-off of ~1e-16 x AAOD(k), so
                                                     ! layers with a zero RH derivative (dust; sea salt where w = 1) need
                                                     ! FD_TOL x floor >> 1e-16 x column AAOD; real signals are >= 1e-7 x column AAOD

  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS TEST ****
  !
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 3
  ! ...but only ONE Sensor at a time
  INTEGER, PARAMETER :: N_SENSORS = 1
  ! ============================================================================

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: k, l, m, n
  LOGICAL :: failed
  REAL(fp) :: col_aod, col_aaod, euler, fd, lin, denom, err_max_euler, err_max_fd
  REAL(fp) :: fd_rh, lin_rh, err_max_fd_rh, rh_signal, fac, err
  INTEGER  :: n_worst, worst(3)
  REAL(fp) :: worst_val(5)
  ! Declarations for Jacobian comparisons
  INTEGER :: n_l, n_m
  CHARACTER(256) :: atmk_File
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_k(:,:)


  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)

  ! Define the FORWARD variables
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm_p(N_PROFILES), Atm_m(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_OD(:,:), RTSolution_AOD(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_p(:,:), RTSolution_m(:,:)

  ! Define the K-MATRIX variables
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:), Atmosphere_K2(:,:), Atmosphere_KOD(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:), RTSolution_K2(:,:), RTSolution_KOD(:,:)
  ! ============================================================================

  !First, make sure the right number of inputs have been provided
  IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
     WRITE(*,*) TRIM(PROGRAM_NAME)//': ERROR, ONLY one command-line argument required, returning'
     STOP 1
  ENDIF
  CALL GET_COMMAND_ARGUMENT(1,Sensor_Id)   !read in the value

  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Test program for the CRTM aerosol Absorption Optical Depth (AAOD) K-matrix function.', &
    'CRTM Version: '//TRIM(Version) )


  ! Get sensor id from user
  ! -----------------------
  Sensor_Id = ADJUSTL(Sensor_Id)
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)



  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. Initialise for the requested sensor
  ! ---------------------------------------
  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  Error_Status = CRTM_Init( (/Sensor_Id/), &
                            ChannelInfo, &
                            Aerosol_Model = 'GOCART-GEOS5', &
                            AerosolCoeff_Format = 'netCDF', &
                            AerosolCoeff_File = AEROSOLCOEFF_FILE, &
                            File_Path=COEFFICIENTS_PATH)
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------
  n_Channels = CRTM_ChannelInfo_n_Channels(ChannelInfo(1))
  ! ============================================================================




  ! ============================================================================
  ! 3. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 3a. Allocate the ARRAYS
  ! -----------------------
  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), RTSolution_OD( n_Channels, N_PROFILES ), &
            RTSolution_AOD( n_Channels, N_PROFILES ), &
            RTSolution_p( n_Channels, N_PROFILES ), RTSolution_m( n_Channels, N_PROFILES ), &
            Atmosphere_K( n_Channels, N_PROFILES ), Atmosphere_K2( n_Channels, N_PROFILES ), &
            Atmosphere_KOD( n_Channels, N_PROFILES ), &
            RTSolution_K( n_Channels, N_PROFILES ), RTSolution_K2( n_Channels, N_PROFILES ), &
            RTSolution_KOD( n_Channels, N_PROFILES ), &
            STAT = Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 3b. Allocate the STRUCTURES
  ! ---------------------------
  ! The INPUT structures
  ! ...Forward variables
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere forward structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ...K-matrix variables
  CALL CRTM_RTSolution_Create( RTSolution_K  , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K2 , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_KOD, N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_K)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_K2)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_KOD)) ) THEN
    Message = 'Error allocating CRTM RTSolution K-matrix structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! The OUTPUT structures
  ! ...Forward variables
  CALL CRTM_RTSolution_Create( RTSolution   , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_OD, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_AOD, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_p , N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_m , N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_OD)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_AOD)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_p)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_m)) ) THEN
    Message = 'Error allocating CRTM RTSolution forward structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ...K-matrix variables
  CALL CRTM_Atmosphere_Create( Atmosphere_K  , N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_K2 , N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_KOD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_K)) .OR. &
       ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_K2)) .OR. &
       ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_KOD)) ) THEN
    Message = 'Error allocating CRTM Atmosphere K-matrix structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ============================================================================




  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! NOTE: Load_Atm_Data.inc is the test_AOD profile set and uses the CRTM-scheme
  !       aerosol type constants (DUST_AEROSOL=1, SEASALT_SSAM/SSCM1/SSCM2=2,3,4,
  !       SEASALT_SSCM3=5, SULFATE=8).  Against the GOCART-GEOS5 LUT loaded here
  !       those IDs select Dust 1-5 and Sea salt 3.  Relative_Humidity is not set
  !       by the profile (0 in every layer = the first LUT RH node); an interior,
  !       off-node RH profile is added here so that the RH Jacobians use interior
  !       interpolation stencils (nodes: 0.05 spacing up to 0.80, 0.01 above).
  CALL Load_Atm_Data()
  DO k = 1, N_LAYERS
    Atm(1)%Relative_Humidity(k) = Off_Node( 0.30_fp + 0.62_fp*REAL(k-1,fp)/REAL(N_LAYERS-1,fp) )
    Atm(2)%Relative_Humidity(k) = Off_Node( 0.92_fp - 0.62_fp*REAL(k-1,fp)/REAL(N_LAYERS-1,fp) )
  END DO
  ! ============================================================================




  ! ============================================================================
  ! 5. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
  !
  ! 5a. Zero the K-matrix OUTPUT structures
  ! ---------------------------------------
  CALL CRTM_Atmosphere_Zero( Atmosphere_K )
  CALL CRTM_Atmosphere_Zero( Atmosphere_K2 )
  CALL CRTM_Atmosphere_Zero( Atmosphere_KOD )

  ! 5b. Initialize the K-matrix INPUT
  ! ---------------------------------
  !     RTSolution_K  : AAOD input only  -> CRTM_AAOD_K  (the test subject)
  !     RTSolution_K2 : AOD  input only  -> CRTM_AAOD_K  (cross-check)
  !     RTSolution_KOD: AOD  input only  -> CRTM_AOD_K   (cross-check reference)
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      RTSolution_K(l,m)%Layer_Optical_Depth              = ZERO
      RTSolution_K(l,m)%Layer_Absorption_Optical_Depth   = ONE
      RTSolution_K2(l,m)%Layer_Optical_Depth             = ONE
      RTSolution_K2(l,m)%Layer_Absorption_Optical_Depth  = ZERO
      RTSolution_KOD(l,m)%Layer_Optical_Depth            = ONE
      RTSolution_KOD(l,m)%Layer_Absorption_Optical_Depth = ZERO
    END DO
  END DO
  ! ============================================================================




  ! ============================================================================
  ! 6. **** CALL THE CRTM AAOD K-MATRIX FUNCTION ****
  !
  Error_Status = CRTM_AAOD_K( Atm         , &
                              RTSolution_K, &
                              ChannelInfo , &
                              RTSolution  , &
                              Atmosphere_K  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM AAOD K-Matrix function'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 6b. Parallel-module cross-check: CRTM_AAOD_K driven by Layer_Optical_Depth
  !     only must reproduce the untouched CRTM_AOD_K bitwise
  ! ---------------------------------------------------------------------------
  Error_Status = CRTM_AAOD_K( Atm, RTSolution_K2 , ChannelInfo, RTSolution_OD, Atmosphere_K2  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM AAOD K-Matrix function (AOD input)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Error_Status = CRTM_AOD_K ( Atm, RTSolution_KOD, ChannelInfo, RTSolution_AOD, Atmosphere_KOD )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM AOD K-Matrix function (cross-check)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  failed = .FALSE.
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      ! ...K-matrix outputs of the two calls (RH, concentration, Reff) bitwise
      IF ( ANY(Atmosphere_K2(l,m)%Relative_Humidity /= Atmosphere_KOD(l,m)%Relative_Humidity) ) THEN
        WRITE( *,'(5x,"Relative_Humidity K mismatch: profile ",i0,", channel ",i0)' ) m, RTSolution(l,m)%Sensor_Channel
        failed = .TRUE.
      END IF
      DO n = 1, N_AEROSOLS
        IF ( ANY(Atmosphere_K2(l,m)%Aerosol(n)%Concentration    /= Atmosphere_KOD(l,m)%Aerosol(n)%Concentration   ) .OR. &
             ANY(Atmosphere_K2(l,m)%Aerosol(n)%Effective_Radius /= Atmosphere_KOD(l,m)%Aerosol(n)%Effective_Radius) ) THEN
          WRITE( *,'(5x,"Aerosol ",i0," K mismatch: profile ",i0,", channel ",i0)' ) n, m, RTSolution(l,m)%Sensor_Channel
          failed = .TRUE.
        END IF
      END DO
      ! ...forward outputs: CRTM_AOD_K leaves AAOD = 0, all three calls give the same OD,
      !    both CRTM_AAOD_K calls give the same AAOD
      IF ( ANY(RTSolution_AOD(l,m)%Layer_Optical_Depth /= RTSolution(l,m)%Layer_Optical_Depth) .OR. &
           ANY(RTSolution_OD(l,m)%Layer_Optical_Depth  /= RTSolution(l,m)%Layer_Optical_Depth) .OR. &
           ANY(RTSolution_OD(l,m)%Layer_Absorption_Optical_Depth /= RTSolution(l,m)%Layer_Absorption_Optical_Depth) .OR. &
           ANY(RTSolution_AOD(l,m)%Layer_Absorption_Optical_Depth /= ZERO) ) THEN
        WRITE( *,'(5x,"Forward OD/AAOD mismatch between the K calls: profile ",i0,", channel ",i0)' ) &
               m, RTSolution(l,m)%Sensor_Channel
        failed = .TRUE.
      END IF
      IF ( failed ) EXIT
    END DO
    IF ( failed ) EXIT
  END DO
  IF ( failed ) THEN
    Message = 'CRTM_AAOD_K with Layer_Optical_Depth input is not bitwise identical to CRTM_AOD_K!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Message = 'CRTM_AAOD_K (Layer_Optical_Depth input) bitwise identical to CRTM_AOD_K.'
  CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ! ============================================================================




  ! ============================================================================
  ! 7. **** OUTPUT AND CHECK THE RESULTS ****
  !
  ! 7a. Bounds: 0 <= AAOD <= AOD and 0 <= dAAOD/dC <= dAOD/dC
  ! -----------------------------------------------------------
  failed = .FALSE.
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      IF ( .NOT. ALL( RTSolution(l,m)%Layer_Absorption_Optical_Depth >= ZERO .AND. &
                      RTSolution(l,m)%Layer_Absorption_Optical_Depth <= RTSolution(l,m)%Layer_Optical_Depth ) ) failed = .TRUE.
      DO n = 1, N_AEROSOLS
        IF ( .NOT. ALL( Atmosphere_K(l,m)%Aerosol(n)%Concentration >= ZERO .AND. &
                        Atmosphere_K(l,m)%Aerosol(n)%Concentration <= Atmosphere_KOD(l,m)%Aerosol(n)%Concentration ) ) &
          failed = .TRUE.
      END DO
    END DO
  END DO
  IF ( failed ) THEN
    Message = 'AAOD or dAAOD/dC outside [0, AOD counterpart]!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Message = '0 <= AAOD <= AOD and 0 <= dAAOD/dC <= dAOD/dC for all profiles/channels/layers.'
  CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )

  ! 7b. Euler identity (AAOD linear homogeneous in C) and K-matrix vs central
  !     finite differences of CRTM_AAOD along
  !       (i)  a NON-uniform concentration direction dC(k,n) = EPS_FD.C(k,n).f(k,n),
  !            f = 1 + 0.5(-1)^k + 0.2 n  (a uniform relative direction would only
  !            repeat the Euler identity, since AAOD is homogeneous of degree 1 in C)
  !       (ii) an RH direction dRH(k) = H_RH.(1 + 0.5(-1)^k)
  ! -----------------------------------------------------------------------------
  Atm_p = Atm
  Atm_m = Atm
  DO m = 1, N_PROFILES
    DO n = 1, N_AEROSOLS
      DO k = 1, N_LAYERS
        fac = ONE + 0.5_fp*REAL((-1)**k,fp) + 0.2_fp*REAL(n,fp)
        Atm_p(m)%Aerosol(n)%Concentration(k) = Atm(m)%Aerosol(n)%Concentration(k) * (ONE + EPS_FD*fac)
        Atm_m(m)%Aerosol(n)%Concentration(k) = Atm(m)%Aerosol(n)%Concentration(k) * (ONE - EPS_FD*fac)
      END DO
    END DO
  END DO
  Error_Status = CRTM_AAOD( Atm_p, ChannelInfo, RTSolution_p )
  IF ( Error_Status /= SUCCESS ) STOP 1
  Error_Status = CRTM_AAOD( Atm_m, ChannelInfo, RTSolution_m )
  IF ( Error_Status /= SUCCESS ) STOP 1
  failed = .FALSE.
  err_max_euler = ZERO
  err_max_fd    = ZERO
  DO m = 1, N_PROFILES
    WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)
    WRITE( *,'(/5x,"Channel",7x,"Column AOD",14x,"Column AAOD",9x,"sum K_C.C (Euler)",6x,"K.dC/eps",16x,"FD/eps")')
    DO l = 1, n_Channels
      col_aod  = SUM( RTSolution(l,m)%Layer_Optical_Depth(1:N_LAYERS) )
      col_aaod = SUM( RTSolution(l,m)%Layer_Absorption_Optical_Depth(1:N_LAYERS) )
      euler = ZERO
      lin   = ZERO
      DO n = 1, N_AEROSOLS
        DO k = 1, N_LAYERS
          fac = ONE + 0.5_fp*REAL((-1)**k,fp) + 0.2_fp*REAL(n,fp)
          euler = euler + Atmosphere_K(l,m)%Aerosol(n)%Concentration(k) * Atm(m)%Aerosol(n)%Concentration(k)
          lin   = lin   + Atmosphere_K(l,m)%Aerosol(n)%Concentration(k) * Atm(m)%Aerosol(n)%Concentration(k) * fac
        END DO
      END DO
      lin = EPS_FD * lin
      fd  = ( SUM(RTSolution_p(l,m)%Layer_Absorption_Optical_Depth(1:N_LAYERS)) - &
              SUM(RTSolution_m(l,m)%Layer_Absorption_Optical_Depth(1:N_LAYERS)) ) / TWO
      WRITE( *,'(5x,i7,5(2x,es22.15))') RTSolution(l,m)%Sensor_Channel, col_aod, col_aaod, euler, lin/EPS_FD, fd/EPS_FD
      denom = MAX( col_aaod, FLOOR_REL*col_aod )
      err_max_euler = MAX( err_max_euler, ABS(euler-col_aaod)/denom )
      IF ( .NOT. ( ABS(euler-col_aaod) <= REL_TOL*denom ) ) failed = .TRUE.
      denom = MAX( ABS(lin), ABS(fd), FLOOR_REL*EPS_FD*col_aod )
      err_max_fd = MAX( err_max_fd, ABS(lin-fd)/denom )
      IF ( .NOT. ( ABS(lin-fd) <= FD_TOL*denom ) ) failed = .TRUE.
    END DO
  END DO
  WRITE( *,'(/5x,"Max rel. error: Euler identity ",es10.3,", K vs finite difference (concentration direction) ",es10.3)') &
         err_max_euler, err_max_fd
  IF ( failed ) THEN
    Message = 'K-matrix inconsistent with the forward model (Euler identity or concentration finite difference)!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ...RH direction
  Atm_p = Atm
  Atm_m = Atm
  DO m = 1, N_PROFILES
    DO k = 1, N_LAYERS
      fac = ONE + 0.5_fp*REAL((-1)**k,fp)
      Atm_p(m)%Relative_Humidity(k) = Atm(m)%Relative_Humidity(k) + H_RH*fac
      Atm_m(m)%Relative_Humidity(k) = Atm(m)%Relative_Humidity(k) - H_RH*fac
    END DO
  END DO
  Error_Status = CRTM_AAOD( Atm_p, ChannelInfo, RTSolution_p )
  IF ( Error_Status /= SUCCESS ) STOP 1
  Error_Status = CRTM_AAOD( Atm_m, ChannelInfo, RTSolution_m )
  IF ( Error_Status /= SUCCESS ) STOP 1
  failed = .FALSE.
  err_max_fd_rh = ZERO
  rh_signal     = ZERO
  n_worst = 0
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      col_aod  = SUM( RTSolution(l,m)%Layer_Optical_Depth(1:N_LAYERS) )
      col_aaod = SUM( RTSolution(l,m)%Layer_Absorption_Optical_Depth(1:N_LAYERS) )
      ! ...per layer: the aerosol optics of layer k depend on RH(k) only, so the
      !    per-layer difference of the jointly perturbed profiles is the exact
      !    central FD of that layer (and avoids the round-off of a column sum)
      DO k = 1, N_LAYERS
        fac = ONE + 0.5_fp*REAL((-1)**k,fp)
        lin_rh = H_RH * Atmosphere_K(l,m)%Relative_Humidity(k) * fac
        fd_rh  = ( RTSolution_p(l,m)%Layer_Absorption_Optical_Depth(k) - &
                   RTSolution_m(l,m)%Layer_Absorption_Optical_Depth(k) ) / TWO
        denom = MAX( ABS(lin_rh), ABS(fd_rh), FLOOR_RH*MAX(col_aaod, FLOOR_REL*col_aod) )
        err = ABS(lin_rh-fd_rh)/denom
        rh_signal = MAX( rh_signal, ABS(lin_rh)/(H_RH*MAX(col_aaod,FLOOR_REL*col_aod)) )
        IF ( err > err_max_fd_rh ) THEN
          err_max_fd_rh = err
          worst = (/ m, RTSolution(l,m)%Sensor_Channel, k /)
          worst_val = (/ Atm(m)%Relative_Humidity(k), Atmosphere_K(l,m)%Relative_Humidity(k), lin_rh, fd_rh, &
                         RTSolution(l,m)%Layer_Absorption_Optical_Depth(k) /)
        END IF
        IF ( .NOT. ( ABS(lin_rh-fd_rh) <= FD_TOL*denom ) ) THEN
          failed = .TRUE.
          IF ( n_worst < 12 ) THEN
            n_worst = n_worst + 1
            WRITE( *,'(5x,"RH FD mismatch: profile ",i0," channel ",i0," layer ",i0," RH=",f7.4," K_RH=",es13.5, &
                   &" K.dRH=",es13.5," FD=",es13.5," rel=",es9.2," AAOD(k)=",es12.4," OD(k)=",es12.4, &
                   &" C(1:3)=",3(es10.3,1x))') m, RTSolution(l,m)%Sensor_Channel, k, Atm(m)%Relative_Humidity(k), &
                   Atmosphere_K(l,m)%Relative_Humidity(k), lin_rh, fd_rh, err, &
                   RTSolution(l,m)%Layer_Absorption_Optical_Depth(k), RTSolution(l,m)%Layer_Optical_Depth(k), &
                   (Atm(m)%Aerosol(n)%Concentration(k), n = 1, N_AEROSOLS)
          END IF
        END IF
      END DO
    END DO
  END DO
  WRITE( *,'(5x,"K vs finite difference (RH direction, per layer): max rel. error ",es10.3," at profile ",i0, &
             &" channel ",i0," layer ",i0," (RH=",f7.4," K_RH=",es13.5," K.dRH=",es13.5," FD=",es13.5, &
             &" AAOD(k)=",es12.4,"); max |K_RH.d| / column AAOD ",es10.3)') &
         err_max_fd_rh, worst, worst_val, rh_signal
  IF ( failed .OR. .NOT. ( rh_signal > RH_SIGNAL_MIN ) ) THEN
    Message = 'K-matrix inconsistent with the forward model in the RH direction (or RH sensitivity absent everywhere)!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Message = 'K-matrix consistent with CRTM_AAOD: Euler identity to 1e-12, finite differences to 1e-6 (concentration and RH).'
  CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ! ============================================================================

  ! ============================================================================
  ! 9. **** COMPARE Atmosphere_K RESULTS TO SAVED VALUES ****
  !
  WRITE( *, '( /5x, "Comparing calculated results with saved ones..." )' )

  ! 9a. Create the output file if it does not exist
  ! -----------------------------------------------
  ! ...Generate filename
  atmk_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.Atmosphere.bin'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(atmk_File) ) THEN
    Message = 'Atmosphere_K save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write Atmosphere_K structure to file
    Error_Status = CRTM_Atmosphere_WriteFile( atmk_file, Atmosphere_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating Atmosphere_K save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

  ! 9b. Inquire the saved file
  ! --------------------------
  Error_Status = CRTM_Atmosphere_InquireFile( atmk_File, &
                                              n_Channels = n_l, &
                                              n_Profiles = n_m )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring Atmosphere_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9c. Compare the dimensions
  ! --------------------------
  IF ( n_l /= n_Channels .OR. n_m /= N_PROFILES ) THEN
    Message = 'Dimensions of saved data different from that calculated!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9d. Read the saved data
  ! -----------------------
  Error_Status = CRTM_Atmosphere_ReadFile( atmk_File, atm_k, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading Atmosphere_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9e. Compare the Jacobians with the saved reference (drift guard: the file
  !     was created by the first run of this test, it is not an independent value)
  ! ----------------------------------------------------------------------------
  IF ( ALL(CRTM_Atmosphere_Compare(Atmosphere_K, atm_k)) ) THEN
    Message = 'Atmosphere_K Jacobians are the same as the saved reference (drift guard).'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'Atmosphere_K Jacobians are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    ! Write the current Atmosphere_K results to file
    atmk_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.Atmosphere.bin'
    Error_Status = CRTM_Atmosphere_WriteFile( atmk_file, Atmosphere_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating temporary Atmosphere_K save file for failed comparison'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
    STOP 1
  END IF
  ! ============================================================================

  ! ============================================================================
  ! 8. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ============================================================================

  ! ============================================================================
  ! 10. **** CLEAN UP ****
  !
  ! 10a. Deallocate the structures
  ! ------------------------------
  CALL CRTM_Atmosphere_Destroy(atm_K)
  CALL CRTM_Atmosphere_Destroy(Atmosphere_K)
  CALL CRTM_Atmosphere_Destroy(Atmosphere_K2)
  CALL CRTM_Atmosphere_Destroy(Atmosphere_KOD)
  CALL CRTM_RTSolution_Destroy(RTSolution_K)
  CALL CRTM_RTSolution_Destroy(RTSolution_K2)
  CALL CRTM_RTSolution_Destroy(RTSolution_KOD)
  CALL CRTM_RTSolution_Destroy(RTSolution)
  CALL CRTM_RTSolution_Destroy(RTSolution_OD)
  CALL CRTM_RTSolution_Destroy(RTSolution_AOD)
  CALL CRTM_RTSolution_Destroy(RTSolution_p)
  CALL CRTM_RTSolution_Destroy(RTSolution_m)
  CALL CRTM_Atmosphere_Destroy(Atm)
  CALL CRTM_Atmosphere_Destroy(Atm_p)
  CALL CRTM_Atmosphere_Destroy(Atm_m)

  ! 10b. Deallocate the arrays
  ! --------------------------
  DEALLOCATE(RTSolution, RTSolution_OD, RTSolution_AOD, RTSolution_p, RTSolution_m, &
             RTSolution_K, RTSolution_K2, RTSolution_KOD, &
             Atmosphere_K, Atmosphere_K2, Atmosphere_KOD, atm_k, &
             STAT = Allocate_Status)
  ! ============================================================================

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'

  ! Nudge a relative humidity off the GOCART-GEOS5 LUT RH nodes (0.05 spacing up
  ! to 0.80, 0.01 above) by at least 0.002, so that finite differences of +-1e-5
  ! stay inside one interpolation stencil
  FUNCTION Off_Node( rh_in ) RESULT( rh )
    REAL(fp), INTENT(IN) :: rh_in
    REAL(fp) :: rh, spacing, node
    rh = rh_in
    IF ( rh <= 0.8_fp ) THEN
      spacing = 0.05_fp
    ELSE
      spacing = 0.01_fp
    END IF
    node = REAL( NINT( rh / spacing ), fp ) * spacing
    IF ( ABS( rh - node ) < 0.002_fp ) rh = node + 0.0025_fp
  END FUNCTION Off_Node

END PROGRAM test_AAOD
