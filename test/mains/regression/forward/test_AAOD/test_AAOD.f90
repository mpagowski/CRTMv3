!
! test_AAOD
!
! Test program for the CRTM aerosol Absorption Optical Depth (AAOD) Forward
! function CRTM_AAOD (CRTM_AAOD_Module, parallel to CRTM_AOD_Module), which
! returns both RTSolution%Layer_Optical_Depth and
! RTSolution%Layer_Absorption_Optical_Depth.  The AOD part is cross-checked
! bitwise against the untouched CRTM_AOD.
!
!

PROGRAM test_AAOD

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  USE File_Utility, ONLY: File_Exists
  USE Compare_Float_Numbers, ONLY: Compares_Within_Tolerance, DEFAULT_N_SIGFIG
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_AAOD'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/forward/'
  CHARACTER(*), PARAMETER :: AEROSOLCOEFF_FILE = 'AerosolCoeff.GOCART-GEOS5.nc4'

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
  INTEGER :: l, m
  LOGICAL :: failed
  ! Declarations for RTSolution comparison
  INTEGER :: n_l, n_m, n_k
  CHARACTER(256) :: rts_File
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)


  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  ! ...Output of the untouched CRTM_AOD, for the parallel-module cross-check
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_AOD(:,:)
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
    'Test program for the CRTM Aerosol Optical Depth (AOD) Forward function '//&
    'with the absorption optical depth (AAOD) output.', &
    'CRTM Version: '//TRIM(Version) )


  ! Get sensor id from user
  ! -----------------------
  Sensor_Id = ADJUSTL(Sensor_Id)
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)



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
  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), &
            RTSolution_AOD( n_Channels, N_PROFILES ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 3b. Allocate the STRUCTURES
  ! ---------------------------
  ! The input structure
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! The output structure
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_AOD, N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_AOD)) ) THEN
    Message = 'Error allocating CRTM RTSolution structures'
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
  !       those IDs select Dust 1-5 and Sea salt 3, so the 'sea salt' profile 2 is
  !       in fact dust bins 2-4 and its AAOD/AOD (~0.1) is a dust value.
  CALL Load_Atm_Data()
  ! ============================================================================




  ! ============================================================================
  ! 5. **** CALL THE CRTM AAOD MODEL ****
  !
  Error_Status = CRTM_AAOD( Atm        , &
                            ChannelInfo, &
                            RTSolution   )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM AAOD Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 5b. Parallel-module cross-check: the untouched CRTM_AOD on the same inputs
  !     must give bitwise identical Layer_Optical_Depth and Backscat_Coefficient
  ! ---------------------------------------------------------------------------
  Error_Status = CRTM_AOD( Atm           , &
                           ChannelInfo   , &
                           RTSolution_AOD  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM AOD Model (cross-check)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  failed = .FALSE.
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      IF ( ANY(RTSolution(l,m)%Layer_Optical_Depth  /= RTSolution_AOD(l,m)%Layer_Optical_Depth ) .OR. &
           ANY(RTSolution(l,m)%Backscat_Coefficient /= RTSolution_AOD(l,m)%Backscat_Coefficient) .OR. &
           ANY(RTSolution_AOD(l,m)%Layer_Absorption_Optical_Depth /= ZERO) ) THEN
        WRITE( *,'(5x,"CRTM_AAOD vs CRTM_AOD mismatch: profile ",i0,", channel ",i0)' ) m, RTSolution(l,m)%Sensor_Channel
        failed = .TRUE.
      END IF
    END DO
  END DO
  IF ( failed ) THEN
    Message = 'CRTM_AAOD Layer_Optical_Depth is not bitwise identical to CRTM_AOD!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Message = 'CRTM_AAOD Layer_Optical_Depth/Backscat_Coefficient bitwise identical to CRTM_AOD (which leaves AAOD = 0).'
  CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ! ============================================================================




  ! ============================================================================
  ! 6. **** OUTPUT THE RESULTS TO SCREEN ****
  !
  ! 6a. Per-channel column AOD and AAOD
  ! -----------------------------------
  DO m = 1, N_PROFILES
    WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(Sensor_Id)
    WRITE( *,'(/5x,"Channel",7x,"Column AOD",14x,"Column AAOD")')
    DO l = 1, n_Channels
      WRITE( *,'(5x,i7,2(2x,es22.15))') RTSolution(l,m)%Sensor_Channel, &
                                        SUM(RTSolution(l,m)%Layer_Optical_Depth), &
                                        SUM(RTSolution(l,m)%Layer_Absorption_Optical_Depth)
    END DO
  END DO

  ! 6b. Check the physical bounds 0 <= AAOD <= AOD in every layer
  ! -------------------------------------------------------------
  WRITE( *, '( /5x, "Checking 0 <= AAOD <= AOD for every layer..." )' )
  failed = .FALSE.
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      IF ( ANY(RTSolution(l,m)%Layer_Absorption_Optical_Depth < ZERO) .OR. &
           ANY(RTSolution(l,m)%Layer_Absorption_Optical_Depth > RTSolution(l,m)%Layer_Optical_Depth) ) THEN
        WRITE( *,'(5x,"AAOD bounds violated: profile ",i0,", channel ",i0)' ) m, RTSolution(l,m)%Sensor_Channel
        failed = .TRUE.
      END IF
    END DO
  END DO
  IF ( failed ) THEN
    Message = 'Layer AAOD outside [0,AOD]!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Message = 'Layer AAOD within [0,AOD] for all profiles/channels/layers.'
  CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ! ============================================================================

  ! ============================================================================
  ! 8. **** COMPARE RTSolution RESULTS TO SAVED VALUES ****
  !
  WRITE( *, '( /5x, "Comparing calculated results with saved ones..." )' )

  ! 8a. Create the output file if it does not exist
  ! -----------------------------------------------
  ! ...Generate a filename
  rts_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution.nc'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(rts_File) ) THEN
    Message = 'RTSolution save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write RTSolution structure to file
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution, NetCDF=.TRUE., Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating RTSolution save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

  ! 8b. Inquire the saved file
  ! --------------------------
  Error_Status = CRTM_RTSolution_InquireFile( rts_File, &
                                              NetCDF=.TRUE.,    &
                                              n_Channels = n_l, &
                                              n_Layers   = n_k, &
                                              n_Profiles = n_m )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring RTSolution save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 8c. Compare the dimensions
  ! --------------------------
  IF ( n_l /= n_Channels .OR. n_m /= N_PROFILES .OR. n_k /= N_LAYERS ) THEN
    Message = 'Dimensions of saved data different from that calculated!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 8d. Allocate the structure to read in saved data
  ! ------------------------------------------------
  ALLOCATE( rts( n_l, n_m ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating RTSolution saved data array'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  CALL CRTM_RTSolution_Create( rts, n_k )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(rts)) ) THEN
    Message = 'Error allocating CRTM RTSolution saved data structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 8e. Read the saved data
  ! -----------------------
  Error_Status = CRTM_RTSolution_ReadFile( rts_File, rts, NetCDF=.TRUE., Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading RTSolution save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 8f. Compare the per-layer AOD and AAOD explicitly
  !     (CRTM_RTSolution_Compare does not include Layer_Absorption_Optical_Depth)
  ! ---------------------------------------------------------------------------
  failed = .FALSE.
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      IF ( .NOT. ALL(Compares_Within_Tolerance( RTSolution(l,m)%Layer_Optical_Depth(1:N_LAYERS), &
                                                rts(l,m)%Layer_Optical_Depth(1:N_LAYERS), &
                                                DEFAULT_N_SIGFIG )) ) THEN
        WRITE( *,'(5x,"Layer_Optical_Depth differs: profile ",i0,", channel ",i0)' ) &
               m, RTSolution(l,m)%Sensor_Channel
        failed = .TRUE.
      END IF
      IF ( .NOT. ALL(Compares_Within_Tolerance( RTSolution(l,m)%Layer_Absorption_Optical_Depth(1:N_LAYERS), &
                                                rts(l,m)%Layer_Absorption_Optical_Depth(1:N_LAYERS), &
                                                DEFAULT_N_SIGFIG )) ) THEN
        WRITE( *,'(5x,"Layer_Absorption_Optical_Depth differs: profile ",i0,", channel ",i0)' ) &
               m, RTSolution(l,m)%Sensor_Channel
        failed = .TRUE.
      END IF
    END DO
  END DO
  IF ( .NOT. failed ) THEN
    Message = 'RTSolution AOD and AAOD results are the same!'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'RTSolution AOD and/or AAOD results are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    ! Write the current RTSolution results to file
    rts_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution.nc'
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution, NetCDF=.TRUE., Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating temporary RTSolution save file for failed comparison'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
    STOP 1
  END IF
  ! ============================================================================

  ! ============================================================================
  ! 7. **** DESTROY THE CRTM ****
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
  ! 9. **** CLEAN UP ****
  !
  ! 9a. Deallocate the structures
  ! -----------------------------
  CALL CRTM_RTSolution_Destroy(RTSolution)
  CALL CRTM_RTSolution_Destroy(RTSolution_AOD)
  CALL CRTM_RTSolution_Destroy(rts)
  CALL CRTM_Atmosphere_Destroy(Atm)

  ! 9b. Deallocate the arrays
  ! -------------------------
  DEALLOCATE(RTSolution, RTSolution_AOD, rts, STAT=Allocate_Status)
  ! ============================================================================

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_AAOD
