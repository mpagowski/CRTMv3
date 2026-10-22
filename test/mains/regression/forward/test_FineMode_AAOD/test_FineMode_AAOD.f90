!
! test_FineMode_AAOD
!
! Test program for the CRTM fine-mode Absorption Aerosol Optical Depth
! (AAOD), computed by the AOD Forward function with FINE_MODE = .TRUE. and
! AAOD = .TRUE.
! The concentration of each aerosol type is scaled by its fmass1(RH), the
! fraction of its mass in particles with diameter < 1 um. The result is
! returned in RTSolution%Layer_Optical_Depth.
!
! Also checks that 0 <= fine-mode AAOD <= total AAOD in every layer, and
! that the column fine-mode AAOD < total AAOD for every channel and profile.
!
!
! CREATION HISTORY:
!       Written by:     Mariusz Pagowski, 2026/10/22
!                       mariusz.pagowski@colorado.edu
!                       Adapted from test_AAOD for the fine-mode AAOD
!

PROGRAM test_FineMode_AAOD

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
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_FineMode_AAOD'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/forward/'

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
  ! Relative tolerance for fine-mode <= total
  REAL(fp), PARAMETER :: FM_RTOL = 1.0e-12_fp
  ! ============================================================================



  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  CHARACTER(256) :: AerosolCoeff_File
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: l, m
  LOGICAL :: has_new_coeff
  LOGICAL :: has_old_coeff
  ! Declarations for RTSolution comparison
  INTEGER :: n_l, n_m
  CHARACTER(256) :: rts_File
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)


  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_Total(:,:)
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
    'Test program for the CRTM fine-mode AAOD Forward function.', &
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
  has_new_coeff = File_Exists(COEFFICIENTS_PATH//'AerosolCoeff.GOCART-GEOS5.BRC.kb.v2.nc')
  has_old_coeff = File_Exists(COEFFICIENTS_PATH//'AerosolCoeff.GOCART-GEOS5.nc4')
  IF ( has_new_coeff ) THEN
    AerosolCoeff_File = 'AerosolCoeff.GOCART-GEOS5.BRC.kb.v2.nc'
  ELSE
    AerosolCoeff_File = 'AerosolCoeff.GOCART-GEOS5.nc4'
  END IF
  Error_Status = CRTM_Init( (/Sensor_Id/), &
                            ChannelInfo, &
                            Aerosol_Model = 'GOCART-GEOS5', &
                            AerosolCoeff_Format = 'netCDF', &
                            AerosolCoeff_File = TRIM(AerosolCoeff_File), &
                            File_Path=COEFFICIENTS_PATH)
  IF ( Error_Status /= SUCCESS .AND. has_old_coeff .AND. TRIM(AerosolCoeff_File) /= 'AerosolCoeff.GOCART-GEOS5.nc4' ) THEN
    AerosolCoeff_File = 'AerosolCoeff.GOCART-GEOS5.nc4'
    Error_Status = CRTM_Init( (/Sensor_Id/), &
                              ChannelInfo, &
                              Aerosol_Model = 'GOCART-GEOS5', &
                              AerosolCoeff_Format = 'netCDF', &
                              AerosolCoeff_File = TRIM(AerosolCoeff_File), &
                              File_Path=COEFFICIENTS_PATH)
  END IF
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
            RTSolution_Total( n_Channels, N_PROFILES ), STAT=Allocate_Status )
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
  CALL CRTM_RTSolution_Create( RTSolution_Total, N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_Total)) ) THEN
    Message = 'Error allocating CRTM RTSolution structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ============================================================================




  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  CALL Load_Atm_Data()
  CALL Set_FineMode_Aerosols()
  ! ============================================================================




  ! ============================================================================
  ! 5. **** CALL THE CRTM AOD MODEL ****
  !
  Error_Status = CRTM_AOD( Atm        , &
                           ChannelInfo, &
                           RTSolution , &
                           AAOD = .TRUE., FINE_MODE = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM AOD Model with fine-mode AAOD'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! ...and the total AAOD for the consistency check
  Error_Status = CRTM_AOD( Atm             , &
                           ChannelInfo     , &
                           RTSolution_Total, &
                           AAOD = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM AOD Model for the total AAOD'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! Check 0 <= fine-mode <= total, and column fine-mode < total
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      IF ( ANY(RTSolution(l,m)%Layer_Optical_Depth < ZERO) .OR. &
           ANY(RTSolution(l,m)%Layer_Optical_Depth > &
               RTSolution_Total(l,m)%Layer_Optical_Depth*(ONE+FM_RTOL)) ) THEN
        WRITE( Message,'("Fine-mode AAOD outside [0,total] for channel ",i0,", profile ",i0)' ) &
               RTSolution(l,m)%Sensor_Channel, m
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP 1
      END IF
      IF ( SUM(RTSolution(l,m)%Layer_Optical_Depth) >= &
           SUM(RTSolution_Total(l,m)%Layer_Optical_Depth) ) THEN
        WRITE( Message,'("Column fine-mode AAOD not < total for channel ",i0,", profile ",i0)' ) &
               RTSolution(l,m)%Sensor_Channel, m
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP 1
      END IF
    END DO
  END DO
  Message = 'Fine-mode AAOD is within [0, total AAOD]'
  CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ! ============================================================================




  ! ============================================================================
  ! 6. **** OUTPUT THE RESULTS TO SCREEN ****
  !
  DO m = 1, N_PROFILES
    WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(Sensor_Id)
    DO l = 1, n_Channels
      WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
      CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
    END DO
  END DO
  ! ============================================================================

  ! ============================================================================
  ! 8. **** COMPARE RTSolution RESULTS TO SAVED VALUES ****
  !
  WRITE( *, '( /5x, "Comparing calculated results with saved ones..." )' )

  ! 8a. Create the output file if it does not exist
  ! -----------------------------------------------
  ! ...Generate a filename
  rts_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution.bin'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(rts_File) ) THEN
    Message = 'RTSolution save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write RTSolution structure to file
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating RTSolution save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

  ! 8b. Inquire the saved file
  ! --------------------------
  Error_Status = CRTM_RTSolution_InquireFile( rts_File, &
                                              n_Channels = n_l, &
                                              n_Profiles = n_m )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring RTSolution save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 8c. Compare the dimensions
  ! --------------------------
  IF ( n_l /= n_Channels .OR. n_m /= N_PROFILES ) THEN
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

  ! 8e. Read the saved data
  ! -----------------------
  Error_Status = CRTM_RTSolution_ReadFile( rts_File, rts, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading RTSolution save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 8f. Compare the structures
  ! --------------------------
  IF ( ALL(CRTM_RTSolution_Compare(RTSolution, rts)) ) THEN
    Message = 'RTSolution results are the same!'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'RTSolution results are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    ! Write the current RTSolution results to file
    rts_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution.bin'
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution, Quiet=.TRUE. )
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
  CALL CRTM_RTSolution_Destroy(RTSolution_Total)
  CALL CRTM_Atmosphere_Destroy(Atm)

  ! 9b. Deallocate the arrays
  ! -------------------------
  DEALLOCATE(RTSolution, RTSolution_Total, rts, STAT=Allocate_Status)
  ! ============================================================================

CONTAINS

  ! Fine-mode test setup for the GOCART-GEOS5 table: assign aerosol types to
  ! the Load_Atm_Data aerosol profiles, and set relative humidity profiles
  ! that reach 0.98, so the RH-dependent fmass1 of sea salt bin 2, nitrate
  ! bin 1, sulfate and organic carbon 2 is exercised
  SUBROUTINE Set_FineMode_Aerosols()
    ! GOCART-GEOS5 aerosol type IDs (AerosolCoeff.GOCART-GEOS5 table order)
    INTEGER, PARAMETER :: DUST_1           = 1
    INTEGER, PARAMETER :: DUST_2           = 2
    INTEGER, PARAMETER :: SEASALT_2        = 7
    INTEGER, PARAMETER :: ORGANIC_CARBON_2 = 12
    INTEGER, PARAMETER :: SULFATE_1        = 15
    INTEGER, PARAMETER :: NITRATE_1        = 17
    INTEGER :: k
    Atm(1)%Aerosol(1)%Type = DUST_1
    Atm(1)%Aerosol(2)%Type = SEASALT_2
    Atm(1)%Aerosol(3)%Type = SULFATE_1
    Atm(2)%Aerosol(1)%Type = NITRATE_1
    Atm(2)%Aerosol(2)%Type = ORGANIC_CARBON_2
    Atm(2)%Aerosol(3)%Type = DUST_2
    DO k = 1, N_LAYERS
      Atm(1)%Relative_Humidity(k) = 0.30_fp + 0.68_fp*REAL(k-1,fp)/REAL(N_LAYERS-1,fp)
      Atm(2)%Relative_Humidity(k) = 0.98_fp - 0.68_fp*REAL(k-1,fp)/REAL(N_LAYERS-1,fp)
    END DO
  END SUBROUTINE Set_FineMode_Aerosols

  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_FineMode_AAOD
