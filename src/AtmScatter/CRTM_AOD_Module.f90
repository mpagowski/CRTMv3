!
! CRTM_AOD_Module
!
! Module containing the CRTM Aerosol Optical Depth (AOD) at nadir
! functions.
!
!
! CREATION HISTORY:
!       Written by:     Quanhua Liu, 29-Jun-2010
!                       Quanhua.Liu@noaa.gov
!       Modified by     Yingtao Ma, 2020/6/11
!                       yingtao.ma@noaa.gov
!                       Implemented CMAQ aerosol
!       Modified by     Mariusz Pagowski, 2026/9/26
!                       mariusz.pagowski@colorado.edu
!                       Added optional AAOD argument to CRTM_AOD, CRTM_AOD_TL,
!                       CRTM_AOD_AD and CRTM_AOD_K to compute the absorption
!                       aerosol optical depth (AAOD)
!       Modified by     Mariusz Pagowski, 2026/10/22
!                       mariusz.pagowski@colorado.edu
!                       Added optional FINE_MODE argument to compute the fine-mode
!                       AOD/AAOD, with each aerosol type scaled by its fmass1(RH)

MODULE CRTM_AOD_Module


  ! ------------
  ! Module usage
  ! ------------
  USE Type_Kinds,               ONLY: fp
  USE Message_Handler,          ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,          ONLY: ZERO, ONE, &
                                      MAX_N_PHASE_ELEMENTS, &
                                      MAX_N_LEGENDRE_TERMS
  USE CRTM_Atmosphere_Define,   ONLY: CRTM_Atmosphere_type, &
                                      CRTM_Atmosphere_IsValid, &
                                      AerosolCoeff_BYPASS_AEROSOL
  USE CRTM_ChannelInfo_Define,  ONLY: CRTM_ChannelInfo_type, &
                                      CRTM_ChannelInfo_n_Channels
  USE CRTM_Options_Define,      ONLY: CRTM_Options_type
  USE CRTM_AtmOptics_Define,    ONLY: CRTM_AtmOptics_type      , &
                                      CRTM_AtmOptics_Associated, &
                                      CRTM_AtmOptics_Create    , &
                                      CRTM_AtmOptics_Destroy   , &
                                      CRTM_AtmOptics_Zero
  USE CRTM_AerosolScatter,      ONLY: CRTM_Compute_AerosolScatter, &
                                      CRTM_Compute_AerosolScatter_TL, &
                                      CRTM_Compute_AerosolScatter_AD
  USE CRTM_RTSolution_Define,   ONLY: CRTM_RTSolution_type, &
                                      CRTM_RTSolution_Associated
  USE CRTM_AerosolCoeff,        ONLY: CRTM_AerosolCoeff_IsLoaded, &
                                      AeroC
  USE AerosolCoeff_Define,      ONLY: AerosolCoeff_typeID_to_index

  ! Internal variable definition modules
  ! ...AerosolScatter
  USE ASvar_Define, ONLY: ASvar_type, &
                          ASvar_Associated, &
                          ASvar_Destroy   , &
                          ASvar_Create


  ! -----------------------
  ! Disable implicit typing
  ! -----------------------
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Public procedures
  PUBLIC :: CRTM_AOD
  PUBLIC :: CRTM_AOD_TL
  PUBLIC :: CRTM_AOD_AD
  PUBLIC :: CRTM_AOD_K

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


CONTAINS


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_AOD
!
! PURPOSE:
!       Function that calculates layer total optical depth profile at nadir.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_AOD( Atmosphere       , &
!                                ChannelInfo      , &
!                                RTSolution       , &
!                                Options = Options, &
!                                AAOD      = AAOD     , &
!                                FINE_MODE = FINE_MODE  )
!
! INPUTS:
!       Atmosphere:     Structure containing the Atmosphere data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Rank-1 (n_Profiles)
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelInfo:    Structure returned from the CRTM_Init() function
!                       that contains the satellite/sensor channel index
!                       information.
!                       UNITS:      N/A
!                       TYPE:       CRTM_ChannelInfo_type
!                       DIMENSION:  Rank-1 (n_Sensors)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       RTSolution:     Structure containing the layer aerosol optical
!                       profile for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       Options:        Options structure containing the optional arguments
!                       for the CRTM.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Options_type
!                       DIMENSION:  Same as input Atmosphere structure
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       AAOD:           Set to .TRUE. to compute the aerosol absorption optical
!                       depth (AAOD) instead of the AOD. The AAOD is returned
!                       in the same Layer_Optical_Depth components that
!                       otherwise hold the AOD. Aerosol scattering is switched
!                       off, so the aerosol extinction coefficient is replaced
!                       by the absorption coefficient, ke*(1-w).
!                       If not specified, or .FALSE., the AOD is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       FINE_MODE:      Set to .TRUE. to compute the fine-mode AOD (or AAOD,
!                       if AAOD is also .TRUE.). The concentration of each
!                       aerosol type is scaled by its fmass1(RH), the fraction
!                       of its mass in particles with diameter < 1 um,
!                       linearly interpolated to the layer relative humidity.
!                       Only available for the GOCART-GEOS5 aerosol scheme.
!                       If not specified, or .FALSE., the total AOD (or AAOD)
!                       is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       - Many of the components of the Options optional input structure
!         are not used in this function. Consult the CRTM User Guide for
!         which Options components are usable for AOD calculations.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_AOD( &
    Atmosphere , &  ! Input, M
    ChannelInfo, &  ! Input, M
    RTSolution , &  ! Output, L x M
    Options    , &  ! Optional input, M
    AAOD       , &  ! Optional input
    FINE_MODE  ) &  ! Optional input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),        INTENT(IN)     :: Atmosphere(:)     ! M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)    ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)   ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)        ! M
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: AAOD
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: FINE_MODE
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_AOD'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    LOGICAL :: Compute_AAOD
    LOGICAL :: Compute_FM
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: ln
    ! Component variables
    TYPE(CRTM_AtmOptics_type) :: AtmOptics
    TYPE(ASVar_type) :: ASvar
    ! Fine-mode variables
    TYPE(CRTM_Atmosphere_type) :: Atm           ! Atmosphere with fine-mode concentrations
    REAL(fp), ALLOCATABLE :: f(:,:), f_RH(:,:)  ! fmass1 and d(fmass1)/d(RH), L x n_Aerosols


    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS


    ! Compute absorption AOD (AAOD) instead of AOD if requested
    Compute_AAOD = .FALSE.
    IF ( PRESENT(AAOD) ) Compute_AAOD = AAOD
    ! Compute the fine-mode AOD/AAOD if requested
    Compute_FM = .FALSE.
    IF ( PRESENT(FINE_MODE) ) Compute_FM = FINE_MODE


    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN


    ! Check the number of channels
    IF ( SIZE(RTSolution,DIM=1) < n_Channels ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure array too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    ! ...Number of atmospheric profiles.
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(RTSolution,DIM=2) /= n_Profiles ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for RTSolution argument.'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF


    ! Check the RTSolution structure has been allocated
    IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) ) THEN
      Error_Status = FAILURE
      Message = 'RTSolution output structure components have not been allocated'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles

      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF


      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input

        ! Check whether to skip this profile
        IF ( Options(m)%Skip_Profile ) CYCLE Profile_Loop
      END IF


      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF


      ! Check the RTSolution layer dimension
      IF ( ANY(RTSolution(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF


      ! Allocate AtmOptics based on Atmosphere dimension
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS, &
                                  AeroC%N_PHASE_ELEMENTS  )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structure for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms = 4
      ! ...AAOD: without scattering the aerosol extinction coefficient
      !    is replaced by the absorption coefficient, ke*(1-w), in
      !    CRTM_Compute_AerosolScatter, so the optical depth is the AAOD
      AtmOptics%Include_Scattering = .NOT. Compute_AAOD


      ! Allocate the aerosol scattering internal variable if necessary
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS    , &
                           MAX_N_PHASE_ELEMENTS    , &
                           Atmosphere(m)%n_Layers  , &
                           Atmosphere(m)%n_Aerosols  )
      END IF


      ! Fine-mode: scale the aerosol concentrations by fmass1
      Atm = Atmosphere(m)
      IF ( Compute_FM .AND. Atmosphere(m)%n_Aerosols > 0 ) THEN
        Error_Status = FineMode_Factors( Atmosphere(m), f, f_RH )
        IF ( Error_Status /= SUCCESS ) THEN
          WRITE( Message,'("Error computing fine-mode factors for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
        CALL FineMode_Scale( Atmosphere(m), f, Atm )
      END IF


      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for channel(l)/sensor(n) count
      ln = 0

      Sensor_Loop: DO n = 1, n_Sensors

        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index


        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )


          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status = CRTM_Compute_AerosolScatter( Atm          , &  ! Input
                                                        SensorIndex  , &  ! Input
                                                        ChannelIndex , &  ! Input
                                                        AtmOptics    , &  ! In/Output
                                                        ASVar          )  ! Internal variable output
            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF
          END IF

          ! Save the nadir optical depth
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth

          ! Save the Backscat_Coefficient, only non-zero for GOCART-GEOS5
          ! CD: double check if should be divided by 2PI
          RTSolution(ln,m)%Backscat_Coefficient(1:Atmosphere(m)%n_Layers) = AtmOptics%Backscat_Coefficient
          
        END DO Channel_Loop

      END DO Sensor_Loop


      ! Deallocate local sensor independent data structures
      CALL CRTM_AtmOptics_Destroy( AtmOptics )

    END DO Profile_Loop

  END FUNCTION CRTM_AOD


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_AOD_TL
!
! PURPOSE:
!       Function that calculates tangent-linear layer total optical depth.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_AOD_TL( Atmosphere       , &
!                                   Atmosphere_TL    , &
!                                   ChannelInfo      , &
!                                   RTSolution       , &
!                                   RTSolution_TL    , &
!                                   Options = Options, &
!                                   AAOD      = AAOD     , &
!                                   FINE_MODE = FINE_MODE  )
!
! INPUTS:
!       Atmosphere:     Structure containing the Atmosphere data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Rank-1 (n_Profiles)
!                       ATTRIBUTES: INTENT(IN)
!
!       Atmosphere_TL:  Structure containing the tangent-linear Atmosphere data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Same as input Atmosphere structure
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelInfo:    Structure returned from the CRTM_Init() function
!                       that contains the satellite/sensor channel index
!                       information.
!                       UNITS:      N/A
!                       TYPE:       CRTM_ChannelInfo_type
!                       DIMENSION:  Rank-1 (n_Sensors)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       RTSolution:     Structure containing the layer aerosol optical
!                       profile for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       RTSolution_TL:  Structure containing the tangent-linear aerosol
!                       optical depth profile for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Same as RTSolution output
!                       ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       Options:        Options structure containing the optional arguments
!                       for the CRTM.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Options_type
!                       DIMENSION:  Same as input Atmosphere structure
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       AAOD:           Set to .TRUE. to compute the aerosol absorption optical
!                       depth (AAOD) instead of the AOD. The AAOD is returned
!                       in the same Layer_Optical_Depth components that
!                       otherwise hold the AOD. Aerosol scattering is switched
!                       off, so the aerosol extinction coefficient is replaced
!                       by the absorption coefficient, ke*(1-w).
!                       If not specified, or .FALSE., the AOD is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       FINE_MODE:      Set to .TRUE. to compute the fine-mode AOD (or AAOD,
!                       if AAOD is also .TRUE.). The concentration of each
!                       aerosol type is scaled by its fmass1(RH), the fraction
!                       of its mass in particles with diameter < 1 um,
!                       linearly interpolated to the layer relative humidity.
!                       Only available for the GOCART-GEOS5 aerosol scheme.
!                       If not specified, or .FALSE., the total AOD (or AAOD)
!                       is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       - Many of the components of the Options optional input structure
!         are not used in this function. Consult the CRTM User Guide for
!         which Options components are usable for AOD calculations.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_AOD_TL( &
    Atmosphere   , &  ! Input, M
    Atmosphere_TL, &  ! Input, M
    ChannelInfo  , &  ! Input, M
    RTSolution   , &  ! Output, L x M
    RTSolution_TL, &  ! Output, L x M
    Options      , &  ! Optional FWD input, M
    AAOD         , &  ! Optional input
    FINE_MODE    ) &  ! Optional input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),        INTENT(IN)     :: Atmosphere(:)      ! M
    TYPE(CRTM_Atmosphere_type),        INTENT(IN)     :: Atmosphere_TL(:)   ! M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)     ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)    ! L x M
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution_TL(:,:) ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)         ! M
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: AAOD
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: FINE_MODE
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_AOD_TL'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    LOGICAL :: Compute_AAOD
    LOGICAL :: Compute_FM
    INTEGER :: Status_FWD, Status_TL
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: ln
    ! Component variables
    TYPE(CRTM_AtmOptics_type)   :: AtmOptics, AtmOptics_TL
    TYPE(ASVar_type) :: ASvar
    ! Fine-mode variables
    TYPE(CRTM_Atmosphere_type) :: Atm, Atm_TL   ! Atmospheres with fine-mode concentrations
    REAL(fp), ALLOCATABLE :: f(:,:), f_RH(:,:)  ! fmass1 and d(fmass1)/d(RH), L x n_Aerosols


    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS


    ! Compute absorption AOD (AAOD) instead of AOD if requested
    Compute_AAOD = .FALSE.
    IF ( PRESENT(AAOD) ) Compute_AAOD = AAOD
    ! Compute the fine-mode AOD/AAOD if requested
    Compute_FM = .FALSE.
    IF ( PRESENT(FINE_MODE) ) Compute_FM = FINE_MODE


    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN


    ! Check the number of channels
    IF ( SIZE(RTSolution,   DIM=1) < n_Channels .OR. &
         SIZE(RTSolution_TL,DIM=1) < n_Channels      ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure arrays too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(Atmosphere_TL)       /= n_Profiles .OR. &
         SIZE(RTSolution,   DIM=2) /= n_Profiles .OR. &
         SIZE(RTSolution_TL,DIM=2) /= n_Profiles      ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for input arguments.'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF


    ! Check the RTSolution structures have been allocated
    IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) .OR. &
         ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_TL)) ) THEN
      Error_Status = FAILURE
      Message = 'RTSolution output structure components have not been allocated'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles


      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF


      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
      END IF


      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF


      ! Check the RTSolution layer dimensions
      IF ( ANY(RTSolution(:,m)%n_Layers    < Atmosphere(m)%n_Layers) .OR. &
           ANY(RTSolution_TL(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF


      ! Allocate AtmOptics based on Atmosphere dimensions
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS  , &
                                  AeroC%N_PHASE_ELEMENTS    )
      CALL CRTM_AtmOptics_Create( AtmOptics_TL, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS  , &
                                  AeroC%N_PHASE_ELEMENTS    )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) .OR. &
           .NOT. CRTM_AtmOptics_Associated( Atmoptics_TL ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structures for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms    = 4
      AtmOptics_TL%n_Legendre_Terms = AtmOptics%n_Legendre_Terms
      ! ...AAOD: without scattering the aerosol extinction coefficient
      !    is replaced by the absorption coefficient, ke*(1-w), in
      !    CRTM_Compute_AerosolScatter, so the optical depth is the AAOD
      AtmOptics%Include_Scattering    = .NOT. Compute_AAOD
      AtmOptics_TL%Include_Scattering = .NOT. Compute_AAOD


      ! Allocate the aerosol scattering internal variable if necessary
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS    , &
                           MAX_N_PHASE_ELEMENTS    , &
                           Atmosphere(m)%n_Layers  , &
                           Atmosphere(m)%n_Aerosols  )
      END IF


      ! Fine-mode: scale the aerosol concentrations by fmass1
      Atm = Atmosphere(m)
      IF ( Compute_FM .AND. Atmosphere(m)%n_Aerosols > 0 ) THEN
        Error_Status = FineMode_Factors( Atmosphere(m), f, f_RH )
        IF ( Error_Status /= SUCCESS ) THEN
          WRITE( Message,'("Error computing fine-mode factors for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
        CALL FineMode_Scale( Atmosphere(m), f, Atm )
      END IF
      Atm_TL = Atmosphere_TL(m)
      IF ( Compute_FM .AND. Atmosphere(m)%n_Aerosols > 0 ) &
        CALL FineMode_Scale_TL( Atmosphere(m), Atmosphere_TL(m), f, f_RH, Atm_TL )


      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for channel(l)/sensor(n) count
      ln = 0

      Sensor_Loop: DO n = 1, n_Sensors

        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index


        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)
          RTSolution_TL(ln,m)%Sensor_Id        = RTSolution(ln,m)%Sensor_Id
          RTSolution_TL(ln,m)%WMO_Satellite_Id = RTSolution(ln,m)%WMO_Satellite_Id
          RTSolution_TL(ln,m)%WMO_Sensor_Id    = RTSolution(ln,m)%WMO_Sensor_Id
          RTSolution_TL(ln,m)%Sensor_Channel   = RTSolution(ln,m)%Sensor_Channel


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
          CALL CRTM_AtmOptics_Zero( AtmOptics_TL )


          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Status_FWD = CRTM_Compute_AerosolScatter( Atm          , &  ! Input
                                                      SensorIndex  , &  ! Input
                                                      ChannelIndex , &  ! Input
                                                      AtmOptics    , &  ! In/Output
                                                      ASVar          )  ! Internal variable output
            Status_TL = CRTM_Compute_AerosolScatter_TL( Atm             , &  ! FWD Input
                                                        AtmOptics       , &  ! FWD Input
                                                        Atm_TL          , &  ! TL  Input
                                                        SensorIndex     , &  ! Input
                                                        ChannelIndex    , &  ! Input
                                                        AtmOptics_TL    , &  ! TL  Output
                                                        ASVar             )  ! Internal variable
            IF ( Status_FWD /= SUCCESS .OR. Status_TL /= SUCCESS) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
              RETURN
            END IF

          END IF

          ! Save the nadir optical depths
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth
          RTSolution_TL(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics_TL%Optical_Depth

          ! Save the backscatter coefficient, only non-zero for GOCART-GEOS5
          ! CD: double check if should be divided by 2PI
          RTSolution(ln,m)%Backscat_Coefficient(1:Atmosphere(m)%n_Layers) = AtmOptics%Backscat_Coefficient
          RTSolution_TL(ln,m)%Backscat_Coefficient(1:Atmosphere(m)%n_Layers) = AtmOptics_TL%Backscat_Coefficient

        END DO Channel_Loop

      END DO Sensor_Loop


      ! Deallocate local sensor independent data structures
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_AtmOptics_Destroy( AtmOptics_TL )

    END DO Profile_Loop

  END FUNCTION CRTM_AOD_TL


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_AOD_AD
!
! PURPOSE:
!       Function that calculates the adjoint nadir aerosol optical depth.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_AOD_AD( Atmosphere       , &
!                                   RTSolution_AD    , &
!                                   ChannelInfo      , &
!                                   RTSolution       , &
!                                   Atmosphere_AD    , &
!                                   Options = Options, &
!                                   AAOD      = AAOD     , &
!                                   FINE_MODE = FINE_MODE  )
!
! INPUTS:
!       Atmosphere:     Structure containing the Atmosphere data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Rank-1 (n_Profiles)
!                       ATTRIBUTES: INTENT(IN)
!
!       RTSolution_AD:  Structure containing the RT solution adjoint inputs.
!                       **NOTE: On EXIT from this function, the contents of
!                               this structure may be modified (e.g. set to
!                               zero.)
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       ChannelInfo:    Structure returned from the CRTM_Init() function
!                       that contains the satellite/sensor channel index
!                       information.
!                       UNITS:      N/A
!                       TYPE:       CRTM_ChannelInfo_type
!                       DIMENSION:  Rank-1 (n_Sensors)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       RTSolution:     Structure containing the soluition to the RT equation
!                       for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       Atmosphere_AD:  Structure containing the adjoint Atmosphere data.
!                       **NOTE: On ENTRY to this function, the contents of
!                               this structure should be defined (e.g.
!                               initialized to some value based on the
!                               position of this function in the call chain.)
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Same as input Atmosphere argument
!                       ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       Options:        Options structure containing the optional arguments
!                       for the CRTM.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Options_type
!                       DIMENSION:  Same as input Atmosphere structure
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       AAOD:           Set to .TRUE. to compute the aerosol absorption optical
!                       depth (AAOD) instead of the AOD. The AAOD is returned
!                       in the same Layer_Optical_Depth components that
!                       otherwise hold the AOD. Aerosol scattering is switched
!                       off, so the aerosol extinction coefficient is replaced
!                       by the absorption coefficient, ke*(1-w).
!                       If not specified, or .FALSE., the AOD is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       FINE_MODE:      Set to .TRUE. to compute the fine-mode AOD (or AAOD,
!                       if AAOD is also .TRUE.). The concentration of each
!                       aerosol type is scaled by its fmass1(RH), the fraction
!                       of its mass in particles with diameter < 1 um,
!                       linearly interpolated to the layer relative humidity.
!                       Only available for the GOCART-GEOS5 aerosol scheme.
!                       If not specified, or .FALSE., the total AOD (or AAOD)
!                       is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       - Many of the components of the Options optional input structure
!         are not used in this function. Consult the CRTM User Guide for
!         which Options components are usable for AOD calculations.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_AOD_AD( &
    Atmosphere   , &  ! Input, M
    RTSolution_AD, &  ! Input, M
    ChannelInfo  , &  ! Input, M
    RTSolution   , &  ! Output, L x M
    Atmosphere_AD, &  ! Output, L x M
    Options      , &  ! Optional input, M
    AAOD         , &  ! Optional input
    FINE_MODE    ) &  ! Optional input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),        INTENT(IN)     :: Atmosphere(:)      ! M
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution_AD(:,:) ! L x M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)     ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)    ! L x M
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT) :: Atmosphere_AD(:)   ! M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)         ! M
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: AAOD
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: FINE_MODE
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_AOD_AD'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    LOGICAL :: Compute_AAOD
    LOGICAL :: Compute_FM
    INTEGER :: Status_FWD, Status_AD
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: na
    INTEGER :: ln
    ! Component variables
    TYPE(CRTM_AtmOptics_type) :: AtmOptics, AtmOptics_AD
    TYPE(ASVar_type) :: ASvar
    ! Fine-mode variables
    TYPE(CRTM_Atmosphere_type) :: Atm           ! Atmosphere with fine-mode concentrations
    REAL(fp), ALLOCATABLE :: f(:,:), f_RH(:,:)  ! fmass1 and d(fmass1)/d(RH), L x n_Aerosols
    REAL(fp), ALLOCATABLE :: Conc_AD_in(:,:)     ! Input AD concentrations, L x n_Aerosols


    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS


    ! Compute absorption AOD (AAOD) instead of AOD if requested
    Compute_AAOD = .FALSE.
    IF ( PRESENT(AAOD) ) Compute_AAOD = AAOD
    ! Compute the fine-mode AOD/AAOD if requested
    Compute_FM = .FALSE.
    IF ( PRESENT(FINE_MODE) ) Compute_FM = FINE_MODE


    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN


    ! Check the number of channels
    IF ( SIZE(RTSolution,   DIM=1) < n_Channels .OR. &
         SIZE(RTSolution_AD,DIM=1) < n_Channels      ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure arrays too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(Atmosphere_AD)       /= n_Profiles .OR. &
         SIZE(RTSolution,   DIM=2) /= n_Profiles .OR. &
         SIZE(RTSolution_AD,DIM=2) /= n_Profiles      ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for input arguments.'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF


    ! Check the RTSolution structures have been allocated
    IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) .OR. &
         ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_AD)) ) THEN
      Error_Status = FAILURE
      Message = 'RTSolution output structure components have not been allocated'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles


      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF


      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
      END IF


      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF


      ! Check the RTSolution layer dimensions
      IF ( ANY(RTSolution(:,m)%n_Layers    < Atmosphere(m)%n_Layers) .OR. &
           ANY(RTSolution_AD(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF


      ! Allocate AtmOptics based on Atmosphere dimensions
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS  , &
                                  AeroC%N_PHASE_ELEMENTS    )
      CALL CRTM_AtmOptics_Create( AtmOptics_AD, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS  , &
                                  AeroC%N_PHASE_ELEMENTS    )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) .OR. &
           .NOT. CRTM_AtmOptics_Associated( Atmoptics_AD ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structures for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms    = 4
      AtmOptics_AD%n_Legendre_Terms = AtmOptics%n_Legendre_Terms
      ! ...AAOD: without scattering the aerosol extinction coefficient
      !    is replaced by the absorption coefficient, ke*(1-w), in
      !    CRTM_Compute_AerosolScatter, so the optical depth is the AAOD
      AtmOptics%Include_Scattering    = .NOT. Compute_AAOD
      AtmOptics_AD%Include_Scattering = .NOT. Compute_AAOD


      ! Allocate the aerosol scattering internal variable if necessary
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS    , &
                           MAX_N_PHASE_ELEMENTS    , &
                           Atmosphere(m)%n_Layers  , &
                           Atmosphere(m)%n_Aerosols  )
      END IF


      ! Fine-mode: scale the aerosol concentrations by fmass1
      Atm = Atmosphere(m)
      IF ( Compute_FM .AND. Atmosphere(m)%n_Aerosols > 0 ) THEN
        Error_Status = FineMode_Factors( Atmosphere(m), f, f_RH )
        IF ( Error_Status /= SUCCESS ) THEN
          WRITE( Message,'("Error computing fine-mode factors for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
        CALL FineMode_Scale( Atmosphere(m), f, Atm )
      END IF


      ! Copy over atmosphere info to adjoint output
      ! ...Climatology
      Atmosphere_AD(m)%Climatology = Atmosphere(m)%Climatology
      ! ...Absorber info
      Atmosphere_AD(m)%Absorber_Id    = Atmosphere(m)%Absorber_Id
      Atmosphere_AD(m)%Absorber_Units = Atmosphere(m)%Absorber_Units
      ! ...Aerosol info
      DO na = 1, Atmosphere(m)%n_Aerosols
        Atmosphere_AD(m)%Aerosol(na)%Type = Atmosphere(m)%Aerosol(na)%Type
      END DO


      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for sensor(n)/channel(l) count
      ln = 0

      Sensor_Loop: DO n = 1, n_Sensors


        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index


        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)
          RTSolution_AD(ln,m)%Sensor_Id        = RTSolution(ln,m)%Sensor_Id
          RTSolution_AD(ln,m)%WMO_Satellite_Id = RTSolution(ln,m)%WMO_Satellite_Id
          RTSolution_AD(ln,m)%WMO_Sensor_Id    = RTSolution(ln,m)%WMO_Sensor_Id
          RTSolution_AD(ln,m)%Sensor_Channel   = RTSolution(ln,m)%Sensor_Channel


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
          CALL CRTM_AtmOptics_Zero( AtmOptics_AD )
          AtmOptics_AD%Optical_Depth = RTSolution_AD(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers)


          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Status_FWD = CRTM_Compute_AerosolScatter( Atm          , &  ! Input
                                                      SensorIndex  , &  ! Input
                                                      ChannelIndex , &  ! Input
                                                      AtmOptics    , &  ! In/Output
                                                      ASVar          )  ! Internal variable output
            ! ...Fine-mode: accumulate the adjoint of the scaled concentrations
            !    from zero, then map it back to the input concentrations and RH
            IF ( Compute_FM ) CALL FineMode_Save_AD( Atmosphere_AD(m), Conc_AD_in )
            Status_AD = CRTM_Compute_AerosolScatter_AD( Atm             , &  ! FWD Input
                                                        AtmOptics       , &  ! FWD Input
                                                        AtmOptics_AD    , &  ! AD  Input
                                                        SensorIndex     , &  ! Input
                                                        ChannelIndex    , &  ! Input
                                                        Atmosphere_AD(m), &  ! AD  Output
                                                        ASVar             )  ! Internal variable input
            IF ( Compute_FM ) &
              CALL FineMode_Scale_AD( Atmosphere(m), f, f_RH, Conc_AD_in, Atmosphere_AD(m) )
            IF ( Status_FWD /= SUCCESS .OR. Status_AD /= SUCCESS) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF
          END IF


          ! Save the nadir optical depths
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth

        END DO Channel_Loop

      END DO Sensor_Loop


      ! Deallocate local sensor independent data structures
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_AtmOptics_Destroy( AtmOptics_AD )

    END DO Profile_Loop

  END FUNCTION CRTM_AOD_AD


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_AOD_K
!
! PURPOSE:
!       Function that calculates the K-matrix nadir aerosol optical depth.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_AOD_K( Atmosphere        , &
!                                  RTSolution_K      , &
!                                  ChannelInfo       , &
!                                  RTSolution        , &
!                                  Atmosphere_K      , &
!                                  Opttions = Options, &
!                                  AAOD      = AAOD     , &
!                                  FINE_MODE = FINE_MODE  )
!
! INPUTS:
!       Atmosphere:     Structure containing the Atmosphere data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Rank-1 (n_Profiles)
!                       ATTRIBUTES: INTENT(IN)
!
!       RTSolution_K:   Structure containing the aerosol optical depth
!                       profile K-matrix input.
!                       **NOTE: On EXIT from this function, the contents of
!                               this structure may be modified (e.g. set to
!                               zero.)
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       ChannelInfo:    Structure returned from the CRTM_Init() function
!                       that contains the satellite/sensor channel index
!                       information.
!                       UNITS:      N/A
!                       TYPE:       CRTM_ChannelInfo_type
!                       DIMENSION:  Rank-1 (n_Sensors)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       RTSolution:     Structure containing the layer aerosol optical
!                       depth profile for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       Atmosphere_K:   Structure containing the K-matrix Atmosphere data.
!                       **NOTE: On ENTRY to this function, the contents of
!                               this structure should be defined (e.g.
!                               initialized to some value based on the
!                               position of this function in the call chain.)
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Same as input RTSolution_K argument
!                       ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       Options:        Options structure containing the optional arguments
!                       for the CRTM.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Options_type
!                       DIMENSION:  Same as input Atmosphere structure
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       AAOD:           Set to .TRUE. to compute the aerosol absorption optical
!                       depth (AAOD) instead of the AOD. The AAOD is returned
!                       in the same Layer_Optical_Depth components that
!                       otherwise hold the AOD. Aerosol scattering is switched
!                       off, so the aerosol extinction coefficient is replaced
!                       by the absorption coefficient, ke*(1-w).
!                       If not specified, or .FALSE., the AOD is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       FINE_MODE:      Set to .TRUE. to compute the fine-mode AOD (or AAOD,
!                       if AAOD is also .TRUE.). The concentration of each
!                       aerosol type is scaled by its fmass1(RH), the fraction
!                       of its mass in particles with diameter < 1 um,
!                       linearly interpolated to the layer relative humidity.
!                       Only available for the GOCART-GEOS5 aerosol scheme.
!                       If not specified, or .FALSE., the total AOD (or AAOD)
!                       is computed.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       - Many of the components of the Options optional input structure
!         are not used in this function. Consult the CRTM User Guide for
!         which Options components are usable for AOD calculations.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_AOD_K ( &
    Atmosphere  , &  ! Input, M
    RTSolution_K, &  ! Input, M
    ChannelInfo , &  ! Input, M
    RTSolution  , &  ! Output, L x M
    Atmosphere_K, &  ! Output, L x M
    Options     , &  ! Optional input, M
    AAOD        , &  ! Optional input
    FINE_MODE   ) &  ! Optional input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),        INTENT(IN)     :: Atmosphere(:)      ! M
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution_K(:,:)  ! L x M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)     ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)    ! L x M
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT) :: Atmosphere_K(:,:)  ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)         ! M
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: AAOD
    LOGICAL,                 OPTIONAL, INTENT(IN)     :: FINE_MODE
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_AOD_K'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    LOGICAL :: Compute_AAOD
    LOGICAL :: Compute_FM
    INTEGER :: Status_FWD, Status_K
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: na
    INTEGER :: ln
    ! Component variables
    TYPE(CRTM_AtmOptics_type) :: AtmOptics, AtmOptics_K
    TYPE(ASVar_type) :: ASvar
    ! Fine-mode variables
    TYPE(CRTM_Atmosphere_type) :: Atm           ! Atmosphere with fine-mode concentrations
    REAL(fp), ALLOCATABLE :: f(:,:), f_RH(:,:)  ! fmass1 and d(fmass1)/d(RH), L x n_Aerosols
    REAL(fp), ALLOCATABLE :: Conc_K_in(:,:)     ! Input K concentrations, L x n_Aerosols


    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS


    ! Compute absorption AOD (AAOD) instead of AOD if requested
    Compute_AAOD = .FALSE.
    IF ( PRESENT(AAOD) ) Compute_AAOD = AAOD
    ! Compute the fine-mode AOD/AAOD if requested
    Compute_FM = .FALSE.
    IF ( PRESENT(FINE_MODE) ) Compute_FM = FINE_MODE


    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN


    ! Check spectral arrays
    IF ( SIZE(RTSolution  ,DIM=1) < n_Channels .OR. &
         SIZE(Atmosphere_K,DIM=1) < n_Channels .OR. &
         SIZE(RTSolution_K,DIM=1) < n_Channels      ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("RTSolution and K-matrix (Atm,RT) structure arrays too small (",&
             &2(i0,","),i0,") for the number of requested channels (",i0,")")') &
             SIZE(RTSolution  ,DIM=1), &
             SIZE(Atmosphere_K,DIM=1), &
             SIZE(RTSolution_K,DIM=1), &
             n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    ! ...Number of atmospheric profiles.
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(RTSolution_K,DIM=2) /= n_Profiles .OR. &
         SIZE(Atmosphere_K,DIM=2) /= n_Profiles .OR. &
         SIZE(RTSolution  ,DIM=2) /= n_Profiles      ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for input arguments.'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE( Options ) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF


    ! Check the RTSolution structures have been allocated
    IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) .OR. &
         ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_K)) ) THEN
      Error_Status = FAILURE
      Message = 'RTSolution output structure components have not been allocated'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles


      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF


      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
      END IF


      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF


      ! Check the RTSolution layer dimensions
      IF ( ANY(RTSolution(:,m)%n_Layers   < Atmosphere(m)%n_Layers) .OR. &
           ANY(RTSolution_K(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF


      ! Allocate AtmOptics based on Atmosphere dimensions
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS  , &
                                  AeroC%N_PHASE_ELEMENTS    )
      CALL CRTM_AtmOptics_Create( AtmOptics_K, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS  , &
                                  AeroC%N_PHASE_ELEMENTS    )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) .OR. &
           .NOT. CRTM_AtmOptics_Associated( Atmoptics_K ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structures for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms    = 4
      AtmOptics_K%n_Legendre_Terms = AtmOptics%n_Legendre_Terms
      ! ...AAOD: without scattering the aerosol extinction coefficient
      !    is replaced by the absorption coefficient, ke*(1-w), in
      !    CRTM_Compute_AerosolScatter, so the optical depth is the AAOD
      AtmOptics%Include_Scattering    = .NOT. Compute_AAOD
      AtmOptics_K%Include_Scattering  = .NOT. Compute_AAOD


      ! Allocate the aerosol scattering internal variable if necessary
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS    , &
                           MAX_N_PHASE_ELEMENTS    , &
                           Atmosphere(m)%n_Layers  , &
                           Atmosphere(m)%n_Aerosols  )
      END IF


      ! Fine-mode: scale the aerosol concentrations by fmass1
      Atm = Atmosphere(m)
      IF ( Compute_FM .AND. Atmosphere(m)%n_Aerosols > 0 ) THEN
        Error_Status = FineMode_Factors( Atmosphere(m), f, f_RH )
        IF ( Error_Status /= SUCCESS ) THEN
          WRITE( Message,'("Error computing fine-mode factors for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
        CALL FineMode_Scale( Atmosphere(m), f, Atm )
      END IF


      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for sensor(n)/channel(l) count
      ln = 0

      Sensor_Loop: DO n = 1, n_Sensors


        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index


        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)
          RTSolution_K(ln,m)%Sensor_Id        = RTSolution(ln,m)%Sensor_Id
          RTSolution_K(ln,m)%WMO_Satellite_Id = RTSolution(ln,m)%WMO_Satellite_Id
          RTSolution_K(ln,m)%WMO_Sensor_Id    = RTSolution(ln,m)%WMO_Sensor_Id
          RTSolution_K(ln,m)%Sensor_Channel   = RTSolution(ln,m)%Sensor_Channel


          ! Copy over atmosphere info to k-matrix output
          ! ...Climatology
          Atmosphere_K(ln,m)%Climatology = Atmosphere(m)%Climatology
          ! ...Absorber info
          Atmosphere_K(ln,m)%Absorber_Id    = Atmosphere(m)%Absorber_Id
          Atmosphere_K(ln,m)%Absorber_Units = Atmosphere(m)%Absorber_Units
          ! ...Aerosol info
          DO na = 1, Atmosphere(m)%n_Aerosols
            Atmosphere_K(ln,m)%Aerosol(na)%Type = Atmosphere(m)%Aerosol(na)%Type
          END DO


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
          CALL CRTM_AtmOptics_Zero( AtmOptics_K )
          AtmOptics_K%Optical_Depth = RTSolution_K(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers)


          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Status_FWD = CRTM_Compute_AerosolScatter( Atm          , &  ! Input
                                                      SensorIndex  , &  ! Input
                                                      ChannelIndex , &  ! Input
                                                      AtmOptics    , &  ! In/Output
                                                      ASVar          )  ! Internal variable output
            ! ...Fine-mode: accumulate the K-matrix of the scaled concentrations
            !    from zero, then map it back to the input concentrations and RH
            IF ( Compute_FM ) CALL FineMode_Save_AD( Atmosphere_K(ln,m), Conc_K_in )
            Status_K = CRTM_Compute_AerosolScatter_AD( Atm               , &  ! FWD Input
                                                       AtmOptics         , &  ! FWD Input
                                                       AtmOptics_K       , &  ! AD  Input
                                                       SensorIndex       , &  ! Input
                                                       ChannelIndex      , &  ! Input
                                                       Atmosphere_K(ln,m), &  ! AD  Output
                                                       ASVar               )  ! Internal variable input
            IF ( Compute_FM ) &
              CALL FineMode_Scale_AD( Atmosphere(m), f, f_RH, Conc_K_in, Atmosphere_K(ln,m) )
            IF ( Status_FWD /= SUCCESS .OR. Status_K /= SUCCESS) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF
          END IF


          ! Save the nadir optical depths
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth

        END DO Channel_Loop

      END DO Sensor_Loop


      ! Deallocate local sensor independent data structures
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_AtmOptics_Destroy( AtmOptics_K )

    END DO Profile_Loop
   END FUNCTION CRTM_AOD_K


!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

  ! -------------------------------------------------------------
  ! Function to compute the fine-mode factors, f = fmass1(RH), and
  ! their derivatives with respect to RH, f_RH, for every aerosol
  ! and layer. fmass1 is the AerosolCoeff fraction of the aerosol
  ! mass in particles with diameter < 1 um.
  ! -------------------------------------------------------------
  FUNCTION FineMode_Factors( &
    Atm , &  ! Input
    f   , &  ! Output, L x n_Aerosols
    f_RH) &  ! Output, L x n_Aerosols
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)  :: Atm
    REAL(fp), ALLOCATABLE,      INTENT(OUT) :: f(:,:)
    REAL(fp), ALLOCATABLE,      INTENT(OUT) :: f_RH(:,:)
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'FineMode_Factors'
    ! Local variables
    CHARACTER(ML) :: Message
    INTEGER :: k, na, it

    Error_Status = SUCCESS
    ALLOCATE( f(Atm%n_Layers, Atm%n_Aerosols), f_RH(Atm%n_Layers, Atm%n_Aerosols) )
    f    = ONE
    f_RH = ZERO

    IF ( TRIM(AeroC%Scheme) /= 'GOCART-GEOS5' ) THEN
      Error_Status = FAILURE
      Message = 'FINE_MODE is only available for the GOCART-GEOS5 aerosol scheme, not '//&
                TRIM(AeroC%Scheme)
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF

    Aerosol_Loop: DO na = 1, Atm%n_Aerosols
      IF ( Atm%Aerosol(na)%Type == AerosolCoeff_BYPASS_AEROSOL ) CYCLE Aerosol_Loop
      it = AerosolCoeff_typeID_to_index( AeroC, Atm%Aerosol(na)%Type )
      IF ( it < 1 ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Invalid aerosol type ",i0," for aerosol #",i0)' ) &
               Atm%Aerosol(na)%Type, na
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      DO k = 1, Atm%Aerosol(na)%n_Layers
        CALL Interpolate_fmass1( AeroC%RH, AeroC%fmass1(:,it), Atm%Relative_Humidity(k), &
                                 f(k,na), f_RH(k,na) )
      END DO
    END DO Aerosol_Loop

  END FUNCTION FineMode_Factors


  ! -------------------------------------------------------------
  ! Subroutine to linearly interpolate fmass1 in relative humidity.
  ! Outside the RH range of the table fmass1 is held constant, so
  ! its derivative is zero.
  ! -------------------------------------------------------------
  SUBROUTINE Interpolate_fmass1( &
    RH_LUT    , &  ! Input
    fmass1_LUT, &  ! Input
    RH        , &  ! Input
    f         , &  ! Output
    f_RH        )  ! Output, d(f)/d(RH)
    ! Arguments
    REAL(fp), INTENT(IN)  :: RH_LUT(:)
    REAL(fp), INTENT(IN)  :: fmass1_LUT(:)
    REAL(fp), INTENT(IN)  :: RH
    REAL(fp), INTENT(OUT) :: f
    REAL(fp), INTENT(OUT) :: f_RH
    ! Local variables
    INTEGER :: i, n

    n = SIZE(RH_LUT)
    IF ( RH <= RH_LUT(1) ) THEN
      f    = fmass1_LUT(1)
      f_RH = ZERO
    ELSE IF ( RH >= RH_LUT(n) ) THEN
      f    = fmass1_LUT(n)
      f_RH = ZERO
    ELSE
      DO i = 1, n-1
        IF ( RH <= RH_LUT(i+1) ) EXIT
      END DO
      f_RH = (fmass1_LUT(i+1) - fmass1_LUT(i)) / (RH_LUT(i+1) - RH_LUT(i))
      f    = fmass1_LUT(i) + f_RH * (RH - RH_LUT(i))
    END IF

  END SUBROUTINE Interpolate_fmass1


  ! -------------------------------------------------------------
  ! Forward fine-mode scaling of the aerosol concentrations
  !   c' = f.c
  ! -------------------------------------------------------------
  SUBROUTINE FineMode_Scale( &
    Atm   , &  ! Input
    f     , &  ! Input
    Atm_FM  )  ! Output
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    REAL(fp),                   INTENT(IN)     :: f(:,:)
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_FM
    ! Local variables
    INTEGER :: na, nl

    DO na = 1, Atm%n_Aerosols
      nl = Atm%Aerosol(na)%n_Layers
      Atm_FM%Aerosol(na)%Concentration(1:nl) = f(1:nl,na) * Atm%Aerosol(na)%Concentration(1:nl)
    END DO

  END SUBROUTINE FineMode_Scale


  ! -------------------------------------------------------------
  ! Tangent-linear fine-mode scaling of the aerosol concentrations
  !   c'_TL = f.c_TL + f_RH.c.RH_TL
  ! -------------------------------------------------------------
  SUBROUTINE FineMode_Scale_TL( &
    Atm      , &  ! FWD Input
    Atm_TL   , &  ! TL  Input
    f        , &  ! Input
    f_RH     , &  ! Input
    Atm_FM_TL  )  ! TL  Output
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm_TL
    REAL(fp),                   INTENT(IN)     :: f(:,:)
    REAL(fp),                   INTENT(IN)     :: f_RH(:,:)
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_FM_TL
    ! Local variables
    INTEGER :: na, nl

    DO na = 1, Atm%n_Aerosols
      nl = Atm%Aerosol(na)%n_Layers
      Atm_FM_TL%Aerosol(na)%Concentration(1:nl) = &
        f(1:nl,na) * Atm_TL%Aerosol(na)%Concentration(1:nl) + &
        f_RH(1:nl,na) * Atm%Aerosol(na)%Concentration(1:nl) * Atm_TL%Relative_Humidity(1:nl)
    END DO

  END SUBROUTINE FineMode_Scale_TL


  ! -------------------------------------------------------------
  ! Subroutine to save the input adjoint aerosol concentrations and
  ! set them to zero, so that the adjoint of the scaled
  ! concentrations, c'_AD, can be accumulated separately.
  ! -------------------------------------------------------------
  SUBROUTINE FineMode_Save_AD( &
    Atm_AD    , &  ! AD In/Output
    Conc_AD_in  )  ! Output, L x n_Aerosols
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_AD
    REAL(fp), ALLOCATABLE,      INTENT(OUT)    :: Conc_AD_in(:,:)
    ! Local variables
    INTEGER :: na, nl

    ALLOCATE( Conc_AD_in(Atm_AD%n_Layers, Atm_AD%n_Aerosols) )
    Conc_AD_in = ZERO
    DO na = 1, Atm_AD%n_Aerosols
      nl = Atm_AD%Aerosol(na)%n_Layers
      Conc_AD_in(1:nl,na) = Atm_AD%Aerosol(na)%Concentration(1:nl)
      Atm_AD%Aerosol(na)%Concentration(1:nl) = ZERO
    END DO

  END SUBROUTINE FineMode_Save_AD


  ! -------------------------------------------------------------
  ! Adjoint fine-mode scaling of the aerosol concentrations. On
  ! input the Atm_AD concentrations hold c'_AD, the adjoint of the
  ! scaled concentrations:
  !   RH_AD = RH_AD + f_RH.c.c'_AD
  !   c_AD  = c_AD(input) + f.c'_AD
  ! -------------------------------------------------------------
  SUBROUTINE FineMode_Scale_AD( &
    Atm       , &  ! FWD Input
    f         , &  ! Input
    f_RH      , &  ! Input
    Conc_AD_in, &  ! Input
    Atm_AD      )  ! AD In/Output
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    REAL(fp),                   INTENT(IN)     :: f(:,:)
    REAL(fp),                   INTENT(IN)     :: f_RH(:,:)
    REAL(fp),                   INTENT(IN)     :: Conc_AD_in(:,:)
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_AD
    ! Local variables
    INTEGER :: na, nl

    DO na = 1, Atm%n_Aerosols
      nl = Atm%Aerosol(na)%n_Layers
      Atm_AD%Relative_Humidity(1:nl) = Atm_AD%Relative_Humidity(1:nl) + &
        f_RH(1:nl,na) * Atm%Aerosol(na)%Concentration(1:nl) * Atm_AD%Aerosol(na)%Concentration(1:nl)
      Atm_AD%Aerosol(na)%Concentration(1:nl) = Conc_AD_in(1:nl,na) + &
        f(1:nl,na) * Atm_AD%Aerosol(na)%Concentration(1:nl)
    END DO

  END SUBROUTINE FineMode_Scale_AD

END MODULE CRTM_AOD_Module
