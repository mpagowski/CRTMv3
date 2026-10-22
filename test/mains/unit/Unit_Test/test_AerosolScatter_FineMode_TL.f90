!
! test_AerosolScatter_FineMode_TL
!
! Unit test for the fine-mode aerosol optical depth (AOD) and absorption
! aerosol optical depth (AAOD) computed by CRTM_AOD and CRTM_AOD_TL with
! FINE_MODE = .TRUE. (GOCART-GEOS5 table).
! The profile contains one aerosol of every type in the table.
!
! Verifies that
!   - the AerosolCoeff fmass1 values (fraction of mass < 1 um diameter)
!     of every aerosol type match the GMAO optics_DU/SS/NI/OC/BC/SU
!     v2.2.1.legacy tables at RH 0.80 and 0.99
!   - the fine-mode AOD/AAOD equals the total AOD/AAOD of a profile whose
!     aerosol concentrations are scaled by fmass1(RH)
!   - 0 <= fine-mode AOD/AAOD <= total AOD/AAOD
!   - the fine-mode TL agrees with a central finite difference of the
!     forward model for perturbations of concentration, relative humidity,
!     and both.
!
! CREATION HISTORY:
!       Written by:     Mariusz Pagowski, 2026/10/22
!                       mariusz.pagowski@colorado.edu
!

PROGRAM test_AerosolScatter_FineMode_TL

  USE CRTM_Module
  USE CRTM_AerosolCoeff, ONLY: AeroC

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_FineMode_TL'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  INTEGER,  PARAMETER :: N_LAYERS = 5
  ! Relative humidities away from the LUT RH nodes. At RH 0.935 the fmass1
  ! of sea salt bin 2, the nitrate bins, sulfate, OC2 and BC2 all depend on RH.
  REAL(fp), PARAMETER :: RH(N_LAYERS) = (/ 0.43_fp, 0.57_fp, 0.66_fp, 0.72_fp, 0.935_fp /)
  REAL(fp), PARAMETER :: FWD_RTOL = 1.0e-12_fp
  REAL(fp), PARAMETER :: FD_STEP  = 1.0e-4_fp
  REAL(fp), PARAMETER :: FD_RTOL  = 1.0e-5_fp
  REAL(fp), PARAMETER :: LUT_RTOL = 1.0e-12_fp
  ! GMAO fmass1 values at RH 0.80 and 0.99, in AerosolCoeff.GOCART-GEOS5
  ! type order: Dust 1-5, Sea salt 1-5, Organic carbon 1-2, Black carbon 1-2,
  ! Sulfate 1-2, Nitrate 1-3
  INTEGER,  PARAMETER :: N_TYPES = 19
  CHARACTER(16), PARAMETER :: TYPE_NAME(N_TYPES) = (/ CHARACTER(16) :: &
    'Dust 1', 'Dust 2', 'Dust 3', 'Dust 4', 'Dust 5', &
    'Sea salt 1', 'Sea salt 2', 'Sea salt 3', 'Sea salt 4', 'Sea salt 5', &
    'Organic carbon 1', 'Organic carbon 2', 'Black carbon 1', 'Black carbon 2', &
    'Sulfate 1', 'Sulfate 2', 'Nitrate 1', 'Nitrate 2', 'Nitrate 3' /)
  REAL(fp), PARAMETER :: FMASS1_RH80(N_TYPES) = (/ &
    0.2069999818951323_fp, ZERO, ZERO, ZERO, ZERO, &
    ONE, 0.46066838679702593_fp, ZERO, ZERO, ZERO, &
    ONE, ONE, ONE, ONE, &
    ONE, ONE, &
    0.6947407853078148_fp, 0.0007819185556116189_fp, 1.7056273925616346E-07_fp /)
  REAL(fp), PARAMETER :: FMASS1_RH99(N_TYPES) = (/ &
    0.2069999818951323_fp, ZERO, ZERO, ZERO, ZERO, &
    ONE, 0.20059365007168792_fp, ZERO, ZERO, ZERO, &
    ONE, 0.8102859980167078_fp, ONE, 0.9970127444337061_fp, &
    0.7031204303629337_fp, 0.7031204303629337_fp, &
    0.5478982345626923_fp, 0.00019148497150574668_fp, 1.936764165554449E-08_fp /)

  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  INTEGER :: Error_Status, n_Channels, n_Aerosols, l, k, n, id, ia, n_Failures
  LOGICAL :: Compute_AAOD
  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Options_type)     :: Options(1)
  TYPE(CRTM_Atmosphere_type)  :: Atm(1), Atm_FM(1), Atm_TL(1), Atm_p(1), Atm_m(1)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:), RTS_p(:,:), RTS_m(:,:)
  REAL(fp), ALLOCATABLE :: aod_fm(:,:), aod_tot(:,:), ref(:,:), aod_tl(:,:), aod_fd(:,:)
  REAL(fp) :: rel_err
  CHARACTER(4)  :: mname(2) = (/ 'AOD ', 'AAOD' /)
  CHARACTER(20) :: dname(3) = (/ 'concentration       ', 'relative humidity   ', 'concentration + RH  ' /)

  PRINT *, 'Starting test_AerosolScatter_FineMode_TL...'
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
            aod_fm(N_LAYERS,n_Channels), aod_tot(N_LAYERS,n_Channels), ref(N_LAYERS,n_Channels), &
            aod_tl(N_LAYERS,n_Channels), aod_fd(N_LAYERS,n_Channels) )
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
  ! 1. fmass1 values from the GMAO tables
  ! ------------------------------------------
  CALL Check_fmass1()


  Mode_Loop: DO ia = 1, 2
    Compute_AAOD = ( ia == 2 )
    PRINT '(/1x,"Fine-mode ",a)', TRIM(mname(ia))

    ! ------------------------------------------
    ! 2. Forward: fine-mode = total for fmass1-scaled concentrations
    ! ------------------------------------------
    Error_Status = CRTM_AOD( Atm, ChannelInfo, RTSolution, Options=Options, &
                             AAOD=Compute_AAOD, FINE_MODE=.TRUE. )
    IF ( Error_Status /= SUCCESS ) STOP 1
    DO l = 1, n_Channels
      aod_fm(:,l) = RTSolution(l,1)%Layer_Optical_Depth
    END DO

    Error_Status = CRTM_AOD( Atm, ChannelInfo, RTSolution, Options=Options, AAOD=Compute_AAOD )
    IF ( Error_Status /= SUCCESS ) STOP 1
    DO l = 1, n_Channels
      aod_tot(:,l) = RTSolution(l,1)%Layer_Optical_Depth
    END DO

    Atm_FM = Atm
    DO n = 1, n_Aerosols
      DO k = 1, N_LAYERS
        Atm_FM(1)%Aerosol(n)%Concentration(k) = Reference_fmass1( n, RH(k) ) * Atm(1)%Aerosol(n)%Concentration(k)
      END DO
    END DO
    Error_Status = CRTM_AOD( Atm_FM, ChannelInfo, RTSolution, Options=Options, AAOD=Compute_AAOD )
    IF ( Error_Status /= SUCCESS ) STOP 1
    DO l = 1, n_Channels
      ref(:,l) = RTSolution(l,1)%Layer_Optical_Depth
    END DO

    PRINT '(1x,"Channel 1 column ",a,": fine-mode ",es16.8,", reference ",es16.8,", total ",es16.8)', &
          TRIM(mname(ia)), SUM(aod_fm(:,1)), SUM(ref(:,1)), SUM(aod_tot(:,1))
    CALL Check( 'Fine-mode '//TRIM(mname(ia))//' is non-zero', MAXVAL(aod_fm) > ZERO )
    rel_err = MAXVAL(ABS(aod_fm - ref)) / MAXVAL(ABS(ref))
    PRINT '(1x,"Rel. Diff: ",es12.4)', rel_err
    CALL Check( 'Fine-mode '//TRIM(mname(ia))//' = total for fmass1-scaled concentrations', &
                rel_err <= FWD_RTOL )
    CALL Check( '0 <= fine-mode '//TRIM(mname(ia))//' <= total', &
                ALL(aod_fm >= ZERO) .AND. ALL(aod_fm <= aod_tot*(ONE+FWD_RTOL)) )
    CALL Check( 'Fine-mode '//TRIM(mname(ia))//' < total', ALL(aod_fm < aod_tot) )


    ! ------------------------------------------
    ! 3. TL vs central finite difference
    ! ------------------------------------------
    DO id = 1, 3
      CALL Set_Perturbation( id )
      Error_Status = CRTM_AOD_TL( Atm, Atm_TL, ChannelInfo, RTSolution, RTSolution_TL, &
                                  Options=Options, AAOD=Compute_AAOD, FINE_MODE=.TRUE. )
      IF ( Error_Status /= SUCCESS ) STOP 1
      DO l = 1, n_Channels
        aod_tl(:,l) = RTSolution_TL(l,1)%Layer_Optical_Depth
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
      Error_Status = CRTM_AOD( Atm_p, ChannelInfo, RTS_p, Options=Options, &
                               AAOD=Compute_AAOD, FINE_MODE=.TRUE. )
      IF ( Error_Status /= SUCCESS ) STOP 1
      Error_Status = CRTM_AOD( Atm_m, ChannelInfo, RTS_m, Options=Options, &
                               AAOD=Compute_AAOD, FINE_MODE=.TRUE. )
      IF ( Error_Status /= SUCCESS ) STOP 1
      DO l = 1, n_Channels
        aod_fd(:,l) = (RTS_p(l,1)%Layer_Optical_Depth - RTS_m(l,1)%Layer_Optical_Depth) / (2.0_fp*FD_STEP)
      END DO

      PRINT '(1x,"TL vs FD, ",a)', TRIM(dname(id))
      PRINT '(3x,"max |",a," TL| = ",es12.4)', TRIM(mname(ia)), MAXVAL(ABS(aod_tl))
      IF ( MAXVAL(ABS(aod_tl)) > ZERO ) THEN
        rel_err = MAXVAL(ABS(aod_tl - aod_fd)) / MAXVAL(ABS(aod_tl))
        PRINT '(3x,"Rel. Diff: ",es12.4)', rel_err
        CALL Check( 'Fine-mode '//TRIM(mname(ia))//' TL vs FD, '//TRIM(dname(id)), rel_err <= FD_RTOL )
      ELSE
        CALL Check( 'Fine-mode '//TRIM(mname(ia))//' TL vs FD (TL is ZERO), '//TRIM(dname(id)), .FALSE. )
      END IF
    END DO
  END DO Mode_Loop


  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_FM )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_p )
  CALL CRTM_Atmosphere_Destroy( Atm_m )
  CALL CRTM_RTSolution_Destroy( RTSolution )
  CALL CRTM_RTSolution_Destroy( RTSolution_TL )
  CALL CRTM_RTSolution_Destroy( RTS_p )
  CALL CRTM_RTSolution_Destroy( RTS_m )
  DEALLOCATE( RTSolution, RTSolution_TL, RTS_p, RTS_m, aod_fm, aod_tot, ref, aod_tl, aod_fd )
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( n_Failures > 0 ) THEN
    PRINT '(1x,"FAIL: ",i0," fine-mode check(s) failed.")', n_Failures
    STOP 1
  END IF
  PRINT *, 'SUCCESS: All fine-mode forward and TL checks passed.'

CONTAINS

  ! Check the AerosolCoeff fmass1 of every aerosol type, matched by name,
  ! against the GMAO values at RH 0.80 and 0.99
  SUBROUTINE Check_fmass1()
    INTEGER :: n, nt, i, ih80, ih99, n_found
    CHARACTER(80) :: name
    LOGICAL :: ok
    ok = .TRUE.
    n_found = 0
    ih80 = MINLOC(ABS(AeroC%RH - 0.80_fp), DIM=1)
    ih99 = MINLOC(ABS(AeroC%RH - 0.99_fp), DIM=1)
    DO n = 1, n_Aerosols
      ! Type names in the coefficient files are NUL padded
      name = AeroC%Type_Name(n)
      DO i = 1, LEN(name)
        IF ( IACHAR(name(i:i)) < 32 ) name(i:i) = ' '
      END DO
      DO nt = 1, N_TYPES
        IF ( TRIM(ADJUSTL(name)) /= TRIM(TYPE_NAME(nt)) ) CYCLE
        n_found = n_found + 1
        IF ( ABS(AeroC%fmass1(ih80,n) - FMASS1_RH80(nt)) > LUT_RTOL*ABS(FMASS1_RH80(nt)) .OR. &
             ABS(AeroC%fmass1(ih99,n) - FMASS1_RH99(nt)) > LUT_RTOL*ABS(FMASS1_RH99(nt)) ) THEN
          PRINT '(1x,a,": fmass1(RH 0.80, 0.99) = ",2es24.16)', TRIM(TYPE_NAME(nt)), &
                AeroC%fmass1(ih80,n), AeroC%fmass1(ih99,n)
          ok = .FALSE.
        END IF
      END DO
    END DO
    CALL Check( 'All 19 aerosol types found by name', n_found == N_TYPES )
    CALL Check( 'AerosolCoeff fmass1 matches the GMAO tables', ok )
  END SUBROUTINE Check_fmass1


  ! Reference fmass1 for aerosol n, linearly interpolated in RH
  FUNCTION Reference_fmass1( n, h ) RESULT( f )
    INTEGER,  INTENT(IN) :: n
    REAL(fp), INTENT(IN) :: h
    REAL(fp) :: f
    INTEGER :: i
    DO i = 1, AeroC%n_RH-1
      IF ( h <= AeroC%RH(i+1) ) EXIT
    END DO
    f = AeroC%fmass1(i,n) + (AeroC%fmass1(i+1,n) - AeroC%fmass1(i,n)) * &
                            (h - AeroC%RH(i)) / (AeroC%RH(i+1) - AeroC%RH(i))
  END FUNCTION Reference_fmass1


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

END PROGRAM test_AerosolScatter_FineMode_TL
